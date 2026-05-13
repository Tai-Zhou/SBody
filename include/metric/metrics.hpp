/**
 * @file metrics.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-16
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <string>

#include "metric/kerr.hpp"
#include "metric/kerr_newman.hpp"
#include "metric/kerr_taub_nut.hpp"
#include "metric/metric.hpp"
#include "metric/post_newtonian.hpp"
#include "metric/reissner_nordstrom.hpp"
#include "metric/schwarzschild.hpp"

namespace sbody {
enum class MetricType {
    POST_NEWTONIAN,
    SCHWARZSCHILD,
    REISSNER_NORDSTROM,
    KERR,
    KERR_NEWMAN,
    KERR_TAUB_NUT,
};
inline MetricType getMetricType(const std::string &name) {
    std::string lower_name;
    std::remove_copy(name.begin(), name.end(), std::back_inserter(lower_name), '-');
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);
    if (lower_name == "postnewtonian")
        return MetricType::POST_NEWTONIAN;
    else if (lower_name == "schwarzschild")
        return MetricType::SCHWARZSCHILD;
    else if (lower_name == "reissnernordstrom")
        return MetricType::REISSNER_NORDSTROM;
    else if (lower_name == "kerr")
        return MetricType::KERR;
    else if (lower_name == "kerrnewman")
        return MetricType::KERR_NEWMAN;
    else if (lower_name == "kerrtaubnut")
        return MetricType::KERR_TAUB_NUT;
    // else if (lower_name == "hayward")
    //     return MetricType::Hayward;
    throw std::invalid_argument("Invalid metric name: " + name);
}
template <typename T>
struct MetricParameters {
    int PN_ = 1;
    T PN1_ = 1.;
    T charge_ = 0.;
    T spin_ = 0.;
    T NUT_ = 0.;
};
template <typename T>
std::unique_ptr<Metric<T>> createMetric(const MetricType type, const MetricParameters<T> &params) {
    switch (type) {
    case MetricType::POST_NEWTONIAN:
        if (params.PN_ == 1)
            return std::make_unique<PN1<T>>(params.PN1_);
        return std::make_unique<PostNewtonian<T>>(params.PN_);
    case MetricType::SCHWARZSCHILD:
        return std::make_unique<Schwarzschild<T>>();
    case MetricType::REISSNER_NORDSTROM:
        return std::make_unique<ReissnerNordstrom<T>>(params.charge_);
    case MetricType::KERR:
        return std::make_unique<Kerr<T>>(params.spin_);
    case MetricType::KERR_NEWMAN:
        return std::make_unique<KerrNewman<T>>(params.spin_, params.charge_);
    case MetricType::KERR_TAUB_NUT:
        return std::make_unique<KerrTaubNUT<T>>(params.spin_, params.NUT_);
    default:
        throw std::invalid_argument("Invalid metric type");
    }
}
} // namespace sbody
