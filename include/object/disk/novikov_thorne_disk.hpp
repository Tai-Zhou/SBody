/**
 * @file novikov_thorne_disk.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-06
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include "object/disk/disk.hpp"

namespace sbody {
template <typename T>
class NovikovThorneDisk : public Disk<T> {
  protected:
  public:
    NovikovThorneDisk(const Metric<T> &metric, T outer_radius)
        : Disk(metric, metric.innermostCircularOrbitRadius().first, outer_radius) {}
    int Hit(const boost::numeric::ublas::bounded_vector<T, 8> &current,
            const boost::numeric::ublas::bounded_vector<T, 8> &last) {
        if (!OppositeSign(current[2], last[2]))
            return StatusCode::SUCCESS;
        if (std::abs(current[2]) < boost::math::constants::third_pi<T>())
            return StatusCode::SUCCESS;
        return StatusCode::FAILURE;
    }
};
} // namespace sbody
