/**
 * @file disk.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef SBODY_DISK_HPP
#define SBODY_DISK_HPP

#include "object/object.hpp"

namespace sbody {
	template <typename T>
	class Disk : public Object<T> {
	  protected:
		const T inner_radius_;
		const T outer_radius_;

	  public:
		Disk(const Metric<T> &metric, T inner_radius, T outer_radius) : Object(metric), inner_radius_(inner_radius), outer_radius_(outer_radius) {}
		T angularVelocity(T radius) const = 0;
		T temperature(T radius) const = 0;
		T intensity(T radius) const = 0;
	};
} // namespace sbody

#endif
