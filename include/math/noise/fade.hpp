/**
 * @file fade.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef SBODY_FADE_HPP
#define SBODY_FADE_HPP

#include <stdexcept>

#include <boost/numeric/ublas/vector.hpp>

namespace sbody {
	template <typename T, std::size_t N>
	T Smoothstep(T x, T x0, T x1, size_t continuity_order) {
		if (x0 == x1) {
			if (x != x0)
				throw std::invalid_argument("Smoothstep: x is out of range.");
			return T{0.5};
		}
		T t = (x - x0) / (x1 - x0);
		if (t < 0 || t > 1)
			throw std::invalid_argument("Smoothstep: x is out of range.");
		if (continuity_order == 0)
			return t;
		else if (continuity_order == 1)
			return (3. - 2. * t) * power2(t);
		else if (continuity_order == 2)
			return ((6. * t - 15.) * t + 10.) * power3(t);
		else if (continuity_order == 3)
			return -(((20. * t - 70.) * t + 84. * t) - 35.) * power4(t);
		throw std::out_of_range("Smoothstep: continuity_order is not supportted.");
	}
} // namespace sbody

#endif
