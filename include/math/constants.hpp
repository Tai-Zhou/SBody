/**
 * @file constant.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef SBODY_CONSTANT_HPP
#define SBODY_CONSTANT_HPP

#include <boost/math/constants/constants.hpp>

namespace sbody {
	/// Global absolute accuracy.
	inline double absolute_accuracy = 1e-12;

	/// Global relative accuracy.
	inline double relative_accuracy = 1e-12;

	/// \f$\sqrt{27}\f$
	constexpr long double M_SQRT27 = 5.19615242270663188058233902451761710l;

	/// Epsilon, a small value. \f$\varepsilon\f$.
	constexpr long double EPSILON = 1e-10l;

	/// \f$\sin\varepsilon\f$.
	constexpr long double SIN_EPSILON = 1e-10l;

	/// \f$\cos\varepsilon\f$.
	constexpr long double COS_EPSILON = 0.999999999999999999995l;

	/// Area of a circle with radius of SQRT_DBL_EPSILON. \f$\pi\varepsilon^2\f$
	constexpr long double EPSILON_CIRCLE_AREA = boost::math::long_double_constants::pi * EPSILON * EPSILON;

	/// Sample number.
	constexpr int SAMPLE_NUMBER = 100;

	/// The angle corresponding to the sample number.
	constexpr long double ANGLE_INTERVAL = boost::math::long_double_constants::two_pi / SAMPLE_NUMBER;

	/// Area of the regular polygon with SAMPLE_NUMBER edges.
	constexpr long double EPSILON_POLYGON_AREA = 0.06279051952931337 / ANGLE_INTERVAL * EPSILON_CIRCLE_AREA;
	// For C++20
	// constexpr long double EPSILON_POLYGON_AREA = std::sin(boost::math::long_double_constants::two_pi / SAMPLE_NUMBER) / ANGLE_INTERVAL * EPSILON_CIRCLE_AREA;
} // namespace sbody
#endif
