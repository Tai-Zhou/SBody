/**
 * @file geometry_test.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2024-02-21
 *
 * @copyright Copyright (c) 2024
 *
 */

#define BOOST_TEST_MODULE geometry_test

#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/test/included/unit_test.hpp>

#include "math/geometry.hpp"

using namespace std;
using namespace sbody;

#define PRECISION long double

BOOST_AUTO_TEST_CASE(TestRotationAngleOfPerpendicularCoordiante) {
	double sin_theta_observer = M_SQRT1_2;
	double cos_theta_observer = M_SQRT1_2;
	double iota, sin_iota, cos_iota, phi;
	rotationAngleOfPerpendicularCoordiante(sin_theta_observer, cos_theta_observer, 1., 0., M_SQRT1_2, -M_SQRT1_2, iota, sin_iota, cos_iota, phi);
	BOOST_CHECK_CLOSE_FRACTION(iota, std::atan2(0.5, M_SQRT1_2), 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(phi, boost::math::constants::two_thirds_pi<double>(), 1e-15);
}
