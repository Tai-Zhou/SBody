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

BOOST_AUTO_TEST_CASE(TestRotationAngleOfPerpendicularCoordiante) {
	double sin_theta_observer = M_SQRT1_2;
	double cos_theta_observer = M_SQRT1_2;
	double iota, sin_iota, cos_iota, phi;
	rotationAngleOfPerpendicularCoordiante(sin_theta_observer, cos_theta_observer, 1., 0., M_SQRT1_2, -M_SQRT1_2, iota, sin_iota, cos_iota, phi);
	BOOST_CHECK_CLOSE_FRACTION(iota, std::atan2(0.5, M_SQRT1_2), 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(phi, boost::math::constants::two_thirds_pi<double>(), 1e-15);
}

BOOST_AUTO_TEST_CASE(TestSinTheta) {
	BOOST_CHECK_CLOSE_FRACTION(sinTheta(static_cast<long double>(0.)), 0., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(sinTheta(static_cast<long double>(M_PI / 4.)), M_SQRT1_2, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(sinTheta(static_cast<long double>(M_PI / 2.)), 1., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(sinTheta(static_cast<long double>(3. * M_PI / 4.)), M_SQRT1_2, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(sinTheta(static_cast<long double>(M_PI)), 0., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(sinTheta(static_cast<long double>(-M_PI / 4.)), M_SQRT1_2, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestCosTheta) {
	BOOST_CHECK_CLOSE_FRACTION(cosTheta(static_cast<long double>(0.)), 1., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(cosTheta(static_cast<long double>(M_PI / 4.)), M_SQRT1_2, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(cosTheta(static_cast<long double>(M_PI / 2.)), 0., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(cosTheta(static_cast<long double>(3. * M_PI / 4.)), -M_SQRT1_2, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(cosTheta(static_cast<long double>(M_PI)), -1., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(cosTheta(static_cast<long double>(-M_PI / 4.)), M_SQRT1_2, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestDotCross) {
	namespace ublas = boost::numeric::ublas;
	ublas::bounded_vector<long double, 3> x, y, z;
	x(0) = 1.;
	x(1) = 0.;
	x(2) = 0.;
	y(0) = 0.;
	y(1) = 1.;
	y(2) = 0.;
	z(0) = 0.;
	z(1) = 0.;
	z(2) = 1.;
	BOOST_CHECK_CLOSE_FRACTION(dotCross(x, y, z), 1., 1e-15);

	x(0) = 2.;
	x(1) = 3.;
	x(2) = 4.;
	y(0) = 5.;
	y(1) = 6.;
	y(2) = 7.;
	z(0) = 8.;
	z(1) = 9.;
	z(2) = 10.;
	BOOST_CHECK_CLOSE_FRACTION(dotCross(x, y, z), 27., 1e-15);
}

BOOST_AUTO_TEST_CASE(TestTriangleArea) {
	namespace ublas = boost::numeric::ublas;
	ublas::bounded_vector<long double, 3> x, y, z;
	x(0) = 0.;
	x(1) = 0.;
	x(2) = 0.;
	y(0) = 1.;
	y(1) = 0.;
	y(2) = 0.;
	z(0) = 0.;
	z(1) = 1.;
	z(2) = 0.;
	BOOST_CHECK_CLOSE_FRACTION(triangleArea(x, y, z), 0.5, 1e-15);

	x(0) = 0.;
	x(1) = 0.;
	x(2) = 0.;
	y(0) = 3.;
	y(1) = 0.;
	y(2) = 0.;
	z(0) = 0.;
	z(1) = 4.;
	z(2) = 0.;
	BOOST_CHECK_CLOSE_FRACTION(triangleArea(x, y, z), 6., 1e-15);
}

BOOST_AUTO_TEST_CASE(TestPointInTriangle) {
	BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., 0., 0.), true);
	BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., 0.5, 0.5), true);
	BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., 0.25, 0.25), true);
	BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., 1., 1.), false);
	BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., -0.1, -0.1), false);
	BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., 1.1, 0.), false);
}

BOOST_AUTO_TEST_CASE(TestModBy2Pi) {
	BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(static_cast<long double>(0.)), 0., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(static_cast<long double>(M_PI)), M_PI, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(static_cast<long double>(2. * M_PI)), 0., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(static_cast<long double>(3. * M_PI)), M_PI, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(static_cast<long double>(-M_PI)), M_PI, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(static_cast<long double>(5. * M_PI)), M_PI, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(static_cast<long double>(0.5 * M_PI)), 0.5 * M_PI, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestPhiDifference) {
	BOOST_CHECK_CLOSE_FRACTION(phiDifference(static_cast<long double>(0.)), 0., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(phiDifference(static_cast<long double>(M_PI)), M_PI, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(phiDifference(static_cast<long double>(2. * M_PI)), 0., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(phiDifference(static_cast<long double>(3. * M_PI)), M_PI, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(phiDifference(static_cast<long double>(-M_PI)), -M_PI, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(phiDifference(static_cast<long double>(0.5 * M_PI)), 0.5 * M_PI, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(phiDifference(static_cast<long double>(1.5 * M_PI)), -0.5 * M_PI, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(phiDifference(static_cast<long double>(-0.5 * M_PI)), -0.5 * M_PI, 1e-15);
}
