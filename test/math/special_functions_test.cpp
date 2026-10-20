/**
 * @file special_functions_test.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2024-02-21
 *
 * @copyright Copyright (c) 2024
 *
 */

#define BOOST_TEST_MODULE special_functions_test

#include <array>
#include <chrono>
#include <csignal>
#include <memory>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/test/included/unit_test.hpp>

#include "math/special_functions.hpp"

using namespace std;
using namespace sbody;

BOOST_AUTO_TEST_CASE(TestCarlsonRC) {
	BOOST_CHECK_CLOSE_FRACTION(carlsonRC(static_cast<long double>(0.1), static_cast<long double>(-0.2)), 1.2022125790504758, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(carlsonRC(static_cast<long double>(1.0), static_cast<long double>(2.0)), 0.7853981633974483, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(carlsonRC(static_cast<long double>(0.5), static_cast<long double>(1.0)), 1.1107207345395915, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestCarlsonRJ) {
	BOOST_CHECK_CLOSE_FRACTION(carlsonRJ(0.2l, 0.1l, 0.3l, 0.4l), 7.584662876718954, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(carlsonRJ(1.0l, 2.0l, 3.0l, 4.0l), 0.2398480997495678l, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(carlsonRJ(0.5l, 0.5l, 1.0l, 0.5l), 2.260648559158617l, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestEllint_1) {
	const long double phi = 0.5l;
	const long double k = 0.2l;
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(k, phi), 0.5007959922424561, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(k), 1.5868678474541662, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(0.5l, 1.0l), 1.0373561200021773l, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(0.0l), boost::math::constants::pi<long double>() / 2., 1e-15);
}

BOOST_AUTO_TEST_CASE(TestEllint_2) {
	const long double phi = 0.5l;
	const long double k = 0.2l;
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_2(k, phi), 0.49920624167354444, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_2(k), 1.5549685462425293, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_2(0.5l, 1.0l), 0.9648764542686275l, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_2(0.0l), boost::math::constants::pi<long double>() / 2., 1e-15);
}

BOOST_AUTO_TEST_CASE(TestEllipticIntegral) {
	// Integral of (0.5+u)^(p5/2) / sqrt((1+u)(2+u)(3+u)(4+u)) from 0.1 to 0.2;
	// the reference values come from mpmath.quad.
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral(0, 0.1l, 0.2l, 0.5l, 1.0l, 1.0l, 1.0l, 2.0l, 1.0l, 3.0l, 1.0l, 4.0l, 1.0l), 0.017600181810740477l, 1e-12);
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral(-2, 0.1l, 0.2l, 0.5l, 1.0l, 1.0l, 1.0l, 2.0l, 1.0l, 3.0l, 1.0l, 4.0l, 1.0l), 0.027163810536225407l, 1e-12);
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral(-4, 0.1l, 0.2l, 0.5l, 1.0l, 1.0l, 1.0l, 2.0l, 1.0l, 3.0l, 1.0l, 4.0l, 1.0l), 0.042007168764765724l, 1e-12);
}

BOOST_AUTO_TEST_CASE(TestEllipticIntegral2Complex) {
	// xi^2 = (2+u)(3+u), X1^2 = 1+u, X4^2 = 4+u: the same integral as TestEllipticIntegral.
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral2Complex(0, 0.1l, 0.2l, 0.5l, 1.0l, 6.0l, 5.0l, 1.0l, 1.0l, 1.0l, 4.0l, 1.0l), 0.017600181810740477l, 1e-12);
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral2Complex(-2, 0.1l, 0.2l, 0.5l, 1.0l, 6.0l, 5.0l, 1.0l, 1.0l, 1.0l, 4.0l, 1.0l), 0.027163810536225407l, 1e-12);
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral2Complex(-4, 0.1l, 0.2l, 0.5l, 1.0l, 6.0l, 5.0l, 1.0l, 1.0l, 1.0l, 4.0l, 1.0l), 0.042007168764765724l, 1e-12);
}

BOOST_AUTO_TEST_CASE(TestEllipticIntegral4Complex) {
	// xi1^2 = (1+u)(2+u), xi2^2 = (3+u)(4+u): the same integral as TestEllipticIntegral.
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral4Complex(0, 0.1l, 0.2l, 0.5l, 1.0l, 2.0l, 3.0l, 1.0l, 12.0l, 7.0l, 1.0l), 0.017600181810740477l, 1e-12);
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral4Complex(-2, 0.1l, 0.2l, 0.5l, 1.0l, 2.0l, 3.0l, 1.0l, 12.0l, 7.0l, 1.0l), 0.027163810536225407l, 1e-12);
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral4Complex(-4, 0.1l, 0.2l, 0.5l, 1.0l, 2.0l, 3.0l, 1.0l, 12.0l, 7.0l, 1.0l), 0.042007168764765724l, 1e-12);
}

BOOST_AUTO_TEST_CASE(TestEllipticIntegralGeokerrPattern) {
	// The real call pattern: a5 = 0, b5 = 1 and the four roots alternate in sign
	// (u1 = -2, u2 = 5, u3 = 10, u4 = 20).  The exact integral is
	// int_0.1^0.2 du / (u sqrt((2+u)(5-u)(10-u)(20-u))) = 0.1108085561343723.
	// All three routines must agree with that value for p5 = -4.
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral(-4, 0.1l, 0.2l, 0.0l, 1.0l, 2.0l, 1.0l, 5.0l, -1.0l, 10.0l, -1.0l, 20.0l, -1.0l), 0.1108085561343723l, 1e-12);
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral2Complex(-4, 0.1l, 0.2l, 0.0l, 1.0l, 50.0l, -15.0l, 1.0l, 2.0l, 1.0l, 20.0l, -1.0l), 0.1108085561343723l, 1e-12);
	BOOST_CHECK_CLOSE_FRACTION(ellipticIntegral4Complex(-4, 0.1l, 0.2l, 0.0l, 1.0l, 10.0l, 3.0l, -1.0l, 200.0l, -30.0l, 1.0l), 0.1108085561343723l, 1e-12);
}

BOOST_AUTO_TEST_CASE(TestJacobiElliptic) {
	const double k = 0.2l;
	const double u = 0.5l;
	double cn, dn;
	double sn = boost::math::jacobi_elliptic(k, u, &cn, &dn);
	BOOST_CHECK_CLOSE_FRACTION(sn, 0.47873005322805676, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(cn, 0.8779621496034222, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(dn, 0.9954057973738401, 1e-15);
}
