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
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRC(static_cast<long double>(0.1), static_cast<long double>(-0.2)), 1.2022125790504758, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRC(static_cast<long double>(1.0), static_cast<long double>(2.0)), 0.7853981633974483, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRC(static_cast<long double>(0.5), static_cast<long double>(1.0)), 1.1107207345395915, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestCarlsonRJ) {
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRJ(0.2l, 0.1l, 0.3l, 0.4l), 7.584662876718954, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRJ(1.0l, 2.0l, 3.0l, 4.0l), 0.7767570595569398, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRJ(0.5l, 0.5l, 1.0l, 0.5l), 3.141592653589793, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestEllint_1) {
	const long double phi = 0.5l;
	const long double k = 0.2l;
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(k, phi), 0.5007959922424561, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(k), 1.5868678474541662, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(0.5l, 1.0l), 0.881373587019543, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(0.0l), boost::math::constants::pi<long double>() / 2., 1e-15);
}

BOOST_AUTO_TEST_CASE(TestEllint_2) {
	const long double phi = 0.5l;
	const long double k = 0.2l;
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_2(k, phi), 0.49920624167354444, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_2(k), 1.5549685462425293, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_2(0.5l, 1.0l), 1.0, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_2(0.0l), boost::math::constants::pi<long double>() / 2., 1e-15);
}

BOOST_AUTO_TEST_CASE(TestEllipticIntegral) {
	BOOST_CHECK_CLOSE_FRACTION(EllipticIntegral(0, 0.1, 0.2, 0.0, 0.0, 0.2, -1.0, -0.1, 1.0, 0.3, -1.0), 0.1, 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(EllipticIntegral(-2, 0.1, 0.2, 0.0, 1.0, 0.2, -1.0, -0.1, 1.0, 0.3, -1.0), 0.0909090909090909, 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(EllipticIntegral(-4, 0.1, 0.2, 0.0, 1.0, 0.2, -1.0, -0.1, 1.0, 0.3, -1.0), 0.0826446280991736, 1e-10);
}

BOOST_AUTO_TEST_CASE(TestEllipticIntegral2Complex) {
	BOOST_CHECK_CLOSE_FRACTION(EllipticIntegral2Complex(0, 0.1, 0.2, 0.0, 1.0, -0.5, 0.1, 1.0, -0.1, 1.0), 0.1, 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(EllipticIntegral2Complex(-2, 0.1, 0.2, 0.0, 1.0, -0.5, 0.1, 1.0, -0.1, 1.0), 0.0909090909090909, 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(EllipticIntegral2Complex(-4, 0.1, 0.2, 0.0, 1.0, -0.5, 0.1, 1.0, -0.1, 1.0), 0.0826446280991736, 1e-10);
}

BOOST_AUTO_TEST_CASE(TestEllipticIntegral4Complex) {
	BOOST_CHECK_CLOSE_FRACTION(EllipticIntegral4Complex(0, 0.1, 0.2, 0.0, 1.0, -0.5, 0.1, 0.5, -0.5, 0.1), 0.1, 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(EllipticIntegral4Complex(-2, 0.1, 0.2, 0.0, 1.0, -0.5, 0.1, 0.5, -0.5, 0.1), 0.0909090909090909, 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(EllipticIntegral4Complex(-4, 0.1, 0.2, 0.0, 1.0, -0.5, 0.1, 0.5, -0.5, 0.1), 0.0826446280991736, 1e-10);
}

BOOST_AUTO_TEST_CASE(TestJacobiElliptic) {
	const double k = 0.2l;
	const double u = 0.5l;
	double cn, dn;
	double sn = boost::math::jacobi_elliptic(k, u, &cn, &dn);
	BOOST_CHECK_CLOSE_FRACTION(sn, 0.479425538604203, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(cn, 0.8775825618903728, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(dn, 0.9800665778412416, 1e-15);
}
