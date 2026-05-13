/**
 * @file polynomial_test.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-03
 *
 * @copyright Copyright (c) 2026
 *
 */

#define BOOST_TEST_MODULE polynomial_test

#include <boost/test/included/unit_test.hpp>

#include "math/solver/polynomial.hpp"

using namespace std;
using namespace sbody;

BOOST_AUTO_TEST_CASE(PolySolveQuadraticTest) {
	long double roots[2];
	BOOST_CHECK_EQUAL(polySolveQuadratic(static_cast<long double>(1.), static_cast<long double>(1.), static_cast<long double>(-1.), roots), 2);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -1.6180339887498948482045868343656381177203091798057638l, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], 0.61803398874989484820458683436563811772030917980576384l, 1e-15);
}

BOOST_AUTO_TEST_CASE(PolySolveQuadraticZeroATest) {
	long double roots[2];
	BOOST_CHECK_EQUAL(polySolveQuadratic(static_cast<long double>(0.), static_cast<long double>(2.), static_cast<long double>(-4.), roots), 1);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], 2., 1e-15);
}

BOOST_AUTO_TEST_CASE(PolySolveQuadraticZeroBTest) {
	long double roots[2];
	BOOST_CHECK_EQUAL(polySolveQuadratic(static_cast<long double>(1.), static_cast<long double>(0.), static_cast<long double>(-4.), roots), 2);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -2., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], 2., 1e-15);
}

BOOST_AUTO_TEST_CASE(PolySolveQuadraticNoRealRootsTest) {
	long double roots[2];
	BOOST_CHECK_EQUAL(polySolveQuadratic(static_cast<long double>(1.), static_cast<long double>(0.), static_cast<long double>(1.), roots), 0);
}

BOOST_AUTO_TEST_CASE(PolySolveQuadraticDoubleRootTest) {
	long double roots[2];
	BOOST_CHECK_EQUAL(polySolveQuadratic(static_cast<long double>(1.), static_cast<long double>(-2.), static_cast<long double>(1.), roots), 2);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], 1., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], 1., 1e-15);
}

BOOST_AUTO_TEST_CASE(PolySolveCubicTest) {
	long double roots[3];
	BOOST_CHECK_EQUAL(polySolveCubic(static_cast<long double>(5.), static_cast<long double>(1.), static_cast<long double>(-1.), roots), 3);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -4.74482607768192328566, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], -0.604068139818793694131, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[2], 0.348894217500716979788, 1e-15);
}

BOOST_AUTO_TEST_CASE(PolySolveCubicZeroCTest) {
	long double roots[3];
	BOOST_CHECK_EQUAL(polySolveCubic(static_cast<long double>(1.), static_cast<long double>(-6.), static_cast<long double>(0.), roots), 3);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -3., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], 0., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[2], 2., 1e-15);
}

BOOST_AUTO_TEST_CASE(PolySolveCubicOneRealRootTest) {
	long double roots[3];
	BOOST_CHECK_EQUAL(polySolveCubic(static_cast<long double>(0.), static_cast<long double>(1.), static_cast<long double>(1.), roots), 1);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -1., 1e-15);
}

BOOST_AUTO_TEST_CASE(PolySolveQuarticTest) {
	long double roots[4];
	BOOST_CHECK_EQUAL(polySolveQuartic(static_cast<long double>(0.), static_cast<long double>(-10.), static_cast<long double>(35.), static_cast<long double>(-50.), roots), 4);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], 2., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], 2., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[2], 2., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[3], 5., 1e-15);
}

BOOST_AUTO_TEST_CASE(PolySolveQuarticZeroDTest) {
	long double roots[4];
	BOOST_CHECK_EQUAL(polySolveQuartic(static_cast<long double>(0.), static_cast<long double>(-5.), static_cast<long double>(4.), static_cast<long double>(0.), roots), 4);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -2., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], 0., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[2], 1., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[3], 4., 1e-15);
}

BOOST_AUTO_TEST_CASE(PolySolveQuarticNoRealRootsTest) {
	long double roots[4];
	BOOST_CHECK_EQUAL(polySolveQuartic(static_cast<long double>(0.), static_cast<long double>(0.), static_cast<long double>(1.), static_cast<long double>(1.), roots), 0);
}

BOOST_AUTO_TEST_CASE(PolySolveQuarticBiquadraticTest) {
	long double roots[4];
	BOOST_CHECK_EQUAL(polySolveQuartic(static_cast<long double>(0.), static_cast<long double>(-5.), static_cast<long double>(4.), static_cast<long double>(0.), roots), 4);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -2., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], 0., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[2], 1., 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[3], 4., 1e-15);
}
