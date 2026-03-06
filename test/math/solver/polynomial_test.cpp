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

#define PRECISION long double

BOOST_AUTO_TEST_CASE(PolySolveQuadraticTest) {
	PRECISION roots[2];
	BOOST_CHECK_EQUAL(polySolveQuadratic(static_cast<PRECISION>(1.), static_cast<PRECISION>(1.), static_cast<PRECISION>(-1.), roots), 2);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -1.6180339887498948482045868343656381177203091798057638l, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], 0.61803398874989484820458683436563811772030917980576384l, 1e-15);
}

BOOST_AUTO_TEST_CASE(PolySolveCubicTest) {
	PRECISION roots[3];
	BOOST_CHECK_EQUAL(polySolveCubic(static_cast<PRECISION>(5.), static_cast<PRECISION>(1.), static_cast<PRECISION>(-1.), roots), 3);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -4.74482607768192328566, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], -0.604068139818793694131, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[2], 0.348894217500716979788, 1e-15);
}
