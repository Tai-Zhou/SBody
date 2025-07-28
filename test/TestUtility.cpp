/**
 * @file TestUtility.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2024-02-21
 *
 * @copyright Copyright (c) 2024
 *
 */

#define BOOST_TEST_MODULE TestUtility

#include <array>
#include <chrono>
#include <csignal>
#include <memory>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include <boost/test/included/unit_test.hpp>
#include <fmt/core.h>

#include "IO.hpp"
#include "Metric.hpp"
#include "Object.hpp"
#include "Unit.hpp"
#include "Utility.hpp"
#include "View.hpp"

using namespace std;
using namespace SBody;

BOOST_AUTO_TEST_CASE(TestPolySolveQuadratic) {
	long double roots[2];
	BOOST_CHECK_EQUAL(PolySolveQuadratic(1.l, 1.l, -1.l, roots), 2);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -1.6180339887498948482045868343656381177203091798057638l, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], 0.61803398874989484820458683436563811772030917980576384l, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestPolySolveCubic) {
	long double roots[3];
	BOOST_CHECK_EQUAL(PolySolveCubic(5.l, 1.l, -1.l, roots), 3);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -4.74482607768192328566, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], -0.604068139818793694131, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[2], 0.348894217500716979788, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestCarlsonRC) {
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRC(0.1l, -0.2l), 1.2022125790504758, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestCarlsonRJ) {
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRJ(0.2l, 0.1l, 0.3l, 0.4l), 7.584662876718954, 1e-15);
}
