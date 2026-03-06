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
#include <fmt/core.h>

#include "Utility.hpp"
#include "metric/metric.hpp"
#include "view/view.hpp"

using namespace std;
using namespace sbody;

#define PRECISION long double

BOOST_AUTO_TEST_CASE(TestCarlsonRC) {
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRC(static_cast<PRECISION>(0.1), static_cast<PRECISION>(-0.2)), 1.2022125790504758, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestCarlsonRJ) {
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRJ(0.2l, 0.1l, 0.3l, 0.4l), 7.584662876718954, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestEllint_1) {
	const PRECISION phi = 0.5l;
	const PRECISION k = 0.2l;
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(k, phi), 0.5007959922424561, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(k), 1.5868678474541662, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestEllint_2) {
	const PRECISION phi = 0.5l;
	const PRECISION k = 0.2l;
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_2(k, phi), 0.49920624167354444, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_2(k), 1.5549685462425293, 1e-15);
}

// BOOST_AUTO_TEST_CASE(TestEllint_3) {
// 	const double phi = 0.5l;
// 	const double k = 0.2l;
// 	const double n = 0.3l;
// 	BOOST_CHECK_CLOSE_FRACTION(gsl_sf_ellint_P(phi, k, n, GSL_PREC_DOUBLE), boost::math::ellint_3(k, -n, phi), 1e-15);
// 	BOOST_CHECK_CLOSE_FRACTION(gsl_sf_ellint_Pcomp(k, n, GSL_PREC_DOUBLE), boost::math::ellint_3(k, -n), 1e-15);
// }

BOOST_AUTO_TEST_CASE(TestJacobiElliptic) {
	const double k = 0.2l;
	const double u = 0.5l;
	double cn, dn;
	double sn_gsl, cn_gsl, dn_gsl;
	// gsl_sf_elljac_e(u,power2(k), &sn_gsl, &cn_gsl, &dn_gsl);
	// BOOST_CHECK_CLOSE_FRACTION(sn_gsl, boost::math::jacobi_elliptic(k, u, &cn, &dn), 1e-15);
	// BOOST_CHECK_CLOSE_FRACTION(cn_gsl, cn, 1e-15);
	// BOOST_CHECK_CLOSE_FRACTION(dn_gsl, dn, 1e-15);
}
