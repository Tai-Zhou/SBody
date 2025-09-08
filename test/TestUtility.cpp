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

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
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

BOOST_AUTO_TEST_CASE(TestEllint_1) {
	const double phi = 0.5l;
	const double k = 0.2l;
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(k, phi), 0.5007959922424561, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(boost::math::ellint_1(k), 1.5868678474541662, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestEllint_2) {
	const double phi = 0.5l;
	const double k = 0.2l;
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
	// gsl_sf_elljac_e(u, Power2(k), &sn_gsl, &cn_gsl, &dn_gsl);
	// BOOST_CHECK_CLOSE_FRACTION(sn_gsl, boost::math::jacobi_elliptic(k, u, &cn, &dn), 1e-15);
	// BOOST_CHECK_CLOSE_FRACTION(cn_gsl, cn, 1e-15);
	// BOOST_CHECK_CLOSE_FRACTION(dn_gsl, dn, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestLUDecomposition) {
	namespace ublas = boost::numeric::ublas;
	ublas::bounded_matrix<double, 3, 3> A, LU;
	ublas::permutation_matrix permutation(3);
	ublas::bounded_vector<double, 3> b, x;
	A(0, 0) = 4.;
	A(0, 1) = 3.;
	A(0, 2) = 2.;
	A(1, 0) = 3.;
	A(1, 1) = 5.;
	A(1, 2) = 1.;
	A(2, 0) = 2.;
	A(2, 1) = 1.;
	A(2, 2) = 6.;
	b(0) = 1.;
	b(1) = 3.;
	b(2) = 2.;
	LU = A;
	x = b;
	int res = ublas::lu_factorize(LU, permutation);
	BOOST_CHECK_EQUAL(res, 0);
	ublas::lu_substitute(LU, permutation, x);
	x = ublas::prod(A, x);
	for (int i = 0; i < 3; ++i)
		BOOST_CHECK_CLOSE_FRACTION(x(i), b(i), 1e-15);
}

BOOST_AUTO_TEST_CASE(TestTwoSidedJacobian) {
	namespace ublas = boost::numeric::ublas;
	int status = Status::SUCCESS;
	auto func = [](boost::numeric::ublas::bounded_vector<double, 2> &x, boost::numeric::ublas::bounded_vector<double, 2> &f) -> int {
		f(0) = x(0) * x(0) * x(1) + std::sin(x(1));
		f(1) = x(0) * x(1) * x(1) + std::cos(x(0));
		return Status::SUCCESS;
	};
	boost::numeric::ublas::bounded_vector<double, 2> x, f;
	x(0) = 0.5;
	x(1) = 1.;
	f(0) = 0.25 + std::sin(1);
	f(1) = 0.5 + std::cos(0.5);
	boost::numeric::ublas::bounded_matrix<double, 2, 2> jacobian;
	double relative_epsilon = 1e-7;
	ublas::bounded_vector<double, 2>
		x1(x), f_trial_plus, f_trial_minus;
	for (std::size_t j = 0; j < 2; ++j) {
		const double x_j = x(j);
		double dx = std::abs(x_j) >= 1. ? relative_epsilon * x_j : std::copysign(relative_epsilon, x_j);
		while (std::abs(dx) > std::abs(x_j) * boost::math::tools::epsilon<double>()) {
			x1(j) = x_j + dx;
			if (status = func(x1, f_trial_plus); status == Status::SUCCESS)
				break;
			dx *= 0.5;
		}
		if (status != Status::SUCCESS) {
			dx = std::abs(x_j) >= 1. ? relative_epsilon * x_j : std::copysign(relative_epsilon, x_j);
			while (std::abs(dx) > std::abs(x_j) * boost::math::tools::epsilon<double>()) {
				x1(j) = x_j - dx;
				if (status = func(x1, f_trial_minus); status == Status::SUCCESS)
					break;
				dx *= 0.5;
			}
			if (status != Status::SUCCESS)
				ublas::matrix_column(jacobian, j) = (f - f_trial_minus) / dx;
		} else {
			x1(j) = x_j - dx;
			if (status = func(x1, f_trial_minus); status == Status::SUCCESS)
				ublas::matrix_column(jacobian, j) = (f - f_trial_minus) / dx;
			else
				ublas::matrix_column(jacobian, j) = 0.5 * (f_trial_plus - f_trial_minus) / dx;
		}
		x1(j) = x_j;
	}
}

BOOST_AUTO_TEST_CASE(TestCoordinateOrthogonalization) {
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRJ(0.2l, 0.1l, 0.3l, 0.4l), 7.584662876718954, 1e-15);
}
