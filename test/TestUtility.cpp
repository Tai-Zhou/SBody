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

#define PRECISION long double

BOOST_AUTO_TEST_CASE(TestPolySolveQuadratic) {
	PRECISION roots[2];
	BOOST_CHECK_EQUAL(PolySolveQuadratic(static_cast<PRECISION>(1.), static_cast<PRECISION>(1.), static_cast<PRECISION>(-1.), roots), 2);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -1.6180339887498948482045868343656381177203091798057638l, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], 0.61803398874989484820458683436563811772030917980576384l, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestPolySolveCubic) {
	PRECISION roots[3];
	BOOST_CHECK_EQUAL(PolySolveCubic(static_cast<PRECISION>(5.), static_cast<PRECISION>(1.), static_cast<PRECISION>(-1.), roots), 3);
	BOOST_CHECK_CLOSE_FRACTION(roots[0], -4.74482607768192328566, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[1], -0.604068139818793694131, 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(roots[2], 0.348894217500716979788, 1e-15);
}

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
	// gsl_sf_elljac_e(u, Power2(k), &sn_gsl, &cn_gsl, &dn_gsl);
	// BOOST_CHECK_CLOSE_FRACTION(sn_gsl, boost::math::jacobi_elliptic(k, u, &cn, &dn), 1e-15);
	// BOOST_CHECK_CLOSE_FRACTION(cn_gsl, cn, 1e-15);
	// BOOST_CHECK_CLOSE_FRACTION(dn_gsl, dn, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestLUDecomposition) {
	namespace ublas = boost::numeric::ublas;
	ublas::bounded_matrix<PRECISION, 3, 3> A, LU;
	ublas::permutation_matrix permutation(3);
	ublas::bounded_vector<PRECISION, 3> b, x;
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
	auto func = [](const boost::numeric::ublas::bounded_vector<PRECISION, 2> &x, boost::numeric::ublas::bounded_vector<PRECISION, 2> &f, double &_) -> int {
		f(0) = x(0) * x(0) * x(1) + std::sin(x(1));
		f(1) = x(0) * x(1) * x(1) + std::cos(x(0));
		return StatusCode::SUCCESS;
	};
	double no_use = 0.;
	DNewtonMultiFunctionSolver<PRECISION, 2, double> solver(func, no_use);
	boost::numeric::ublas::bounded_vector<PRECISION, 2> x, x_dir, f;
	x(0) = 0.5;
	x(1) = 1.;
	x_dir(0) = 0.5;
	x_dir(1) = 0.5;
	f(0) = 0.25 + std::sin(1);
	f(1) = 0.5 + std::cos(0.5);
	boost::numeric::ublas::bounded_matrix<PRECISION, 2, 2> jacobian, directional_jacobian;
	PRECISION relative_epsilon = 1e-7;
	solver.TwoSidedJacobian(x, f, relative_epsilon, jacobian);
	solver.TwoSidedDirectionalJacobian(x, x_dir, f, relative_epsilon, directional_jacobian);
	BOOST_CHECK_CLOSE_FRACTION(jacobian(0, 0), directional_jacobian(0, 0), 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(jacobian(0, 1), directional_jacobian(0, 1), 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(jacobian(1, 0), directional_jacobian(1, 0), 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(jacobian(1, 1), directional_jacobian(1, 1), 1e-10);
}

BOOST_AUTO_TEST_CASE(TestCoordinateOrthogonalization) {
	BOOST_CHECK_CLOSE_FRACTION(CarlsonRJ(0.2l, 0.1l, 0.3l, 0.4l), 7.584662876718954, 1e-15);
}

BOOST_AUTO_TEST_CASE(TestRotationAngleOfPerpendicularCoordiante) {
	double sin_theta_observer = M_SQRT1_2;
	double cos_theta_observer = M_SQRT1_2;
	double iota, sin_iota, cos_iota, phi;
	int status = RotationAngleOfPerpendicularCoordiante(sin_theta_observer, cos_theta_observer, 1., 0., M_SQRT1_2, -M_SQRT1_2, iota, sin_iota, cos_iota, phi);
	BOOST_CHECK_EQUAL(status, StatusCode::SUCCESS);
	BOOST_CHECK_CLOSE_FRACTION(iota, std::atan2(0.5, M_SQRT1_2), 1e-15);
	BOOST_CHECK_CLOSE_FRACTION(phi, boost::math::constants::two_thirds_pi<double>(), 1e-15);
}
