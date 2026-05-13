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

#include "math/solver/system/gauss_newton.hpp"
#include "math/solver/system/nonlinear_system.hpp"
#include "metric/metric.hpp"
#include "view/view.hpp"

using namespace std;
using namespace sbody;

BOOST_AUTO_TEST_CASE(TestLUDecomposition) {
	namespace ublas = boost::numeric::ublas;
	ublas::bounded_matrix<long double, 3, 3> A, LU;
	ublas::permutation_matrix permutation(3);
	ublas::bounded_vector<long double, 3> b, x;
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
	auto func = [](const boost::numeric::ublas::bounded_vector<long double, 2> &x, boost::numeric::ublas::bounded_vector<long double, 2> &f, double &_) -> int {
		f(0) = x(0) * x(0) * x(1) + std::sin(x(1));
		f(1) = x(0) * x(1) * x(1) + std::cos(x(0));
		return StatusCode::SUCCESS;
	};
	double no_use = 0.;
	GaussNewtonSolver<long double, 2, double> solver(func, no_use);
	boost::numeric::ublas::bounded_vector<long double, 2> x, x_dir, f;
	x(0) = 0.5;
	x(1) = 1.;
	x_dir(0) = 0.5;
	x_dir(1) = 0.5;
	f(0) = 0.25 + std::sin(1);
	f(1) = 0.5 + std::cos(0.5);
	boost::numeric::ublas::bounded_matrix<long double, 2, 2> jacobian, directional_jacobian;
	long double relative_epsilon = 1e-7;
	solver.TwoSidedJacobian(x, f, relative_epsilon, jacobian);
	solver.TwoSidedDirectionalJacobian(x, x_dir, f, relative_epsilon, directional_jacobian);
	BOOST_CHECK_CLOSE_FRACTION(jacobian(0, 0), directional_jacobian(0, 0), 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(jacobian(0, 1), directional_jacobian(0, 1), 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(jacobian(1, 0), directional_jacobian(1, 0), 1e-10);
	BOOST_CHECK_CLOSE_FRACTION(jacobian(1, 1), directional_jacobian(1, 1), 1e-10);
}
