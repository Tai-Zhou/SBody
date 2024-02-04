/**
 * @file Utility.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Utility.h"

#include <cmath>
#include <vector>

#include <fmt/core.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_ellint.h>

#include "IO.h"

namespace SBody {
	double absolute_accuracy = 1e-13, relative_accuracy = 1e-13;
	Integrator::Integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), void *params, const gsl_odeiv2_step_type *type) : control_(gsl_odeiv2_control_y_new(absolute_accuracy, relative_accuracy)), evolve_(gsl_odeiv2_evolve_alloc(8)), step_(gsl_odeiv2_step_alloc(type, 8)) {
		system_ = gsl_odeiv2_system{function, jacobian, 8UL, params};
	}
	Integrator::~Integrator() {
		gsl_odeiv2_control_free(control_);
		gsl_odeiv2_evolve_free(evolve_);
		gsl_odeiv2_step_free(step_);
	}
	int Integrator::Apply(double *t, double t1, double *h, double *y) {
		int status = 0;
		double theta_0 = y[2];
		if (*h > 0)
			while (status <= 0 && *t < t1)
				status = gsl_odeiv2_evolve_apply(evolve_, control_, step_, &system_, t, t1, h, y);
		else
			while (status <= 0 && *t > t1)
				status = gsl_odeiv2_evolve_apply(evolve_, control_, step_, &system_, t, t1, h, y);
		MapTheta(theta_0, y);
		ModBy2Pi(y[3]);
		return status;
	}
	int Integrator::ApplyStep(double *t, double t1, double *h, double *y) {
		double theta_0 = y[2];
		int status = gsl_odeiv2_evolve_apply(evolve_, control_, step_, &system_, t, t1, h, y);
		MapTheta(theta_0, y);
		ModBy2Pi(y[3]);
		return status;
	}
	int Integrator::ApplyFixedStep(double *t, const double h, double *y) {
		double theta_0 = y[2];
		int status = gsl_odeiv2_evolve_apply_fixed_step(evolve_, control_, step_, &system_, t, h, y);
		MapTheta(theta_0, y);
		ModBy2Pi(y[3]);
		return status;
	}
	int Integrator::Reset() {
		if (int status = gsl_odeiv2_evolve_reset(evolve_); status != GSL_SUCCESS)
			return status;
		if (int status = gsl_odeiv2_step_reset(step_); status != GSL_SUCCESS)
			return status;
		return GSL_SUCCESS;
	}
	FunctionSolver::FunctionSolver(const gsl_root_fsolver_type *type) {
		solver_ = gsl_root_fsolver_alloc(type);
	}
	FunctionSolver::FunctionSolver(gsl_function *function, double lower, double upper, const gsl_root_fsolver_type *type) {
		solver_ = gsl_root_fsolver_alloc(type);
		gsl_root_fsolver_set(solver_, function, lower, upper);
	}
	FunctionSolver::~FunctionSolver() {
		gsl_root_fsolver_free(solver_);
	}
	int FunctionSolver::Set(gsl_function *function, double lower, double upper) {
		return gsl_root_fsolver_set(solver_, function, lower, upper);
	}
	int FunctionSolver::Iterate() {
		return gsl_root_fsolver_iterate(solver_);
	}
	int FunctionSolver::Solve() {
		while (gsl_root_test_interval(gsl_root_fsolver_x_lower(solver_), gsl_root_fsolver_x_upper(solver_), absolute_accuracy, relative_accuracy) != GSL_SUCCESS)
			if (int status = gsl_root_fsolver_iterate(solver_); status != GSL_SUCCESS)
				return status;
		return GSL_SUCCESS;
	}
	double FunctionSolver::Root() {
		return gsl_root_fsolver_root(solver_);
	}
	double FunctionSolver::Lower() {
		return gsl_root_fsolver_x_lower(solver_);
	}
	double FunctionSolver::Upper() {
		return gsl_root_fsolver_x_upper(solver_);
	}
	DerivativeSolver::DerivativeSolver(const gsl_root_fdfsolver_type *type) {
		solver_ = gsl_root_fdfsolver_alloc(type);
	}
	DerivativeSolver::DerivativeSolver(gsl_function_fdf *function, double root, const gsl_root_fdfsolver_type *type) {
		solver_ = gsl_root_fdfsolver_alloc(type);
		gsl_root_fdfsolver_set(solver_, function, root);
	}
	int DerivativeSolver::Set(gsl_function_fdf *function, double root) {
		return gsl_root_fdfsolver_set(solver_, function, root);
	}
	DerivativeSolver::~DerivativeSolver() {
		gsl_root_fdfsolver_free(solver_);
	}
	int DerivativeSolver::Iterate() {
		return gsl_root_fdfsolver_iterate(solver_);
	}
	int DerivativeSolver::Solve() {
		return GSL_FAILURE;
	}
	double DerivativeSolver::Root() {
		return gsl_root_fdfsolver_root(solver_);
	}
	GslBlock::~GslBlock() {
		for (auto vector : vectors_)
			gsl_vector_free(vector);
		for (auto matrix : matrices_)
			gsl_matrix_free(matrix);
		for (auto permutation : permutations_)
			gsl_permutation_free(permutation);
	}
	gsl_vector *GslBlock::VectorAlloc(size_t n) {
		vectors_.push_back(gsl_vector_alloc(n));
		return vectors_.back();
	}
	gsl_vector *GslBlock::VectorCalloc(size_t n) {
		vectors_.push_back(gsl_vector_calloc(n));
		return vectors_.back();
	}
	gsl_vector *GslBlock::VectorAllocFromBlock(gsl_block *block, const size_t offset, const size_t n, const size_t stride) {
		vectors_.push_back(gsl_vector_alloc_from_block(block, offset, n, stride));
		return vectors_.back();
	}
	gsl_vector *GslBlock::VectorAllocRowFromMatrix(gsl_matrix *matrix, const size_t i) {
		vectors_.push_back(gsl_vector_alloc_row_from_matrix(matrix, i));
		return vectors_.back();
	}
	gsl_vector *GslBlock::VectorAllocColFromMatrix(gsl_matrix *matrix, const size_t j) {
		vectors_.push_back(gsl_vector_alloc_col_from_matrix(matrix, j));
		return vectors_.back();
	}
	gsl_matrix *GslBlock::MatrixAlloc(size_t n1, size_t n2) {
		matrices_.push_back(gsl_matrix_alloc(n1, n2));
		return matrices_.back();
	}
	gsl_matrix *GslBlock::MatrixCalloc(size_t n1, size_t n2) {
		matrices_.push_back(gsl_matrix_calloc(n1, n2));
		return matrices_.back();
	}
	gsl_permutation *GslBlock::PermutationAlloc(size_t n) {
		permutations_.push_back(gsl_permutation_alloc(n));
		return permutations_.back();
	}
	gsl_permutation *GslBlock::PermutationCalloc(size_t n) {
		permutations_.push_back(gsl_permutation_calloc(n));
		return permutations_.back();
	}
	double Dot(const double x[], const double y[], size_t dimension) {
		if (dimension == 3)
			return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
		double sum = 0;
		while (dimension-- > 0)
			sum += x[dimension] * y[dimension];
		return sum;
	}
	double Dot(const double x[], size_t dimension) {
		if (dimension == 3)
			return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
		double sum = 0;
		while (dimension-- > 0)
			sum += x[dimension] * x[dimension];
		return sum;
	}
	double Norm(const double x[], size_t dimension) {
		if (dimension == 3)
			return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
		double sum = 0;
		while (dimension-- > 0)
			sum += x[dimension] * x[dimension];
		return sqrt(sum);
	}
	int Cross(const double x[], const double y[], double z[]) {
		z[0] = x[1] * y[2] - x[2] * y[1];
		z[1] = x[2] * y[0] - x[0] * y[2];
		z[2] = x[0] * y[1] - x[1] * y[0];
		return GSL_SUCCESS;
	}
	double DotCross(const double x[], const double y[], const double z[]) {
		return x[0] * (y[1] * z[2] - y[2] * z[1]) + x[1] * (y[2] * z[0] - y[0] * z[2]) + x[2] * (y[0] * z[1] - y[1] * z[0]);
	}
	double TriangleArea(const double a, const double b, const double c) {
		const double s = 0.5 * (a + b + c);
		return sqrt(s * (s - a) * (s - b) * (s - c));
	}
	double TriangleArea(const double x[], const double y[], const double z[]) {
		const double a = sqrt(gsl_pow_2(y[0] - x[0]) + gsl_pow_2(y[1] - x[1]) + gsl_pow_2(y[2] - x[2])),
					 b = sqrt(gsl_pow_2(z[0] - y[0]) + gsl_pow_2(z[1] - y[1]) + gsl_pow_2(z[2] - y[2])),
					 c = sqrt(gsl_pow_2(x[0] - z[0]) + gsl_pow_2(x[1] - z[1]) + gsl_pow_2(x[2] - z[2])),
					 s = 0.5 * (a + b + c);
		return sqrt(s * (s - a) * (s - b) * (s - c));
	}
	int RotateAroundAxis(double x[], int axis, double angle) {
		const double sin_angle = sin(angle), cos_angle = cos(angle);
		switch (axis) {
		case 0:
			angle = x[1] * cos_angle - x[2] * sin_angle;
			x[2] = x[1] * sin_angle + x[2] * cos_angle;
			x[1] = angle;
			return GSL_SUCCESS;
		case 1:
			angle = x[2] * cos_angle - x[0] * sin_angle;
			x[0] = x[2] * sin_angle + x[0] * cos_angle;
			x[2] = angle;
			return GSL_SUCCESS;
		case 2:
			angle = x[0] * cos_angle - x[1] * sin_angle;
			x[1] = x[0] * sin_angle + x[1] * cos_angle;
			x[0] = angle;
			return GSL_SUCCESS;
		default:
			return GSL_EINVAL;
		}
	}
	int CartesianToSpherical(double x[], size_t dimension) {
#ifndef GSL_RANGE_CHECK_OFF
		if (dimension != 3 && dimension != 4 && dimension != 8) {
			PrintlnError("CartesianToSpherical dimension = {}", dimension);
			return GSL_EINVAL;
		}
#endif
		if (dimension == 3) {
			const double position[3] = {x[0], x[1], x[2]};
			return CartesianToSpherical(position, x);
		}
		const double position[3] = {x[1], x[2], x[3]};
		if (dimension == 4)
			return CartesianToSpherical(position, x + 1);
		const double velocity[3] = {x[5], x[6], x[7]};
		return CartesianToSpherical(position, velocity, x + 1, x + 5);
	}
	int CartesianToSpherical(const double cartesian[], double spherical[]) {
		if (spherical[0] = Norm(cartesian); spherical[0] < epsilon)
			return GSL_EZERODIV;
		spherical[1] = acos(cartesian[2] / spherical[0]);
		spherical[2] = atan2(cartesian[1], cartesian[0]);
		return GSL_SUCCESS;
	}
	int CartesianToSpherical(const double cartesian_position[], const double cartesian_velocity[], double spherical_position[], double spherical_velocity[]) {
		if (spherical_position[0] = Norm(cartesian_position); spherical_position[0] < epsilon)
			return GSL_EZERODIV;
		spherical_position[1] = acos(cartesian_position[2] / spherical_position[0]);
		spherical_velocity[0] = SBody::Dot(cartesian_position, cartesian_velocity) / spherical_position[0];
		if (const double norm_x_y = Norm(cartesian_position, 2); norm_x_y < epsilon) {
			spherical_position[2] = atan2(cartesian_velocity[1], cartesian_velocity[0]);
			spherical_velocity[1] = GSL_SIGN(cartesian_position[2]) * Norm(cartesian_velocity, 2) / spherical_position[0];
			spherical_velocity[2] = 0;
		} else {
			spherical_position[2] = atan2(cartesian_position[1], cartesian_position[0]);
			spherical_velocity[1] = (-cartesian_velocity[2] + cartesian_position[2] / spherical_position[0] * spherical_velocity[0]) / norm_x_y;
			spherical_velocity[2] = (cartesian_velocity[1] * cartesian_position[0] - cartesian_velocity[0] * cartesian_position[1]) / gsl_pow_2(norm_x_y);
		}
		return GSL_SUCCESS;
	}
	int SphericalToCartesian(double x[], size_t dimension) {
#ifndef GSL_RANGE_CHECK_OFF
		if (dimension != 3 && dimension != 4 && dimension != 8) {
			PrintlnError("SphericalToCartesian dimension = {}", dimension);
			return GSL_EINVAL;
		}
#endif
		if (dimension == 3) {
			const double position[3] = {x[0], x[1], x[2]};
			return SphericalToCartesian(position, x);
		}
		const double position[3] = {x[1], x[2], x[3]};
		if (dimension == 4)
			return SphericalToCartesian(position, x + 1);
		const double velocity[3] = {x[5], x[6], x[7]};
		return SphericalToCartesian(position, velocity, x + 1, x + 5);
	}
	int SphericalToCartesian(const double spherical[], double cartesian[]) {
		const double sin_theta = abs(sin(spherical[1]));
		cartesian[0] = spherical[0] * sin_theta * cos(spherical[2]);
		cartesian[1] = spherical[0] * sin_theta * sin(spherical[2]);
		cartesian[2] = spherical[0] * GSL_SIGN(spherical[1]) * cos(spherical[1]);
		return GSL_SUCCESS;
	}
	int SphericalToCartesian(const double spherical_position[], const double spherical_velocity[], double cartesian_position[], double cartesian_velocity[]) {
		const double sin_theta = abs(sin(spherical_position[1])), cos_theta = GSL_SIGN(spherical_position[1]) * cos(spherical_position[1]), sin_phi = sin(spherical_position[2]), cos_phi = cos(spherical_position[2]);
		if (cartesian_position != nullptr) {
			cartesian_position[0] = spherical_position[0] * sin_theta * cos_phi;
			cartesian_position[1] = spherical_position[0] * sin_theta * sin_phi;
			cartesian_position[2] = spherical_position[0] * cos_theta;
		}
		cartesian_velocity[0] = spherical_velocity[0] * sin_theta * cos_phi + spherical_position[0] * (cos_theta * cos_phi * spherical_velocity[1] - sin_theta * sin_phi * spherical_velocity[2]);
		cartesian_velocity[1] = spherical_velocity[0] * sin_theta * sin_phi + spherical_position[0] * (cos_theta * sin_phi * spherical_velocity[1] + sin_theta * cos_phi * spherical_velocity[2]);
		cartesian_velocity[2] = spherical_velocity[0] * cos_theta - spherical_position[0] * sin_theta * spherical_velocity[1];
		return GSL_SUCCESS;
	}
	int OppositeSign(double x, double y) {
		return (x > 0. && y < 0.) || (x < 0. && y > 0.);
	}
	void MapTheta(const double theta_0, double *y) {
		if (OppositeSign(theta_0, y[2])) {
			y[2] = -y[2];
			y[3] += M_PI;
			y[6] = -y[6];
		} else if (y[2] <= -M_PI_2)
			y[2] += M_PI;
		else if (y[2] > M_PI_2)
			y[2] -= M_PI;
	}
	void ModBy2Pi(double &phi) {
		while (phi < 0)
			phi += M_2PI;
		while (phi >= M_2PI)
			phi -= M_2PI;
	}
	double PhiDifference(double phi) {
		while (phi <= -M_PI)
			phi += M_2PI;
		while (phi > M_PI)
			phi -= M_2PI;
		return phi;
	}
	double EllipticIntegral(double x, double y, double a1, double b1, double a2, double b2, double a3, double b3, double a4, double b4) {
		if (x == y)
			return 0.;
		const double X1 = sqrt(a1 + b1 * x), X2 = sqrt(a2 + b2 * x), X3 = sqrt(a3 + b3 * x), X4 = sqrt(a4 + b4 * x);
		const double Y1 = sqrt(a1 + b1 * y), Y2 = sqrt(a2 + b2 * y), Y3 = sqrt(a3 + b3 * y), Y4 = sqrt(a4 + b4 * y);
		const double U2_12 = gsl_pow_2((X1 * X2 * Y3 * Y4 + Y1 * Y2 * X3 * X4) / (y - x));
		const double U2_13 = U2_12 - (a1 * b4 - a4 * b1) * (a2 * b3 - a3 * b2);
		const double U2_14 = U2_12 - (a1 * b3 - a3 * b1) * (a2 * b4 - a4 * b2);
		return 2. * gsl_sf_ellint_RF(U2_12, U2_13, U2_14, GSL_PREC_DOUBLE);
	}
	double EllipticIntegral2Imaginary(double x, double y, double f, double g, double h, double a1, double b1, double a4, double b4) {
		if (x == y)
			return 0.;
		const double X1 = sqrt(a1 + b1 * x), X4 = sqrt(a4 + b4 * x);
		const double Y1 = sqrt(a1 + b1 * y), Y4 = sqrt(a4 + b4 * y);
		const double xi = sqrt(f + (g + h * y) * y), eta = sqrt(f + (g + h * x) * x);
		const double c11_c44 = 2 * sqrt((f * b1 * b1 - g * a1 * b1 + h * a1 * a1) * (f * b4 * b4 - g * a4 * b4 + h * a4 * a4)), c2_14 = 2. * f * b1 * b4 - g * (a1 * b4 + a4 * b1) + 2. * h * a1 * a4;
		const double M2 = gsl_pow_2(X1 * Y4 + X4 * Y1) * (gsl_pow_2((xi + eta) / (y - x)) - h);
		return 4. * gsl_sf_ellint_RF(M2, M2 + c2_14 - c11_c44, M2 + c2_14 + c11_c44, GSL_PREC_DOUBLE);
	}
	double EllipticIntegral4Imaginary(double x, double y, double f1, double g1, double h1, double f2, double g2, double h2) {
		return GSL_NAN;
	}
	double EllipticIntegral_2(double x, double y, double a5, double b5, double a1, double b1, double a2, double b2, double a3, double b3, double a4, double b4) {
		if (x == y)
			return 0.;
		const double d12 = a1 * b2 - a2 * b1, d13 = a1 * b3 - a3 * b1, d14 = a1 * b4 - a4 * b1, d15_1 = 1. / (a1 * b5 - a5 * b1);
		const double d25 = a2 * b5 - a5 * b2, d35 = a3 * b5 - a5 * b3, d45 = a4 * b5 - a5 * b4;
		const double X1 = sqrt(a1 + b1 * x), X2 = sqrt(a2 + b2 * x), X3 = sqrt(a3 + b3 * x), X4 = sqrt(a4 + b4 * x), X5 = sqrt(a5 + b5 * x);
		const double Y1 = sqrt(a1 + b1 * y), Y2 = sqrt(a2 + b2 * y), Y3 = sqrt(a3 + b3 * y), Y4 = sqrt(a4 + b4 * y), Y5 = sqrt(a5 + b5 * y);
		const double U2_12 = gsl_pow_2((X1 * X2 * Y3 * Y4 + Y1 * Y2 * X3 * X4) / (y - x));
		const double U2_13 = U2_12 - d14 * (a2 * b3 - a3 * b2);
		const double U2_14 = U2_12 - d13 * (a2 * b4 - a4 * b2);
		const double W2 = U2_12 - d13 * d14 * d25 * d15_1;
		const double Q2 = gsl_pow_2(X5 * Y5 / (X1 * Y1)) * W2;
		const double P2 = Q2 + d25 * d35 * d45 * d15_1;
		const double R_C = X1 * Y1 == 0. ? 0. : gsl_sf_ellint_RC(P2, Q2, GSL_PREC_DOUBLE);
		const double I1 = 2. * gsl_sf_ellint_RF(U2_12, U2_13, U2_14, GSL_PREC_DOUBLE);
		const double I3 = 2. * (d12 * d13 * d14 * d15_1 / 3. * gsl_sf_ellint_RJ(U2_12, U2_13, U2_14, W2, GSL_PREC_DOUBLE) + R_C);
		return (b5 * I3 - b1 * I1) * d15_1;
	}
	double EllipticIntegral2Imaginary_2(double x, double y, double a5, double b5, double f, double g, double h, double a1, double b1, double a4, double b4) {
		if (x == y)
			return 0.;
		const double d14 = a1 * b4 - a4 * b1, d15 = a1 * b5 - a5 * b1, d15_1 = 1. / d15, d45 = a4 * b5 - a5 * b4;
		const double X1 = sqrt(a1 + b1 * x), X4 = sqrt(a4 + b4 * x), X5 = sqrt(a5 + b5 * x);
		const double Y1 = sqrt(a1 + b1 * y), Y4 = sqrt(a4 + b4 * y), Y5 = sqrt(a5 + b5 * y);
		const double xi = sqrt(f + (g + h * y) * y), eta = sqrt(f + (g + h * x) * x);
		const double M2 = gsl_pow_2(X1 * Y4 + X4 * Y1) * (gsl_pow_2((xi + eta) / (y - x)) - h);
		const double c2_11 = 2 * (f * b1 * b1 - g * a1 * b1 + h * a1 * a1), c2_44 = 2 * (f * b4 * b4 - g * a4 * b4 + h * a4 * a4), c2_55 = 2 * (f * b5 * b5 - g * a5 * b5 + h * a5 * a5);
		const double c2_14 = 2. * f * b1 * b4 - g * (a1 * b4 + a4 * b1) + 2. * h * a1 * a4, c2_15 = 2. * f * b1 * b5 - g * (a1 * b5 + a5 * b1) + 2. * h * a1 * a5;
		const double c11 = sqrt(c2_11), c44 = sqrt(c2_44), c55 = sqrt(c2_55), c11_c44 = c11 * c44;
		const double L2m = M2 + c2_14 - c11_c44, L2p = M2 + c2_14 + c11_c44;
		const double W2p = M2 + d14 * (c2_15 + c11 * c55) * d15_1;
		const double U = (Y1 * Y4 * eta + X1 * X4 * xi) / (y - x), U2 = gsl_pow_2(U);
		const double W2 = U2 - c2_11 * d45 * d15_1;
		const double Q2 = gsl_pow_2(X5 * Y5 / (X1 * Y1)) * W2;
		const double P2 = Q2 + 0.5 * c2_55 * d45 * d15_1;
		const double I1 = 4. * gsl_sf_ellint_RF(M2, L2m, L2p, GSL_PREC_DOUBLE);
		const double I3 = 2. * (c11 / (3. * c55) * (4. * d14 * d15_1 * (c2_15 + c11 * c55) * gsl_sf_ellint_RJ(M2, L2m, L2p, W2p, GSL_PREC_DOUBLE) - 1.5 * I1 + 3. * gsl_sf_ellint_RC(U2, W2, GSL_PREC_DOUBLE)) + gsl_sf_ellint_RC(P2, Q2, GSL_PREC_DOUBLE));
		return (b5 * I3 - b1 * I1) * d15_1;
	}
	double EllipticIntegral_4(double x, double y, double a5, double b5, double a1, double b1, double a2, double b2, double a3, double b3, double a4, double b4) {
		if (x == y)
			return 0.;
		const double d12 = a1 * b2 - a2 * b1, d13 = a1 * b3 - a3 * b1, d14 = a1 * b4 - a4 * b1, d15 = a1 * b5 - a5 * b1;
		const double d23 = a2 * b3 - a3 * b2, d24 = a2 * b4 - a4 * b2, d34 = a3 * b4 - a4 * b3;
		const double d15_1 = 1. / d15, d25 = a2 * b5 - a5 * b2, d35 = a3 * b5 - a5 * b3, d45 = a4 * b5 - a5 * b4;
		const double X1 = sqrt(a1 + b1 * x), X2 = sqrt(a2 + b2 * x), X3 = sqrt(a3 + b3 * x), X4 = sqrt(a4 + b4 * x), X5 = sqrt(a5 + b5 * x);
		const double Y1 = sqrt(a1 + b1 * y), Y2 = sqrt(a2 + b2 * y), Y3 = sqrt(a3 + b3 * y), Y4 = sqrt(a4 + b4 * y), Y5 = sqrt(a5 + b5 * y);
		const double U2_12 = gsl_pow_2((X1 * X2 * Y3 * Y4 + Y1 * Y2 * X3 * X4) / (y - x));
		const double U2_13 = U2_12 - d14 * d23;
		const double U2_14 = U2_12 - d13 * d24;
		const double W2 = U2_12 - d13 * d14 * d25 * d15_1;
		const double Q2 = gsl_pow_2(X5 * Y5 / (X1 * Y1)) * W2;
		const double P2 = Q2 + d25 * d35 * d45 * d15_1;
		const double R_C = X1 * Y1 == 0. ? 0. : gsl_sf_ellint_RC(P2, Q2, GSL_PREC_DOUBLE);
		const double I1 = 2. * gsl_sf_ellint_RF(U2_12, U2_13, U2_14, GSL_PREC_DOUBLE);
		const double I2 = 2. * (d12 * d13 * gsl_sf_ellint_RD(U2_12, U2_13, U2_14, GSL_PREC_DOUBLE) / 3. + X1 * Y1 / (X4 * Y4 * sqrt(U2_14)));
		const double I3 = 2. * (d12 * d13 * d14 * d15_1 / 3. * gsl_sf_ellint_RJ(U2_12, U2_13, U2_14, W2, GSL_PREC_DOUBLE) + R_C);
		return -0.5 * d15_1 * ((b1 * b5) / d15 + (b2 * b5) / d25 + (b3 * b5) / d35 + (b4 * b5) / d45) * I3 + b5 * b5 * d24 * d34 / (2. * d15 * d25 * d35 * d45) * I2 + gsl_pow_2(b1 * d15_1) * (1. - d12 * d13 * b5 * b5 / (2. * b1 * b1 * d25 * d35)) * I1 - b5 * b5 / (d15 * d25 * d35) * (Y1 * Y2 * Y3 / (Y4 * Y5 * Y5) - X1 * X2 * X3 / (X4 * X5 * X5));
	}
	double EllipticIntegral2Imaginary_4(double x, double y, double a5, double b5, double f, double g, double h, double a1, double b1, double a4, double b4) {
		if (x == y)
			return 0.;
		const double d14 = a1 * b4 - a4 * b1, d15 = a1 * b5 - a5 * b1, d15_1 = 1. / d15, d45 = a4 * b5 - a5 * b4;
		const double X1 = sqrt(a1 + b1 * x), X4 = sqrt(a4 + b4 * x), X5 = sqrt(a5 + b5 * x);
		const double Y1 = sqrt(a1 + b1 * y), Y4 = sqrt(a4 + b4 * y), Y5 = sqrt(a5 + b5 * y);
		const double xi = sqrt(f + (g + h * y) * y), eta = sqrt(f + (g + h * x) * x);
		const double M2 = gsl_pow_2(X1 * Y4 + X4 * Y1) * (gsl_pow_2((xi + eta) / (y - x)) - h);
		const double c2_11 = 2 * (f * b1 * b1 - g * a1 * b1 + h * a1 * a1), c2_44 = 2 * (f * b4 * b4 - g * a4 * b4 + h * a4 * a4), c2_55 = 2 * (f * b5 * b5 - g * a5 * b5 + h * a5 * a5);
		const double c2_14 = 2. * f * b1 * b4 - g * (a1 * b4 + a4 * b1) + 2. * h * a1 * a4, c2_15 = 2. * f * b1 * b5 - g * (a1 * b5 + a5 * b1) + 2. * h * a1 * a5;
		const double c11 = sqrt(c2_11), c44 = sqrt(c2_44), c55 = sqrt(c2_55), c11_c44 = c11 * c44;
		const double L2m = M2 + c2_14 - c11_c44, L2p = M2 + c2_14 + c11_c44;
		const double W2p = M2 + d14 * (c2_15 + c11 * c55) * d15_1;
		const double U = (Y1 * Y4 * eta + X1 * X4 * xi) / (y - x), U2 = gsl_pow_2(U);
		const double W2 = U2 - c2_11 * d45 * d15_1;
		const double Q2 = gsl_pow_2(X5 * Y5 / (X1 * Y1)) * W2;
		const double P2 = Q2 + 0.5 * c2_55 * d45 * d15_1;
		const double I1 = 4. * gsl_sf_ellint_RF(M2, L2m, L2p, GSL_PREC_DOUBLE);
		const double I2 = 2. * c11 / (3. * c44) * (4. * (c2_14 + c11_c44) * gsl_sf_ellint_RD(M2, L2m, L2p, GSL_PREC_DOUBLE) - 1.5 * I1 + 3. / U) + 2. * X1 * Y1 / (X4 * Y4 * U);
		const double I3 = 2. * (c11 / (3. * c55) * (4. * d14 * d15_1 * (c2_15 + c11 * c55) * gsl_sf_ellint_RJ(M2, L2m, L2p, W2p, GSL_PREC_DOUBLE) - 1.5 * I1 + 3. * gsl_sf_ellint_RC(U2, W2, GSL_PREC_DOUBLE)) + gsl_sf_ellint_RC(P2, Q2, GSL_PREC_DOUBLE));
		return -0.5 * d15_1 * ((b1 * b5) / d15 + 2. * b5 * (g * b5 - 2. * h * a5) / c2_55 + (b4 * b5) / d45) * I3 + b5 * b5 * 0.5 * c2_44 / (d15 * d45 * c2_55) * I2 + gsl_pow_2(b1 * d15_1) * (1. - 0.5 * c2_11 * b5 * b5 / (b1 * b1 * c2_55)) * I1 - b5 * b5 / (0.5 * d15 * c2_55) * (Y1 * xi / (Y4 * Y5 * Y5) - X1 * eta / (X4 * X5 * X5));
	}
	double Carlson_RC(double x, double y, gsl_mode_t mode) {
		if (y <= 0.)
			return sqrt(x / (x - y)) * gsl_sf_ellint_RC(x - y, -y, mode);
		return gsl_sf_ellint_RC(x, y, mode);
	}
	double Carlson_RJ(double x, double y, double z, double p, gsl_mode_t mode) {
		if (p <= 0.) {
			const double y_1 = 1. / y, y_p_1 = 1. / (y - p), gamma = y + (z - y) * (y - x) * y_p_1;
			return ((gamma - y) * gsl_sf_ellint_RJ(x, y, z, gamma, mode) - 3. * (gsl_sf_ellint_RF(x, y, z, mode) - Carlson_RC(x * z * y_1, p * gamma * y_1))) * y_p_1;
		}
		return gsl_sf_ellint_RJ(x, y, z, p, mode);
	}
} // namespace SBody
