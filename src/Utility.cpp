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

#include "IO.h"

namespace SBody {
	double absolute_accuracy = 1e-15, relative_accuracy = 1e-15;
	Integrator::Integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), int coordinate, void *params, const gsl_odeiv2_step_type *type) : coordinate_(coordinate), control_(gsl_odeiv2_control_y_new(absolute_accuracy, relative_accuracy)), evolve_(gsl_odeiv2_evolve_alloc(8)), step_(gsl_odeiv2_step_alloc(type, 8)) {
		system_ = gsl_odeiv2_system{function, jacobian, 8UL, params};
	}
	Integrator::~Integrator() {
		gsl_odeiv2_control_free(control_);
		gsl_odeiv2_evolve_free(evolve_);
		gsl_odeiv2_step_free(step_);
	}
	int Integrator::Apply(double *t, double t1, double *h, double *y) {
		int status = 0;
		if (*h > 0)
			while (status <= 0 && *t < t1)
				status = gsl_odeiv2_evolve_apply(evolve_, control_, step_, &system_, t, t1, h, y);
		else
			while (status <= 0 && *t > t1)
				status = gsl_odeiv2_evolve_apply(evolve_, control_, step_, &system_, t, t1, h, y);
		if (coordinate_ == 2) {
			if (y[2] <= -M_PI_2)
				y[2] += M_PI;
			if (y[2] > M_PI_2)
				y[2] -= M_PI;
		}
		if (coordinate_ > 0)
			y[3] = ModBy2Pi(y[3]);
		return status;
	}
	int Integrator::ApplyStep(double *t, double t1, double *h, double *y) {
		int status = gsl_odeiv2_evolve_apply(evolve_, control_, step_, &system_, t, t1, h, y);
		if (coordinate_ == 2) {
			if (y[2] <= -M_PI_2)
				y[2] += M_PI;
			if (y[2] > M_PI_2)
				y[2] -= M_PI;
		}
		if (coordinate_ > 0)
			y[3] = ModBy2Pi(y[3]);
		return status;
	}
	int Integrator::ApplyFixedStep(double *t, const double h, double *y) {
		int status = gsl_odeiv2_evolve_apply_fixed_step(evolve_, control_, step_, &system_, t, h, y);
		if (coordinate_ == 2) {
			if (y[2] <= -M_PI_2)
				y[2] += M_PI;
			if (y[2] > M_PI_2)
				y[2] -= M_PI;
		}
		if (coordinate_ > 0)
			y[3] = ModBy2Pi(y[3]);
		/*if (!cartesian) {
			y[2] = ModBy2Pi(y[2]);
			if (y[2] >= M_PI) {
				y[2] = M_2PI - y[2];
				y[6] = -y[6];
				y[3] += M_PI;
			}
			y[3] = ModBy2Pi(y[3]);
		}*/
		return status;
	}
	int Integrator::Reset() {
		if (int status = gsl_odeiv2_evolve_reset(evolve_); status != GSL_SUCCESS)
			return status;
		if (int status = gsl_odeiv2_step_reset(step_); status != GSL_SUCCESS)
			return status;
		return GSL_SUCCESS;
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
		spherical[2] = ModBy2Pi(atan2(cartesian[1], cartesian[0]));
		return GSL_SUCCESS;
	}
	int CartesianToSpherical(const double cartesian_position[], const double cartesian_velocity[], double spherical_position[], double spherical_velocity[]) {
		if (spherical_position[0] = Norm(cartesian_position); spherical_position[0] < epsilon)
			return GSL_EZERODIV;
		spherical_position[1] = acos(cartesian_position[2] / spherical_position[0]);
		spherical_velocity[0] = SBody::Dot(cartesian_position, cartesian_velocity) / spherical_position[0];
		if (const double norm_x_y = Norm(cartesian_position, 2); norm_x_y < epsilon) {
			spherical_position[2] = ModBy2Pi(atan2(cartesian_velocity[1], cartesian_velocity[0]));
			spherical_velocity[1] = GSL_SIGN(cartesian_position[2]) * Norm(cartesian_velocity, 2) / spherical_position[0];
			spherical_velocity[2] = 0;
		} else {
			spherical_position[2] = ModBy2Pi(atan2(cartesian_position[1], cartesian_position[0]));
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
		const double sin_theta = sin(spherical[1]);
		cartesian[0] = spherical[0] * sin_theta * cos(spherical[2]);
		cartesian[1] = spherical[0] * sin_theta * sin(spherical[2]);
		cartesian[2] = spherical[0] * cos(spherical[1]);
		return GSL_SUCCESS;
	}
	int SphericalToCartesian(const double spherical_position[], const double spherical_velocity[], double cartesian_position[], double cartesian_velocity[]) {
		const double sin_theta = sin(spherical_position[1]), cos_theta = cos(spherical_position[1]), sin_phi = sin(spherical_position[2]), cos_phi = cos(spherical_position[2]);
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
		return (x >= 0. && y <= 0.) || (x <= 0. && y >= 0.);
	}
	double ModBy2Pi(double x) {
		while (x < 0)
			x += M_2PI;
		while (x >= M_2PI)
			x -= M_2PI;
		return x;
	}
	double PhiDifference(double x) {
		while (x <= -M_PI)
			x += M_2PI;
		while (x > M_PI)
			x -= M_2PI;
		return x;
	}
	double _0x(double x) {
		return x > 0 ? x : 0;
	}
	double _0x1(double x) {
		if (x < 0)
			return 0;
		return x < 1 ? x : 1;
	}
} // namespace SBody
