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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

namespace SBody {
	double absolute_accuracy = 1e-15, relative_accuracy = 1e-15;
	const double sin_epsilon = sin(epsilon);
	const double cos_epsilon = cos(epsilon);
	Integrator::Integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), int coordinate, void *params, const gsl_odeiv2_step_type *type) : coordinate_(coordinate), control_(gsl_odeiv2_control_y_new(absolute_accuracy, relative_accuracy)), evolve_(gsl_odeiv2_evolve_alloc(8)), step_(gsl_odeiv2_step_alloc(type, 8)) {
		system_ = gsl_odeiv2_system{function, jacobian, 8UL, params};
	}
	Integrator::~Integrator() {
		gsl_odeiv2_control_free(control_);
		gsl_odeiv2_evolve_free(evolve_);
		gsl_odeiv2_step_free(step_);
	}
	int Integrator::Apply(double *t, double t1, double *h, double *y) {
		int status = gsl_odeiv2_evolve_apply(evolve_, control_, step_, &system_, t, t1, h, y);
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
		if (int s = gsl_odeiv2_evolve_reset(evolve_); s)
			return s;
		if (int s = gsl_odeiv2_step_reset(step_); s)
			return s;
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
		if (dimension == 3) {
			if (y == nullptr)
				return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
			return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
		}
		double sum = 0;
		while (dimension--)
			sum += x[dimension] * y[dimension];
		return sum;
	}
	double Norm(const double x[], size_t dimension) {
		if (dimension == 3)
			return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
		double sum = 0;
		while (dimension--)
			sum += x[dimension] * x[dimension];
		return sqrt(sum);
	}
	int Cross(const double x[], const double y[], double z[]) {
		z[0] = x[1] * y[2] - x[2] * y[1];
		z[1] = x[2] * y[0] - x[0] * y[2];
		z[2] = x[0] * y[1] - x[1] * y[0];
		return 0;
	}
	int RotateAroundAxis(double x[], int axis, double angle) {
		const double sin_angle = sin(angle), cos_angle = cos(angle);
		if (axis == 0) {
			angle = x[1] * cos_angle - x[2] * sin_angle;
			x[2] = x[1] * sin_angle + x[2] * cos_angle;
			x[1] = angle;
		} else if (axis == 1) {
			angle = x[2] * cos_angle - x[0] * sin_angle;
			x[0] = x[2] * sin_angle + x[0] * cos_angle;
			x[2] = angle;
		} else if (axis == 2) {
			angle = x[0] * cos_angle - x[1] * sin_angle;
			x[1] = x[0] * sin_angle + x[1] * cos_angle;
			x[0] = angle;
		}
		return 0;
	}
	int CartesianToSpherical(double x[], size_t dimension) {
		std::vector<double> y(x, x + dimension);
		return CartesianToSpherical(y.data(), x, dimension);
	}
	int CartesianToSpherical(const double cartesian[], double spherical[], size_t dimension) {
		if (dimension < 4)
			return 1;
		spherical[1] = Norm(cartesian + 1);
		if (spherical[1] < epsilon)
			return 1;
		spherical[2] = acos(cartesian[3] / spherical[1]);
		if (dimension < 8) {
			spherical[3] = ModBy2Pi(atan2(cartesian[2], cartesian[1]));
		} else {
			spherical[5] = SBody::Dot(cartesian + 1, cartesian + 5) / spherical[1];
			if (const double norm_x_y = Norm(cartesian + 1, 2); norm_x_y < epsilon) {
				spherical[3] = ModBy2Pi(atan2(cartesian[6], cartesian[5]));
				spherical[6] = GSL_SIGN(cartesian[3]) * Norm(cartesian + 5, 2) / spherical[1];
				spherical[7] = 0;
			} else {
				spherical[3] = ModBy2Pi(atan2(cartesian[2], cartesian[1]));
				spherical[6] = (-cartesian[7] + cartesian[3] / spherical[1] * spherical[5]) / norm_x_y;
				spherical[7] = (cartesian[6] * cartesian[1] - cartesian[5] * cartesian[2]) / gsl_pow_2(norm_x_y);
			}
		}
		return 0;
	}
	int CartesianToSpherical(const double cartesian_position[], const double cartesian_velocity[], double spherical_position[], double spherical_velocity[]) {
		spherical_position[0] = Norm(cartesian_position);
		if (spherical_position[0] < epsilon)
			return 1;
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
		return 0;
	}
	int SphericalToCartesian(double x[], size_t dimension) {
		std::vector<double> y(x, x + dimension);
		return SphericalToCartesian(y.data(), x, dimension);
	}
	int SphericalToCartesian(const double spherical[], double cartesian[], size_t dimension) {
		if (dimension < 4)
			return 1;
		const double sin_theta = sin(spherical[2]), cos_theta = cos(spherical[2]), sin_phi = sin(spherical[3]), cos_phi = cos(spherical[3]);
		cartesian[1] = spherical[1] * sin_theta * cos_phi;
		cartesian[2] = spherical[1] * sin_theta * sin_phi;
		cartesian[3] = spherical[1] * cos_theta;
		if (dimension >= 8) {
			cartesian[5] = spherical[5] * sin_theta * cos_phi + spherical[1] * (cos_theta * cos_phi * spherical[6] - sin_theta * sin_phi * spherical[7]);
			cartesian[6] = spherical[5] * sin_theta * sin_phi + spherical[1] * (cos_theta * sin_phi * spherical[6] + sin_theta * cos_phi * spherical[7]);
			cartesian[7] = spherical[5] * cos_theta - spherical[1] * sin_theta * spherical[6];
		}
		return 0;
	}
	int SphericalToCartesian(const double spherical_position[], const double spherical_velocity[], double cartesian_position[], double cartesian_velocity[]) {
		const double sin_theta = sin(spherical_position[1]), cos_theta = cos(spherical_position[1]), sin_phi = sin(spherical_position[2]), cos_phi = cos(spherical_position[2]);
		cartesian_position[0] = spherical_position[0] * sin_theta * cos_phi;
		cartesian_position[1] = spherical_position[0] * sin_theta * sin_phi;
		cartesian_position[2] = spherical_position[0] * cos_theta;
		cartesian_velocity[0] = spherical_velocity[0] * sin_theta * cos_phi + spherical_position[0] * cos_theta * cos_phi * spherical_velocity[1] - spherical_position[0] * sin_theta * sin_phi * spherical_velocity[2];
		cartesian_velocity[1] = spherical_velocity[0] * sin_theta * sin_phi + spherical_position[0] * cos_theta * sin_phi * spherical_velocity[1] + spherical_position[0] * sin_theta * cos_phi * spherical_velocity[2];
		cartesian_velocity[2] = spherical_velocity[0] * cos_theta - spherical_position[0] * sin_theta * spherical_velocity[1];
		return 0;
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
	double _0x(double x) {
		return x > 0 ? x : 0;
	}
	double _0x1(double x) {
		if (x < 0)
			return 0;
		return x < 1 ? x : 1;
	}
} // namespace SBody
