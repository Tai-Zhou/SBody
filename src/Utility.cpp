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

#include <gsl/gsl_errno.h>

namespace SBody {
	double absAcc = 1e-15, relAcc = 1e-15;
	Integrator::Integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), int coordinate, void *params, const gsl_odeiv2_step_type *type) : coordinate_(coordinate), control_(gsl_odeiv2_control_y_new(absAcc, relAcc)), evolve_(gsl_odeiv2_evolve_alloc(8)), step_(gsl_odeiv2_step_alloc(type, 8)) {
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
	void Cross(const double x[], const double y[], double z[]) {
		z[0] = x[1] * y[2] - x[2] * y[1];
		z[1] = x[2] * y[0] - x[0] * y[2];
		z[2] = x[0] * y[1] - x[1] * y[0];
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
