#include "Utility.h"

#include <cmath>

#include <gsl/gsl_math.h>

#include "Constant.h"

namespace SBody {
	double absAcc = 1e-15, relAcc = 1e-15;
	integrator::integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), int coordinate, void *params, const gsl_odeiv2_step_type *type) : coordinate(coordinate), type(type), control(gsl_odeiv2_control_y_new(absAcc, relAcc)), evolve(gsl_odeiv2_evolve_alloc(8)), step(gsl_odeiv2_step_alloc(type, 8)) {
		system = {function, jacobian, 8, params};
	}
	integrator::~integrator() {
		gsl_odeiv2_evolve_free(evolve);
	}
	int integrator::apply(double *t, double t1, double *h, double *y) {
		int status = gsl_odeiv2_evolve_apply(evolve, control, step, &system, t, t1, h, y);
		if (coordinate == 2) {
			if (y[2] <= -M_PI_2)
				y[2] += M_PI;
			if (y[2] > M_PI_2)
				y[2] -= M_PI;
		}
		if (coordinate > 0)
			y[3] = mod2Pi(y[3]);
		/*if (!cartesian) {
			y[2] = mod2Pi(y[2]);
			if (y[2] >= M_PI) {
				y[2] = M_2PI - y[2];
				y[6] = -y[6];
				y[3] += M_PI;
			}
			y[3] = mod2Pi(y[3]);
		}*/
		return status;
	}
	int integrator::apply_fixed(double *t, const double h, double *y) {
		int status = gsl_odeiv2_evolve_apply_fixed_step(evolve, control, step, &system, t, h, y);
		if (coordinate == 2) {
			if (y[2] <= -M_PI_2)
				y[2] += M_PI;
			if (y[2] > M_PI_2)
				y[2] -= M_PI;
		}
		if (coordinate > 0)
			y[3] = mod2Pi(y[3]);
		/*if (!cartesian) {
			y[2] = mod2Pi(y[2]);
			if (y[2] >= M_PI) {
				y[2] = M_2PI - y[2];
				y[6] = -y[6];
				y[3] += M_PI;
			}
			y[3] = mod2Pi(y[3]);
		}*/
		return status;
	}
	int integrator::reset() {
		return gsl_odeiv2_evolve_reset(evolve);
	}
	double dot(const double x[], const double y[], size_t dimension) {
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
	double norm(const double x[], size_t dimension) {
		if (dimension == 3)
			return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
		double sum = 0;
		while (dimension--)
			sum += x[dimension] * x[dimension];
		return sqrt(sum);
	}
	void cross(const double x[], const double y[], double z[]) {
		z[0] = x[1] * y[2] - x[2] * y[1];
		z[1] = x[2] * y[0] - x[0] * y[2];
		z[2] = x[0] * y[1] - x[1] * y[0];
	}
	int sign(int x) {
		return x > 0 ? 1 : -1;
	}
	double sign(double x) {
		return x > 0. ? 1. : -1.;
	}
	int oppositeSign(double x, double y) {
		return (x >= 0. && y <= 0.) || (x <= 0. && y >= 0.);
	}
	double mod2Pi(double x) {
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
