#include "Utility.h"

#include <cmath>

#include "Constant.h"

namespace SBody {
	double absAcc = 1e-15, relAcc = 1e-15;
	integrator::integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), int coordinate, double hstart, void *params, const gsl_odeiv2_step_type *type) : coordinate(coordinate) {
		system = gsl_odeiv2_system{function, jacobian, 8UL, params};
		driver = gsl_odeiv2_driver_alloc_y_new(&system, type, hstart, absAcc, relAcc);
	}
	integrator::~integrator() {
		gsl_odeiv2_driver_free(driver);
	}
	int integrator::apply(double *t, double t1, double *y) {
		int status = gsl_odeiv2_driver_apply(driver, t, t1, y);
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
	int integrator::apply_fixed(double *t, const double h, double *y, const unsigned long int n) {
		int status = gsl_odeiv2_driver_apply_fixed_step(driver, t, h, n, y);
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
		gsl_odeiv2_driver_reset(driver);
		return 0;
	}
	int integrator::resetHstart(const double hstart) {
		gsl_odeiv2_driver_reset_hstart(driver, hstart);
		return 0;
	}
	double integrator::getH() {
		return driver->h;
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
		return x >= 0. ? 1. : -1.;
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
