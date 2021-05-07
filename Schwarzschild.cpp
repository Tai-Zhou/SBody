#include "Schwarzschild.h"

#include <cmath>

#include <gsl/gsl_errno.h>

#include "Constant.h"
#include "Utility.h"

namespace schwarzschild {
	const int dimension = 8;
	// from cartesian to spherical
	int c2s(const double x[], double r[]) {
		r[0] = x[0];
		r[4] = x[4];
		return ::c2s(x + 1, x + 5, r + 1, r + 5);
	}
	// from spherical to cartesian
	int s2c(const double r[], double x[]) {
		x[0] = r[0];
		x[4] = r[4];
		return ::s2c(r + 1, r + 5, x + 1, x + 5);
	}
	int function(double t, const double y[], double dydt[], void *params) {
		dydt[0] = y[4]; //d\tau/dt
		dydt[1] = y[5]; //dr/dt
		dydt[2] = y[6]; //d\theta/dt
		dydt[3] = y[7]; //d\phi/dt
		const double m = ((source *)params)->mass, r = y[1], sint = sin(y[2]), cost = cos(y[2]);
		const double r2m = r - 2. * m, r3m = r - 3. * m;
		const double r2mr = r2m * r;
		//d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = 2. * m * y[5] / r2mr * y[4];
		//d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = -m * r2m / cub(r) + 3. * m / r2mr * sqr(y[5]) + r2m * (sqr(y[6]) + sqr(sint * y[7]));
		//d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = -2. * r3m / r2mr * y[5] * y[6] + sint * cost * sqr(y[7]);
		//d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = -2. * (r3m / r2mr * y[5] + cost / sint * y[6]) * y[7];
		return GSL_SUCCESS;
	}
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
		return GSL_SUCCESS;
	}
	double energy(const double r[], void *params) {
		return (2. * ((source *)params)->mass / r[1] - 1.) / r[4];
	}
	double angularMomentum(const double r[], void *params) {
		return sqr(r[1]) * r[7] / r[4];
	}
	namespace particle {
		int normalization(double y[], void *params) {
			const double g00 = 1. - 2. * ((source *)params)->mass / y[1];
			if (g00 <= 0)
				return 1;
			y[4] = sqrt(g00 - (sqr(y[5]) / g00 + sqr(y[1] * y[6]) + sqr(y[1] * sin(y[2]) * y[7])));
			return std::isnan(y[4]);
		}
	} // namespace particle
	namespace light {
		int normalization(double y[], void *params) {
			const double g00 = 1. - 2. * ((source *)params)->mass / y[1];
			if (g00 <= 0)
				return 1;
			const double eff = g00 / sqrt(sqr(y[5]) + g00 * (sqr(y[1] * y[6]) + sqr(y[1] * sin(y[2]) * y[7])));
			y[4] = 1.; //frequency
			y[5] *= eff;
			y[6] *= eff;
			y[7] *= eff;
			return 0;
		}
	} // namespace light
} // namespace schwarzschild
