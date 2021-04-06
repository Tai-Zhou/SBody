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
		double M = constant::G * ((source *)params)->mass * constant::M_sun / constant::c2;
		double dydt4 = 2 * M * y[5] / (y[1] - 2 * M) / y[1];
		//d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = dydt4 * y[4];
		//d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = -M * (1 - 2 * M / y[1]) / sqr(y[1]) * constant::c2 + M / (y[1] - 2 * M) / y[1] * sqr(y[5]) + (y[1] - 2 * M) * (sqr(y[6]) + sqr(sin(y[2]) * y[7])) + y[5] * dydt4;
		//d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = -2 / y[1] * y[5] * y[6] + sin(y[2]) * cos(y[2]) * sqr(y[7]) + y[6] * dydt4;
		//d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = -2 * y[7] * (y[5] / y[1] + cos(y[2]) / sin(y[2]) * y[6]) + y[7] * dydt4;
		return GSL_SUCCESS;
	}
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
		return GSL_SUCCESS;
	}
	double energy(const double r[], void *params) {
		return (1 - 2 * constant::G * ((source *)params)->mass * constant::M_sun / constant::c2 / r[1]) / r[4];
	}
	double angularMomentum(const double r[], void *params) {
		return sqr(r[1]) * r[7] / r[4];
	}
	namespace particle {
		int normalization(double y[], void *params) {
			double g00 = 1 - 2 * constant::G * ((source *)params)->mass * constant::M_sun / constant::c2 / y[1];
			if (g00 <= 0)
				return 1;
			y[4] = sqrt(g00 - (sqr(y[5]) / g00 + sqr(y[1] * y[6]) + sqr(y[1] * sin(y[2]) * y[7])) / constant::c2);
			return std::isnan(y[4]);
		}
	} // namespace particle
	namespace light {
		int normalization(double y[], void *params) {
			double g00 = 1 - 2 * constant::G * ((source *)params)->mass * constant::M_sun / constant::c2 / y[1];
			if (g00 <= 0)
				return 1;
			double eff = constant::c * g00 / sqrt(sqr(y[5]) + g00 * (sqr(y[1] * y[6]) + sqr(y[1] * sin(y[2]) * y[7])));
			y[4] = 1; //frequency
			y[5] *= eff;
			y[6] *= eff;
			y[7] *= eff;
			return 0;
		}
	} // namespace light
} // namespace schwarzschild
