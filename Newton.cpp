#include "Newton.h"

#include <gsl/gsl_errno.h>

#include "Constant.h"
#include "Utility.h"

namespace newton {
	const int dimension = 6;
	// from cartesian to spherical
	int c2s(const double x[], double r[]) {
		return ::c2s(x, x + 3, r, r + 3);
	}
	// from spherical to cartesian
	int s2c(const double r[], double x[]) {
		return ::s2c(r, r + 3, x, x + 3);
	}
	int function(double t, const double y[], double dydt[], void *params) {
		dydt[0] = y[3];
		dydt[1] = y[4];
		dydt[2] = y[5];
		double F = constant::G * ((source *)params)->mass * constant::M_sun / cub(norm(y));
		dydt[3] = -F * y[0];
		dydt[4] = -F * y[1];
		dydt[5] = -F * y[2];
		return GSL_SUCCESS;
	}
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
		return GSL_SUCCESS;
	}
	double energy(const double x[], void *params) {
		return dot(x + 3) / 2 - constant::G * ((source *)params)->mass * constant::M_sun / norm(x);
	}
	double angularMomentum(const double x[], void *params) {
		double J[3];
		cross(x, x + 3, J);
		return norm(J);
	}
} // namespace newton
