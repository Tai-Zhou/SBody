#include "PostNewtonian.h"

#include <cmath>

#include <gsl/gsl_errno.h>

#include "Constant.h"
#include "Utility.h"

namespace postnewtonian {
	const int dimension = 6;
	int PN = 1;
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
		double r = norm(x);
		double rdot = dot(x, x + 3) / r;
		double vsqr = dot(x + 3);
		double E = vsqr / 2 - constant::G * ((source *)params)->mass * constant::M_sun / r;
		if (PN & 1)
			;
		if (PN & 2)
			;
		if (PN & 8)
			;
		return E;
	}
	double angularMomentum(const double x[], void *params) {
		double J[3], JPN[3];
		cross(x + 1, x + 5, J);
		if (PN & 1) {
		}
		if (PN & 2) {
		}
		if (PN & 8) {
		}
		return norm(J);
	}
} // namespace postnewtonian
