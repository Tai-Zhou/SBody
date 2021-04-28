#include "PostNewtonian.h"

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
		double r = norm(y);
		double rdot = dot(y, y + 3) / r;
		double vsqr = dot(y + 3);
		double SM = ((source *)params)->mass;
		double F = SM / cub(norm(y));
		dydt[0] = y[3];
		dydt[1] = y[4];
		dydt[2] = y[5];
		dydt[3] = -F * y[0];
		dydt[4] = -F * y[1];
		dydt[5] = -F * y[2];
		double A, B, A1 = 0, B1 = 0, A2 = 0, B2 = 0, A25 = 0, B25 = 0, A3 = 0, B3 = 0, A35 = 0, B35 = 0;
		if (PN & 1) {
			A1 = vsqr - 4 * SM / r;
			B1 = -4 * rdot;
		}
		if (PN & 2) {
			A2 = SM / r * (-2 * sqr(rdot) + 9 * SM / r);
			B2 = 2 * SM / r * rdot;
		}
		if (PN & 8) {
			A3 = -16 * cub(SM / r) + sqr(rdot * SM / r);
			B3 = -4 * rdot * sqr(SM / r);
		}
		A = A1 + A2 + A25 + A3 + A35;
		B = B1 + B2 + B25 + B3 + B35;
		dydt[3] -= SM * (A * y[0] / r + B * y[3]) / sqr(r);
		dydt[4] -= SM * (A * y[1] / r + B * y[4]) / sqr(r);
		dydt[5] -= SM * (A * y[2] / r + B * y[5]) / sqr(r);
		return GSL_SUCCESS;
	}
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
		return GSL_SUCCESS;
	}
	double energy(const double x[], void *params) {
		double r = norm(x);
		double rdot = dot(x, x + 3) / r;
		double vsqr = dot(x + 3);
		double SM = ((source *)params)->mass;
		double E = vsqr / 2 - SM / r;
		if (PN & 1)
			E += 0.5 * sqr(SM / r) + 0.375 * sqr(vsqr) + 1.5 * vsqr * SM / r;
		if (PN & 2)
			E += -0.5 * cub(SM / r) + 5 / 16 * cub(vsqr) + 1.75 * sqr(SM / r) * vsqr + 0.5 * sqr(SM / r) * sqr(rdot) + 2.625 * SM / r * sqr(vsqr);
		if (PN & 8)
			E += 0.375 * quad(SM / r) + 1.25 * cub(SM / r) * vsqr + 1.5 * cub(SM / r) * sqr(rdot) + 35 / 128 * quad(vsqr) + 135 / 16 * sqr(SM / r) * sqr(vsqr) + 0.75 * sqr(SM / r) * vsqr * sqr(rdot) + 55 / 16 * SM / r * cub(vsqr);
		return E;
	}
	double angularMomentum(const double x[], void *params) {
		double r = norm(x);
		double rdot = dot(x, x + 3) / r;
		double vsqr = dot(x + 3);
		double SM = ((source *)params)->mass;
		double J[3], eff = 0;
		cross(x, x + 3, J);
		if (PN & 1)
			eff += 3 * SM / r + 0.5 * vsqr;
		if (PN & 2)
			eff += 3.5 * sqr(SM / r) + 0.375 * sqr(vsqr) + 3.5 * SM / r * vsqr;
		if (PN & 8)
			eff += 2.5 * cub(SM / r) + 5 / 16 * cub(vsqr) + 11.25 * sqr(SM / r) * vsqr + 0.5 * sqr(SM / r) * sqr(rdot) + 4.125 * SM / r * sqr(vsqr);
		for (int i = 0; i < 3; ++i)
			J[i] += J[i] * eff;
		return norm(J);
	}
} // namespace postnewtonian
