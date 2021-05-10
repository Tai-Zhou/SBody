#include "PostNewtonian.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

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
		const double m = ((source *)params)->mass, r = norm(y), vsqr = dot(y + 3);
		const double mr = m / r, r2 = gsl_pow_2(r), rdot = dot(y, y + 3) / r;
		const double F = mr / r2, rdot2 = gsl_pow_2(rdot);
		dydt[0] = y[3];
		dydt[1] = y[4];
		dydt[2] = y[5];
		double A1 = 0, B1 = 0, A2 = 0, B2 = 0, A25 = 0, B25 = 0, A3 = 0, B3 = 0, A35 = 0, B35 = 0;
		if (PN & 1) {
			A1 = vsqr - 4. * mr;
			B1 = -4. * rdot;
		}
		if (PN & 2) {
			A2 = mr * (-2. * gsl_pow_2(rdot) + 9. * mr);
			B2 = 2. * mr * rdot;
		}
		if (PN & 8) {
			A3 = -16. * gsl_pow_3(mr) + rdot2 * gsl_pow_2(mr);
			B3 = -4. * rdot * gsl_pow_2(mr);
		}
		const double A = A1 + A2 + A25 + A3 + A35, B = B1 + B2 + B25 + B3 + B35;
		dydt[3] = -F * (y[0] + A * y[0] + B * y[3] * r);
		dydt[4] = -F * (y[1] + A * y[1] + B * y[4] * r);
		dydt[5] = -F * (y[2] + A * y[2] + B * y[5] * r);
		return GSL_SUCCESS;
	}
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
		return GSL_SUCCESS;
	}
	double energy(const double y[], void *params) {
		const double m = ((source *)params)->mass, r = norm(y), vsqr = dot(y + 3);
		const double mr = m / r, rdot = dot(y, y + 3) / r, vsqr2 = gsl_pow_2(vsqr), vsqr3 = gsl_pow_3(vsqr), vsqr4 = gsl_pow_4(vsqr);
		const double mr2 = gsl_pow_2(mr), mr3 = gsl_pow_3(mr), mr4 = gsl_pow_4(mr), rdot2 = gsl_pow_2(rdot);
		double E = 0.5 * vsqr - mr;
		if (PN & 1)
			E += 0.5 * mr2 + 0.375 * vsqr2 + 1.5 * vsqr * mr;
		if (PN & 2)
			E += -0.5 * mr3 + 0.3125 * vsqr3 + 1.75 * mr2 * vsqr + 0.5 * mr2 * rdot2 + 2.625 * mr * vsqr2;
		if (PN & 8)
			E += 0.375 * mr4 + 1.25 * mr3 * vsqr + 1.5 * mr3 * rdot2 + 0.2734375 * vsqr4 + 8.4375 * mr2 * vsqr2 + 0.75 * mr2 * vsqr * rdot2 + 3.4375 * mr * vsqr3;
		return E;
	}
	double angularMomentum(const double y[], void *params) {
		const double m = ((source *)params)->mass, r = norm(y), vsqr = dot(y + 3);
		const double mr = m / r, rdot = dot(y, y + 3) / r, vsqr2 = gsl_pow_2(vsqr), vsqr3 = gsl_pow_3(vsqr);
		const double mr2 = gsl_pow_2(mr), mr3 = gsl_pow_3(mr), rdot2 = gsl_pow_2(rdot);
		double J[3], eff = 0;
		cross(y, y + 3, J);
		if (PN & 1)
			eff += 3. * mr + 0.5 * vsqr;
		if (PN & 2)
			eff += 3.5 * mr2 + 0.375 * vsqr2 + 3.5 * mr * vsqr;
		if (PN & 8)
			eff += 2.5 * mr3 + 0.3125 * vsqr3 + 11.25 * mr2 * vsqr + 0.5 * mr2 * rdot2 + 4.125 * mr * vsqr2;
		for (int i = 0; i < 3; ++i)
			J[i] += J[i] * eff;
		return norm(J);
	}
} // namespace postnewtonian
