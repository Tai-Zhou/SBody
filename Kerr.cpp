#include "Kerr.h"

#include <cmath>

#include <gsl/gsl_errno.h>

#include "Constant.h"
#include "Utility.h"

namespace kerr {
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
		double a = ((source *)params)->spin * M, r = y[1], sint = sin(y[2]), cost = cos(y[2]);
		double a2 = sqr(a), a4 = quad(a), r2 = sqr(r), r4 = quad(r), sint2 = sqr(sint), sint4 = quad(sint), cost2 = sqr(cost), cost4 = quad(cost);
		double Delta = r2 - 2 * M * r + a2, rho2 = r2 + a2 * cost2, a2r2 = a2 + r2;
		double rho4 = sqr(rho2), rho6 = cub(rho2), r2rho2 = 2 * r2 - rho2;
		double dydt4 = 2 * M / Delta / rho4 * (a2r2 * r2rho2 * y[5] - 2 * Delta * a2 * r * sint * cost * y[6] * (1 - a / constant::c * sint2 * y[7]) - a / constant::c * (2 * r4 + r2 * rho2 + a2 * r2rho2) * sint2 * y[5] * y[7]);
		//d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = dydt4 * y[4];
		//d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = -Delta * M * r2rho2 / rho6 * constant::c2 - (r / rho2 - (r - M) / Delta) * sqr(y[5]) + 2 * a2 * sint * cost / rho2 * y[5] * y[6] + r * Delta / rho2 * sqr(y[6]) + 2 * Delta * M * a * r2rho2 * sint2 / rho6 * y[7] * constant::c - Delta * sint2 / rho6 * (M * a2 * sint2 * r2rho2 - r * rho4) * sqr(y[7]) + dydt4 * y[5];
		//d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = 2 * M * a2 * r * sint * cost / rho6 * constant::c2 - a2 * sint * cost / Delta / rho2 * sqr(y[5]) - 2 * r / rho2 * y[5] * y[6] + a2 * sint * cost / rho2 * sqr(y[6]) - 4 * M * a * r * sint * cost * a2r2 / rho6 * y[7] * constant::c + sint * cost / rho6 * (2 * M * a4 * r * sint4 + 4 * M * a2 * r * sint2 * rho2 + a2r2 * rho4) * sqr(y[7]) + dydt4 * y[6];
		//d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = (-2 * M * a * r2rho2 / Delta / rho4 * y[5] + 4 * M * a * r * cost / rho4 / sint * y[6]) * constant::c - 2 / Delta / rho4 * (r * rho4 - 2 * M * r2 * rho2 - r2rho2 * M * a2 * sint2) * y[5] * y[7] - 2 * cost / Delta / rho6 / sint * ((rho2 - 2 * M * r) * (a2r2 * rho4 + 4 * M * a2 * r * rho2 * sint2 + 2 * M * a4 * r * sint4) + 4 * sqr(M) * a2 * r2 * a2r2 * sint2) * y[6] * y[7] + dydt4 * y[7];
		return GSL_SUCCESS;
	}
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
		return GSL_SUCCESS;
	}
	double energy(const double r[], void *params) {
		double M = constant::G * ((source *)params)->mass * constant::M_sun / constant::c2;
		double a = ((source *)params)->spin * M;
		double rho2 = sqr(r[1]) + sqr(a) * sqr(cos(r[2]));
		return (2 * M * r[1] / rho2 * (1 - a * sqr(sin(r[2])) * r[7] / constant::c) - 1) / r[4];
	}
	double angularMomentum(const double r[], void *params) {
		double M = constant::G * ((source *)params)->mass * constant::M_sun / constant::c2;
		double a = ((source *)params)->spin * M, sint = sin(r[2]);
		double a2 = sqr(a), r2 = sqr(r[1]), sint2 = sqr(sint);
		double rho2 = r2 + a2 * sqr(cos(r[2]));
		double mrrho2 = 2 * M * r[1] / rho2;
		return (-mrrho2 * a * constant::c + (a2 + r2 + mrrho2 * a2 * sint2) * r[7]) * sint2 / r[4];
	}
	namespace particle {
		int normalization(double y[], void *params) {
			double M = constant::G * ((source *)params)->mass * constant::M_sun / constant::c2;
			double a = ((source *)params)->spin * M, r = y[1], sint = sin(y[2]);
			double a2 = sqr(a), r2 = sqr(r), sint2 = sqr(sint), sint4 = quad(sint);
			double rho2 = r2 + a2 * sqr(cos(y[2]));
			double mrrho2 = 2 * M * r / rho2;
			y[4] = sqrt(1 - mrrho2 + 2 * mrrho2 * a * sint2 * y[7] / constant::c - (rho2 / (r2 - 2 * M * r + a2) * sqr(y[5]) + rho2 * sqr(y[6]) + ((a2 + r2) * sint2 + mrrho2 * a2 * sint4) * sqr(y[7])) / constant::c2);
			return std::isnan(y[4]);
		}
	} // namespace particle
	namespace light {
		int normalization(double y[], void *params) {
			double M = constant::G * ((source *)params)->mass * constant::M_sun / constant::c2;
			double a = ((source *)params)->spin * M, r = y[1], sint = sin(y[2]);
			double a2 = sqr(a), r2 = sqr(r), sint2 = sqr(sint), sint4 = quad(sint);
			double rho2 = r2 + a2 * sqr(cos(y[2]));
			double mrrho2 = 2 * M * r / rho2;
			double effa = rho2 / (r2 - 2 * M * r + a2) * sqr(y[5]) + rho2 * sqr(y[6]) + ((a2 + r2) * sint2 + mrrho2 * a2 * sint4) * sqr(y[7]);
			double effb = -2 * mrrho2 * a * sint2 * y[7] * constant::c;
			double effc = (mrrho2 - 1) * constant::c2;
			double eff = (-effb + sqrt(sqr(effb) - 4 * effa * effc)) / 2 / effa;
			y[4] = 1; //frequency
			y[5] *= eff;
			y[6] *= eff;
			y[7] *= eff;
			return 0;
		}
	} // namespace light
} // namespace kerr
