#include "KerrH.h"

#include <cmath>

#include <gsl/gsl_errno.h>

#include "Constant.h"
#include "Utility.h"

namespace kerrH {
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
	int qdq2qp(const double r[], double u[], void *params) {
		//[u^t,u^r,u^theta,u^phi]
		u[0] = r[0];
		u[1] = r[1];
		u[2] = r[2];
		u[3] = r[3];
		double m = ((source *)params)->mass;
		double a = ((source *)params)->spin * m;
		double r2 = sqr(r[1]), a2 = sqr(a), sint2 = sqr(sin(r[2]));
		double Delta = r2 - 2. * m * r[1] + a2, rho2 = r2 + a2 * sqr(cos(r[2]));
		double mrrho2 = 2. * m * r[1] / rho2;
		u[4] = (mrrho2 - 1. - mrrho2 * a * sint2 * r[7]) / r[4];
		u[5] = rho2 * r[5] / (Delta * r[4]);
		u[6] = rho2 * r[6] / r[4];
		u[7] = (-mrrho2 * a + (a2 + r2 + mrrho2 * a2 * sint2) * r[7]) * sint2 / r[4];
		return 0;
	}
	int qp2qdq(const double u[], double r[], void *params) {
		r[0] = u[0];
		r[1] = u[1];
		r[2] = u[2];
		r[3] = u[3];
		double m = ((source *)params)->mass;
		double a = ((source *)params)->spin * m;
		double a2 = sqr(a), r2 = sqr(r[1]);
		double Delta = r2 - 2. * m * r[1] + a2, rho2 = r2 + a2 * sqr(cos(r[2]));
		double mrrho2 = 2. * m * r[1] / rho2;
		r[4] = -Delta / ((Delta + mrrho2 * (a2 + r2)) * u[4] + mrrho2 * a * u[7]);
		r[5] = Delta / rho2 * u[5] * r[4];
		r[6] = u[6] / rho2 * r[4];
		r[7] = (-mrrho2 * a * u[4] + (1. - mrrho2) / sqr(sin(r[2])) * u[7]) / Delta * r[4];
		return 0;
	}
	int function(double t, const double y[], double dydt[], void *params) {
		double m = ((source *)params)->mass;
		double a = ((source *)params)->spin * m, r = y[1], sint = sin(y[2]), cost = cos(y[2]), E = -y[4], Theta = sqr(y[6]), L = y[7];
		double a2 = sqr(a), a4 = quad(a), r2 = sqr(y[1]), r4 = quad(y[1]), sint2 = sqr(sint), sint3 = cub(sint), sint4 = quad(sint), cost2 = sqr(cost), cost3 = cub(cost), cost4 = quad(cost);
		double Delta = r2 - 2. * m * r + a2, rho2 = r2 + a2 * cost2, a2r2 = a2 + r2;
		double rho4 = sqr(rho2);
		double Q = Theta + cost2 * (a2 * (1. - sqr(y[4])) + sqr(y[7]) / sint2);
		double R = sqr(y[4] * a2r2 + a * y[7]) - Delta * (r2 + sqr(y[7] + a * y[4]) + Q);
		double dRdr = 4. * y[4] * r * (y[4] * a2r2 + a * y[7]) - 2. * (r - m) * (r2 + sqr(y[7] + a * y[4]) + Q) - 2. * Delta * r;
		//[\tau,r,\theta,\phi,p_t,p_r,p_\theta,p_\phi]
		dydt[0] = rho2 / (-a2r2 * (y[4] * a2r2 + a * y[7]) / Delta + a * (y[7] + a * y[4] * sint2));					//d\tau/dt
		dydt[1] = Delta / rho2 * y[5] * dydt[0];																		//dr/dt
		dydt[2] = y[6] / rho2 * dydt[0];																				//d\theta/dt
		dydt[3] = -(y[4] * a * (a2r2 - Delta) + y[7] * (a2 - Delta * (1. + cost2 / sint2))) / (Delta * rho2) * dydt[0]; //d\phi/dt
		dydt[4] = 0.;
		dydt[5] = ((m * rho2 + r * (Delta - rho2)) * sqr(y[5]) + r * sqr(y[6]) + ((sqr(y[4]) - 1.) * a2 + (2. * sqr(y[4]) - 1.) * r2 + 3. * m * r - r2 - sqr(y[7]) - Q) * r * rho2 / Delta - Theta * r + m * rho2 * (sqr(y[7]) + a2 * sqr(y[4]) + 2. * a * y[4] * y[7] + Q) / Delta - R / sqr(Delta) * ((r - m) * rho2 + Delta * r)) / rho4 * dydt[0];
		dydt[6] = (-(Delta * sqr(y[5]) + sqr(y[6]) - (R + Delta * Theta) / Delta) * a2 / rho2 + a2 * (1. - sqr(y[4])) + sqr(y[7]) / sint2 + cost2 / sint4 * sqr(y[7])) * sint * cost / rho2 * dydt[0];
		dydt[7] = 0.;
		return GSL_SUCCESS;
	}
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
		return GSL_SUCCESS;
	}
	double energy(const double r[], void *params) {
		return -r[4];
	}
	double angularMomentum(const double r[], void *params) {
		return r[7];
	}
	double carter(const double r[], void *params) {
		double m = ((source *)params)->mass;
		double a = ((source *)params)->spin * m;
		double a2 = sqr(a), r2 = sqr(r[1]), cost2 = sqr(cos(r[2]));
		double rho2 = r2 + a2 * cost2;
		double mrrho2 = 2. * m * r[1] / rho2;
		return cost2 * a2 + (sqr(r[6]) + cost2 * (sqr(r[7]) / sqr(sin(r[2])) - a2 * sqr(r[4])));
	}
	namespace particle {
		int normalization(double y[], void *params) {
			double m = ((source *)params)->mass;
			double a = ((source *)params)->spin * m, r = y[1], sint = sin(y[2]);
			double a2 = sqr(a), r2 = sqr(r), sint2 = sqr(sint), sint4 = quad(sint);
			double rho2 = r2 + a2 * sqr(cos(y[2]));
			double mrrho2 = 2. * m * r / rho2;
			y[4] = sqrt(1 - mrrho2 + 2 * mrrho2 * a * sint2 * y[7] - (rho2 / (r2 - 2 * m * r + a2) * sqr(y[5]) + rho2 * sqr(y[6]) + ((a2 + r2) * sint2 + mrrho2 * a2 * sint4) * sqr(y[7])));
			return std::isnan(y[4]);
		}
	} // namespace particle
	namespace light {
		int normalization(double y[], void *params) {
			double m = ((source *)params)->mass;
			double a = ((source *)params)->spin * m, r = y[1], sint = sin(y[2]);
			double a2 = sqr(a), r2 = sqr(r), sint2 = sqr(sint), sint4 = quad(sint);
			double rho2 = r2 + a2 * sqr(cos(y[2]));
			double mrrho2 = 2. * m * r / rho2;
			double effa = rho2 / (r2 - 2. * m * r + a2) * sqr(y[5]) + rho2 * sqr(y[6]) + ((a2 + r2) * sint2 + mrrho2 * a2 * sint4) * sqr(y[7]);
			double effb = -2. * mrrho2 * a * sint2 * y[7];
			double effc = mrrho2 - 1.;
			double eff = 0.5 * (-effb + sqrt(sqr(effb) - 4. * effa * effc)) / effa;
			y[4] = 1.; //frequency
			y[5] *= eff;
			y[6] *= eff;
			y[7] *= eff;
			return 0;
		}
	} // namespace light
} // namespace kerrH
