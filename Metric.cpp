#include "Metric.h"

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include "Constant.h"
#include "Utility.h"

namespace SBody {
	source::source(double mass, double spin) : mass(mass), spin(spin) {}
	namespace Metric {
		int c2s(const double x[], const double v[], double r[], double w[]) {
			// x = {x, y, z}
			// v = {v_x, v_y, v_z}
			// r = {r, \theta, \phi}
			// w = {v_r, v_\theta, v_\phi}
			r[0] = norm(x);
			if (r[0] < Constant::epsilon)
				return 1;
			r[1] = acos(x[2] / r[0]);
			w[0] = dot(x, v) / r[0];
			double normXY = norm(x, 2);
			if (normXY < Constant::epsilon) {
				double normVXY = norm(v, 2);
				if (normVXY < Constant::epsilon) {
					r[2] = 0;
					w[1] = 0;
				}
				else {
					if (v[1] >= 0)
						r[2] = acos(v[0] / normVXY);
					else
						r[2] = 2 * M_PI - acos(v[0] / normVXY);
					if (x[2] >= 0)
						w[1] = normVXY / r[0];
					else
						w[1] = -normVXY / r[0];
				}
				w[2] = 0;
			}
			else {
				if (x[1] >= 0)
					r[2] = acos(x[0] / normXY);
				else
					r[2] = 2 * M_PI - acos(x[0] / normXY);
				w[1] = (-v[2] + x[2] / r[0] * w[0]) / normXY;
				w[2] = (v[1] * x[0] - v[0] * x[1]) / gsl_pow_2(normXY);
			}
			return 0;
		}
		int c2s(const double x[], double r[]) {
			r[0] = x[0];
			r[4] = x[4];
			return c2s(x + 1, x + 5, r + 1, r + 5);
		}
		int s2c(const double r[], const double w[], double x[], double v[]) {
			// r = {r, \theta, \phi}
			// w = {v_r, v_\theta, v_\phi}
			// x = {x, y, z}
			// v = {v_x, v_y, v_z}
			x[0] = r[0] * sin(r[1]) * cos(r[2]);
			x[1] = r[0] * sin(r[1]) * sin(r[2]);
			x[2] = r[0] * cos(r[1]);
			v[0] = w[0] * sin(r[1]) * cos(r[2]) + r[0] * cos(r[1]) * cos(r[2]) * w[1] - r[0] * sin(r[1]) * sin(r[2]) * w[2];
			v[1] = w[0] * sin(r[1]) * sin(r[2]) + r[0] * cos(r[1]) * sin(r[2]) * w[1] + r[0] * sin(r[1]) * cos(r[2]) * w[2];
			v[2] = w[0] * cos(r[1]) - r[0] * sin(r[1]) * w[1];
			return 0;
		}
		int s2c(const double r[], double x[]) {
			x[0] = r[0];
			x[4] = r[4];
			return s2c(r + 1, r + 5, x + 1, x + 5);
		}
		std::array<int (*)(double, const double[], double[], void *), 4> function{Newton::function, Schwarzschild::function, Kerr::function, KerrH::function};
		std::array<int (*)(double, const double[], double *, double[], void *), 4> jacobian{Newton::jacobian, Schwarzschild::jacobian, Kerr::jacobian, KerrH::jacobian};
		std::array<double (*)(const double[], void *), 4> energy{Newton::energy, Schwarzschild::energy, Kerr::energy, KerrH::energy};
		std::array<double (*)(const double[], void *), 4> angularMomentum{Newton::angularMomentum, Schwarzschild::angularMomentum, Kerr::angularMomentum, KerrH::angularMomentum};
		std::array<double (*)(const double[], void *params), 4> carter{Newton::carter, Schwarzschild::carter, Kerr::carter, KerrH::carter};
		namespace Newton {
			int PN = 1;
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
			double carter(const double y[], void *params) { //FIXME: not verified!
				double c[3];
				cross(y + 1, y + 5, c);
				return dot(c) / gsl_pow_2(y[4]);
			}
		} // namespace Newton
		namespace Schwarzschild {
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
				dydt[5] = -m * r2m / gsl_pow_3(r) + 3. * m / r2mr * gsl_pow_2(y[5]) + r2m * (gsl_pow_2(y[6]) + gsl_pow_2(sint * y[7]));
				//d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[6] = -2. * r3m / r2mr * y[5] * y[6] + sint * cost * gsl_pow_2(y[7]);
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
				return gsl_pow_2(r[1]) * r[7] / r[4];
			}
			double carter(const double y[], void *params) {
				return gsl_pow_4(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(y[7] * cos(y[2]) * sin(y[2]))) / gsl_pow_2(y[4]);
			}
			int particleNormalization(double y[], void *params) {
				const double g00 = 1. - 2. * ((source *)params)->mass / y[1];
				if (g00 <= 0)
					return 1;
				y[4] = sqrt(g00 - (gsl_pow_2(y[5]) / g00 + gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
				return std::isnan(y[4]);
			}
			int lightNormalization(double y[], void *params) {
				const double g00 = 1. - 2. * ((source *)params)->mass / y[1];
				if (g00 <= 0)
					return 1;
				const double eff = g00 / sqrt(gsl_pow_2(y[5]) + g00 * (gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
				y[4] = 1.; //frequency
				y[5] *= eff;
				y[6] *= eff;
				y[7] *= eff;
				return 0;
			}
		} // namespace Schwarzschild
		namespace Kerr {
			int function(double t, const double y[], double dydt[], void *params) {
				dydt[0] = y[4]; //d\tau/dt
				dydt[1] = y[5]; //dr/dt
				dydt[2] = y[6]; //d\theta/dt
				dydt[3] = y[7]; //d\phi/dt
				const double m = ((source *)params)->mass, a = ((source *)params)->spin, r = y[1], sint = sin(y[2]), cost = cos(y[2]);
				const double a2 = gsl_pow_2(a), a4 = gsl_pow_4(a), r2 = gsl_pow_2(r), r4 = gsl_pow_4(r), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint), cost2 = gsl_pow_2(cost), cott = cost / sint;
				const double Delta = r2 - 2. * m * r + a2, rho2 = r2 + a2 * cost2, a2r2 = a2 + r2;
				const double rho4 = gsl_pow_2(rho2), rho6 = gsl_pow_3(rho2), r2rho2 = 2. * r2 - rho2;
				const double dydt4 = 2. * m / (Delta * rho4) * (a2r2 * r2rho2 * y[5] - 2. * Delta * a2 * r * sint * cost * y[6] * (1. - a * sint2 * y[7]) - a * (2. * r4 + r2 * rho2 + a2 * r2rho2) * sint2 * y[5] * y[7]);
				//d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
				dydt[4] = dydt4 * y[4];
				//d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[5] = (-Delta * m * r2rho2 - (r - (r - m) * rho2 / Delta) * rho4 * gsl_pow_2(y[5]) + 2. * a2 * sint * cost * rho4 * y[5] * y[6] + r * Delta * rho4 * gsl_pow_2(y[6]) + 2. * Delta * m * a * r2rho2 * sint2 * y[7] - Delta * sint2 * (m * a2 * sint2 * r2rho2 - r * rho4) * gsl_pow_2(y[7])) / rho6 + dydt4 * y[5];
				//d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[6] = (2. * m * a2 * r * sint * cost - a2 * sint * cost * rho4 / Delta * gsl_pow_2(y[5]) - 2. * r * rho4 * y[5] * y[6] + a2 * sint * cost * rho4 * gsl_pow_2(y[6]) - 4. * m * a * r * sint * cost * a2r2 * y[7] + sint * cost * (2. * m * a4 * r * sint4 + 4. * m * a2 * r * sint2 * rho2 + a2r2 * rho4) * gsl_pow_2(y[7])) / rho6 + dydt4 * y[6];
				//d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[7] = (-2. * m * a * r2rho2 * rho2 / Delta * y[5] + 4. * m * a * r * cott * rho2 * y[6] - 2. * rho2 / Delta * (r * rho4 - 2. * m * r2 * rho2 - r2rho2 * m * a2 * sint2) * y[5] * y[7] - 2. * cott / Delta * ((rho2 - 2. * m * r) * (a2r2 * rho4 + 4. * m * a2 * r * rho2 * sint2 + 2. * m * a4 * r * sint4) + 4. * gsl_pow_2(m) * a2 * r2 * a2r2 * sint2) * y[6] * y[7]) / rho6 + dydt4 * y[7];
				return GSL_SUCCESS;
			}
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			double energy(const double r[], void *params) {
				const double m = ((source *)params)->mass, a = ((source *)params)->spin;
				const double rho2 = gsl_pow_2(r[1]) + gsl_pow_2(a) * gsl_pow_2(cos(r[2]));
				return (2. * m * r[1] / rho2 * (1. - a * gsl_pow_2(sin(r[2])) * r[7]) - 1.) / r[4];
			}
			double angularMomentum(const double r[], void *params) {
				const double m = ((source *)params)->mass, a = ((source *)params)->spin;
				const double a2 = gsl_pow_2(a), r2 = gsl_pow_2(r[1]), sint2 = gsl_pow_2(sin(r[2]));
				const double mr_rho2 = 2. * m * r[1] / (r2 + a2 * gsl_pow_2(cos(r[2])));
				return (-mr_rho2 * a + (a2 + r2 + mr_rho2 * a2 * sint2) * r[7]) * sint2 / r[4];
			}
			double carter(const double y[], void *params) {
				const double m = ((source *)params)->mass, a = ((source *)params)->spin;
				const double a2 = gsl_pow_2(a), r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2])), cost2 = gsl_pow_2(cos(y[2]));
				const double rho2 = r2 + a2 * cost2;
				const double mr_rho2 = 2. * m * y[1] / rho2;
				return cost2 * a2 + (gsl_pow_2(rho2 * y[6]) + cost2 * (gsl_pow_2(-mr_rho2 * a + (a2 + r2 + mr_rho2 * a2 * sint2) * y[7]) * sint2 - a2 * gsl_pow_2(mr_rho2 * (1. - a * sint2 * y[7]) - 1.))) / gsl_pow_2(y[4]);
			}
			int particleNormalization(double y[], void *params) {
				const double m = ((source *)params)->mass, a = ((source *)params)->spin, r = y[1], sint = sin(y[2]);
				const double a2 = gsl_pow_2(a), r2 = gsl_pow_2(r), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint);
				const double rho2 = r2 + a2 * gsl_pow_2(cos(y[2]));
				const double mr_rho2 = 2. * m * r / rho2;
				y[4] = sqrt(1. - mr_rho2 + 2. * mr_rho2 * a * sint2 * y[7] - (rho2 / (r2 - 2. * m * r + a2) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2 + r2) * sint2 + mr_rho2 * a2 * sint4) * gsl_pow_2(y[7])));
				return std::isnan(y[4]);
			}
			int lightNormalization(double y[], void *params) {
				const double m = ((source *)params)->mass, a = ((source *)params)->spin, r = y[1], sint = sin(y[2]);
				const double a2 = gsl_pow_2(a), r2 = gsl_pow_2(r), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint);
				const double rho2 = r2 + a2 * gsl_pow_2(cos(y[2]));
				const double mr_rho2 = 2. * m * r / rho2;
				const double effa = rho2 / (r2 - 2. * m * r + a2) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2 + r2) * sint2 + mr_rho2 * a2 * sint4) * gsl_pow_2(y[7]);
				const double effb = -2. * mr_rho2 * a * sint2 * y[7];
				const double effc = mr_rho2 - 1.;
				const double eff = 0.5 * (-effb + sqrt(gsl_pow_2(effb) - 4. * effa * effc)) / effa;
				y[4] = 1.; //frequency
				y[5] *= eff;
				y[6] *= eff;
				y[7] *= eff;
				return 0;
			}
		} // namespace Kerr
		namespace KerrH {
			int qdq2qp(const double r[], double u[], void *params) {
				//[u^t,u^r,u^theta,u^phi]
				u[0] = r[0];
				u[1] = r[1];
				u[2] = r[2];
				u[3] = r[3];
				const double m = ((source *)params)->mass, a = ((source *)params)->spin;
				const double r2 = gsl_pow_2(r[1]), a2 = gsl_pow_2(a), sint2 = gsl_pow_2(sin(r[2]));
				const double Delta = r2 - 2. * m * r[1] + a2, rho2 = r2 + a2 * gsl_pow_2(cos(r[2]));
				const double mr_rho2 = 2. * m * r[1] / rho2;
				u[4] = (mr_rho2 - 1. - mr_rho2 * a * sint2 * r[7]) / r[4];
				u[5] = rho2 * r[5] / (Delta * r[4]);
				u[6] = rho2 * r[6] / r[4];
				u[7] = (-mr_rho2 * a + (a2 + r2 + mr_rho2 * a2 * sint2) * r[7]) * sint2 / r[4];
				return 0;
			}
			int qp2qdq(const double u[], double r[], void *params) {
				r[0] = u[0];
				r[1] = u[1];
				r[2] = u[2];
				r[3] = u[3];
				const double m = ((source *)params)->mass, a = ((source *)params)->spin;
				const double a2 = gsl_pow_2(a), r2 = gsl_pow_2(r[1]);
				const double Delta = r2 - 2. * m * r[1] + a2, rho2 = r2 + a2 * gsl_pow_2(cos(r[2]));
				const double mr_rho2 = 2. * m * r[1] / rho2;
				r[4] = -Delta / ((Delta + mr_rho2 * (a2 + r2)) * u[4] + mr_rho2 * a * u[7]);
				r[5] = Delta / rho2 * u[5] * r[4];
				r[6] = u[6] / rho2 * r[4];
				r[7] = (-mr_rho2 * a * u[4] + (1. - mr_rho2) / gsl_pow_2(sin(r[2])) * u[7]) / Delta * r[4];
				return 0;
			}
			int function(double t, const double y[], double dydt[], void *params) {
				const double m = ((source *)params)->mass, a = ((source *)params)->spin, r = y[1], sint = sin(y[2]), cost = cos(y[2]), Theta = gsl_pow_2(y[6]);
				const double a2 = gsl_pow_2(a), r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint), cost2 = gsl_pow_2(cost);
				const double Delta = r2 - 2. * m * r + a2, rho2 = r2 + a2 * cost2, a2r2 = a2 + r2;
				const double rho4 = gsl_pow_2(rho2);
				const double Q = Theta + cost2 * (a2 * (1. - gsl_pow_2(y[4])) + gsl_pow_2(y[7]) / sint2);
				const double R = gsl_pow_2(y[4] * a2r2 + a * y[7]) - Delta * (r2 + gsl_pow_2(y[7] + a * y[4]) + Q);
				//[\tau,r,\theta,\phi,p_t,p_r,p_\theta,p_\phi]
				dydt[0] = rho2 / (-a2r2 * (y[4] * a2r2 + a * y[7]) / Delta + a * (y[7] + a * y[4] * sint2));					//d\tau/dt
				dydt[1] = Delta / rho2 * y[5] * dydt[0];																		//dr/dt
				dydt[2] = y[6] / rho2 * dydt[0];																				//d\theta/dt
				dydt[3] = -(y[4] * a * (a2r2 - Delta) + y[7] * (a2 - Delta * (1. + cost2 / sint2))) / (Delta * rho2) * dydt[0]; //d\phi/dt
				dydt[4] = 0.;
				dydt[5] = ((m * rho2 + r * (Delta - rho2)) * gsl_pow_2(y[5]) + r * gsl_pow_2(y[6]) + ((gsl_pow_2(y[4]) - 1.) * a2 + (2. * gsl_pow_2(y[4]) - 1.) * r2 + 3. * m * r - r2 - gsl_pow_2(y[7]) - Q) * r * rho2 / Delta - Theta * r + m * rho2 * (gsl_pow_2(y[7]) + a2 * gsl_pow_2(y[4]) + 2. * a * y[4] * y[7] + Q) / Delta - R / gsl_pow_2(Delta) * ((r - m) * rho2 + Delta * r)) / rho4 * dydt[0];
				dydt[6] = (-(Delta * gsl_pow_2(y[5]) + gsl_pow_2(y[6]) - (R + Delta * Theta) / Delta) * a2 / rho2 + a2 * (1. - gsl_pow_2(y[4])) + gsl_pow_2(y[7]) / sint2 + cost2 / sint4 * gsl_pow_2(y[7])) * sint * cost / rho2 * dydt[0];
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
			double carter(const double y[], void *params) {
				return gsl_pow_2(y[6]) + gsl_pow_2(cos(y[2])) * (gsl_pow_2(((source *)params)->spin) * (1. - gsl_pow_2(y[4])) + gsl_pow_2(y[7]) / gsl_pow_2(sin(y[2])));
			}
			int particleNormalization(double y[], void *params) {
				const double m = ((source *)params)->mass, a = ((source *)params)->spin, r = y[1], sint = sin(y[2]);
				const double a2 = gsl_pow_2(a), r2 = gsl_pow_2(r), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint);
				const double rho2 = r2 + a2 * gsl_pow_2(cos(y[2]));
				const double mr_rho2 = 2. * m * r / rho2;
				y[4] = sqrt(1 - mr_rho2 + 2 * mr_rho2 * a * sint2 * y[7] - (rho2 / (r2 - 2 * m * r + a2) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2 + r2) * sint2 + mr_rho2 * a2 * sint4) * gsl_pow_2(y[7])));
				return std::isnan(y[4]);
			}
			int lightNormalization(double y[], void *params) {
				const double m = ((source *)params)->mass, a = ((source *)params)->spin, r = y[1], sint = sin(y[2]);
				const double a2 = gsl_pow_2(a), r2 = gsl_pow_2(r), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint);
				const double rho2 = r2 + a2 * gsl_pow_2(cos(y[2]));
				const double mr_rho2 = 2. * m * r / rho2;
				const double effa = rho2 / (r2 - 2. * m * r + a2) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2 + r2) * sint2 + mr_rho2 * a2 * sint4) * gsl_pow_2(y[7]);
				const double effb = -2. * mr_rho2 * a * sint2 * y[7];
				const double effc = mr_rho2 - 1.;
				const double eff = 0.5 * (-effb + sqrt(gsl_pow_2(effb) - 4. * effa * effc)) / effa;
				y[4] = 1.; //frequency
				y[5] *= eff;
				y[6] *= eff;
				y[7] *= eff;
				return 0;
			}
		} // namespace KerrH
	}	  // namespace Metric
} // namespace SBody
