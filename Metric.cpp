#include "Metric.h"

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include "Constant.h"
#include "Utility.h"

namespace SBody {
	namespace Metric {
		double m = 1;
		double a = 0, a2 = 0, a4 = 0;
		double e = 0, e2 = 0, e4 = 0;
		double l = 0, l2 = 0, l4 = 0;
		std::string name = "";
		double (*dot)(const double[], const double[], const double[], const size_t) = nullptr;
		double (*ds2)(const double[], const double[], const size_t) = nullptr;
		int (*qdq2qp)(double[]) = nullptr;
		int (*qp2qdq)(double[]) = nullptr;
		int (*function)(double, const double[], double[], void *) = nullptr;
		int (*functionHamiltonian)(double, const double[], double[], void *) = nullptr;
		int (*jacobian)(double, const double[], double *, double[], void *) = nullptr;
		int (*jacobianHamiltonian)(double, const double[], double *, double[], void *) = nullptr;
		double (*energy)(const double[]) = nullptr;
		double (*energyHamiltonian)(const double[]) = nullptr;
		double (*angularMomentum)(const double[]) = nullptr;
		double (*angularMomentumHamiltonian)(const double[]) = nullptr;
		double (*carter)(const double[], const double) = nullptr;
		double (*carterHamiltonian)(const double[], const double) = nullptr;
		int (*particleNormalization)(double[]) = nullptr;
		int (*lightNormalization)(double[], double) = nullptr;
		void setMetric(int NSK, int PN, double mass, double spin, double charge, double NUT) {
			Newton::PN = PN;
			m = mass;
			a = mass * spin;
			a2 = gsl_pow_2(a);
			a4 = gsl_pow_2(a2);
			e = mass * charge;
			e2 = gsl_pow_2(e);
			e4 = gsl_pow_2(e2);
			l = mass * NUT;
			l2 = gsl_pow_2(l);
			l4 = gsl_pow_2(l2);
			switch (NSK) {
			case 0:
				name = "Newton";
				dot = Newton::dot;
				ds2 = Newton::ds2;
				function = Newton::function;
				jacobian = Newton::jacobian;
				energy = Newton::energy;
				angularMomentum = Newton::angularMomentum;
				carter = Newton::carter;
				particleNormalization = Newton::particleNormalization;
				lightNormalization = Newton::lightNormalization;
				break;
			case 1:
				name = "Schwarzschild";
				dot = Schwarzschild::dot;
				ds2 = Schwarzschild::ds2;
				qdq2qp = Schwarzschild::qdq2qp;
				qp2qdq = Schwarzschild::qp2qdq;
				function = Schwarzschild::function;
				functionHamiltonian = Schwarzschild::functionHamiltonian;
				jacobian = Schwarzschild::jacobian;
				jacobianHamiltonian = Schwarzschild::jacobianHamiltonian;
				energy = Schwarzschild::energy;
				energyHamiltonian = Schwarzschild::energyHamiltonian;
				angularMomentum = Schwarzschild::angularMomentum;
				angularMomentumHamiltonian = Schwarzschild::angularMomentumHamiltonian;
				carter = Schwarzschild::carter;
				carterHamiltonian = Schwarzschild::carterHamiltonian;
				particleNormalization = Schwarzschild::particleNormalization;
				lightNormalization = Schwarzschild::lightNormalization;
				break;
			case 2:
				name = "Kerr";
				dot = Kerr::dot;
				ds2 = Kerr::ds2;
				qdq2qp = Kerr::qdq2qp;
				qp2qdq = Kerr::qp2qdq;
				function = Kerr::function;
				functionHamiltonian = Kerr::functionHamiltonian;
				jacobian = Kerr::jacobian;
				jacobianHamiltonian = Kerr::jacobianHamiltonian;
				energy = Kerr::energy;
				energyHamiltonian = Kerr::energyHamiltonian;
				angularMomentum = Kerr::angularMomentum;
				angularMomentumHamiltonian = Kerr::angularMomentumHamiltonian;
				carter = Kerr::carter;
				carterHamiltonian = Kerr::carterHamiltonian;
				particleNormalization = Kerr::particleNormalization;
				lightNormalization = Kerr::lightNormalization;
				break;
			case 3:
				name = "Kerr-Taub-NUT";
				dot = KerrTaubNUT::dot;
				ds2 = KerrTaubNUT::ds2;
				qdq2qp = KerrTaubNUT::qdq2qp;
				qp2qdq = KerrTaubNUT::qp2qdq;
				function = KerrTaubNUT::function;
				functionHamiltonian = KerrTaubNUT::functionHamiltonian;
				jacobian = KerrTaubNUT::jacobian;
				jacobianHamiltonian = KerrTaubNUT::jacobianHamiltonian;
				energy = KerrTaubNUT::energy;
				energyHamiltonian = KerrTaubNUT::energyHamiltonian;
				angularMomentum = KerrTaubNUT::angularMomentum;
				angularMomentumHamiltonian = KerrTaubNUT::angularMomentumHamiltonian;
				carter = KerrTaubNUT::carter;
				carterHamiltonian = KerrTaubNUT::carterHamiltonian;
				particleNormalization = KerrTaubNUT::particleNormalization;
				lightNormalization = KerrTaubNUT::lightNormalization;
			}
		}
		int c2s(const double x[], const double v[], double r[], double w[]) {
			// x = {x, y, z}
			// v = {v_x, v_y, v_z}
			// r = {r, \theta, \phi}
			// w = {v_r, v_\theta, v_\phi}
			r[0] = norm(x);
			if (r[0] < epsilon)
				return 1;
			r[1] = acos(x[2] / r[0]);
			w[0] = SBody::dot(x, v) / r[0];
			const double normXY = norm(x, 2);
			if (normXY < epsilon) {
				double normVXY = norm(v, 2);
				if (normVXY < epsilon) {
					r[2] = 0;
					w[1] = 0;
				}
				else {
					if (v[1] >= 0)
						r[2] = acos(v[0] / normVXY);
					else
						r[2] = 2. * M_PI - acos(v[0] / normVXY);
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
					r[2] = 2. * M_PI - acos(x[0] / normXY);
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
			const double sint = sin(r[1]), cost = cos(r[1]), sinp = sin(r[2]), cosp = cos(r[2]);
			x[0] = r[0] * sint * cosp;
			x[1] = r[0] * sint * sinp;
			x[2] = r[0] * cost;
			v[0] = w[0] * sint * cosp + r[0] * cost * cosp * w[1] - r[0] * sint * sinp * w[2];
			v[1] = w[0] * sint * sinp + r[0] * cost * sinp * w[1] + r[0] * sint * cosp * w[2];
			v[2] = w[0] * cost - r[0] * sint * w[1];
			return 0;
		}
		int s2c(const double r[], double x[]) {
			x[0] = r[0];
			x[4] = r[4];
			return s2c(r + 1, r + 5, x + 1, x + 5);
		}
		namespace Newton {
			int PN = 1;
			double dot(const double g[], const double x[], const double y[], const size_t dimension) {
				if (dimension == 3)
					return x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
				return -x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
			}
			double ds2(const double x[], const double y[], const size_t dimension) {
				if (dimension == 3)
					return gsl_pow_2(x[1] - y[1]) + gsl_pow_2(x[2] - y[2]) + gsl_pow_2(x[3] - y[3]);
				return -gsl_pow_2(x[0] - y[0]) + gsl_pow_2(x[1] - y[1]) + gsl_pow_2(x[2] - y[2]) + gsl_pow_2(x[3] - y[3]);
			}
			int function(double t, const double y[], double dydt[], void *params) {
				const double r = norm(y), r_1 = 1. / r, v2 = SBody::dot(y + 3);
				const double m_r = m * r_1, rdot = SBody::dot(y, y + 3) * r_1;
				const double F = m_r * gsl_pow_2(r_1), m_r2 = gsl_pow_2(m_r), rdot2 = gsl_pow_2(rdot);
				dydt[0] = y[3];
				dydt[1] = y[4];
				dydt[2] = y[5];
				double A = 0, B = 0;
				if (PN & 1) {
					A += v2 - 4. * m_r;
					B += -4. * rdot;
				}
				if (PN & 2) {
					A += m_r * (-2. * rdot2 + 9. * m_r);
					B += 2. * m_r * rdot;
				}
				if (PN & 8) {
					A += -16. * gsl_pow_3(m_r) + rdot2 * m_r2;
					B += -4. * rdot * m_r2;
				}
				dydt[3] = -F * (y[0] + A * y[0] + B * y[3] * r);
				dydt[4] = -F * (y[1] + A * y[1] + B * y[4] * r);
				dydt[5] = -F * (y[2] + A * y[2] + B * y[5] * r);
				return GSL_SUCCESS;
			}
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			double energy(const double y[]) {
				const double r = norm(y), v2 = SBody::dot(y + 3);
				const double m_r = m / r, rdot = SBody::dot(y, y + 3) / r, v4 = gsl_pow_2(v2), v6 = gsl_pow_3(v2), v8 = gsl_pow_2(v4);
				const double m_r2 = gsl_pow_2(m_r), m_r3 = gsl_pow_3(m_r), m_r4 = gsl_pow_4(m_r), rdot2 = gsl_pow_2(rdot);
				double E = 0.5 * v2 - m_r;
				if (PN & 1)
					E += 0.5 * m_r2 + 0.375 * v4 + 1.5 * v2 * m_r;
				if (PN & 2)
					E += -0.5 * m_r3 + 0.3125 * v6 + 1.75 * m_r2 * v2 + 0.5 * m_r2 * rdot2 + 2.625 * m_r * v4;
				if (PN & 8)
					E += 0.375 * m_r4 + 1.25 * m_r3 * v2 + 1.5 * m_r3 * rdot2 + 0.2734375 * v8 + 8.4375 * m_r2 * v4 + 0.75 * m_r2 * v2 * rdot2 + 3.4375 * m_r * v6;
				return E;
			}
			double angularMomentum(const double y[]) {
				const double r = norm(y), r_1 = 1. / r, v2 = SBody::dot(y + 3);
				const double m_r = m * r_1, rdot = SBody::dot(y, y + 3) * r_1, v4 = gsl_pow_2(v2), v6 = gsl_pow_3(v2);
				const double m_r2 = gsl_pow_2(m_r), m_r3 = gsl_pow_3(m_r), rdot2 = gsl_pow_2(rdot);
				double J[3], eff = 0;
				cross(y, y + 3, J);
				if (PN & 1)
					eff += 3. * m_r + 0.5 * v2;
				if (PN & 2)
					eff += 3.5 * m_r2 + 0.375 * v4 + 3.5 * m_r * v2;
				if (PN & 8)
					eff += 2.5 * m_r3 + 0.3125 * v6 + 11.25 * m_r2 * v2 + 0.5 * m_r2 * rdot2 + 4.125 * m_r * v4;
				for (int i = 0; i < 3; ++i)
					J[i] += J[i] * eff;
				return norm(J);
			}
			double carter(const double y[], const double mu2) { //FIXME: not verified!
				double c[3];
				cross(y + 1, y + 5, c);
				return SBody::dot(c) / gsl_pow_2(y[4]);
			}
			int particleNormalization(double y[]) { //TODO: limit the light speed
				if (SBody::dot(y + 3) >= 1)
					return 1;
				return 0;
			}
			int lightNormalization(double y[], double e) {
				const double v_1 = 1. / SBody::norm(y + 3);
				for (int i = 3; i < 6; ++i)
					y[i] *= v_1;
				return 0;
			}
		} // namespace Newton
		namespace Schwarzschild {
			double dot(const double g[], const double x[], const double y[], const size_t dimension) {
				if (dimension == 3)
					return g[1] * x[1] * y[1] / (g[1] - 2. * m) + gsl_pow_2(g[1]) * x[2] * y[2] + gsl_pow_2(g[1] * sin(g[2])) * x[3] * y[3];
				return -(1. - 2. * m / g[1]) * x[0] * y[0] + g[1] * x[1] * y[1] / (g[1] - 2. * m) + gsl_pow_2(g[1]) * x[2] * y[2] + gsl_pow_2(g[1] * sin(g[2])) * x[3] * y[3];
			}
			double ds2(const double x[], const double y[], const size_t dimension) {
				if (dimension == 3)
					return x[1] * gsl_pow_2(x[1] - y[1]) / (x[1] - 2. * m) + gsl_pow_2(x[1] * (x[2] - y[2])) + gsl_pow_2(x[1] * sin(x[2]) * mod2Pi(x[3] - y[3]));
				return -(1. - 2. * m / x[1]) * gsl_pow_2(x[0] - y[0]) + x[1] * gsl_pow_2(x[1] - y[1]) / (x[1] - 2. * m) + gsl_pow_2(x[1] * (x[2] - y[2])) + gsl_pow_2(x[1] * sin(x[2]) * mod2Pi(x[3] - y[3]));
			}
			int qdq2qp(double y[]) {
				const double dtdtau = 1. / y[4], r2 = gsl_pow_2(y[1]);
				y[4] = (y[4] - 1. + 2. * m / y[1]) * dtdtau; // 1 + p_t
				y[5] *= y[1] / (y[1] - 2. * m) * dtdtau;
				y[6] *= r2 * dtdtau;
				y[7] *= r2 * gsl_pow_2(sin(y[2])) * dtdtau;
				return 0;
			}
			int qp2qdq(double y[]) {
				const double r_1 = 1. / y[1], g00 = 1. - 2. * m * r_1;
				y[4] = -g00 * (y[4] - 1.);
				y[5] *= g00 * y[4];
				y[6] *= gsl_pow_2(r_1) * y[4];
				y[7] *= gsl_pow_2(r_1 / sin(y[2])) * y[4];
				return 0;
			}
			int function(double t, const double y[], double dydt[], void *params) {
				dydt[0] = y[4]; //d\tau/dt
				dydt[1] = y[5]; //dr/dt
				dydt[2] = y[6]; //d\theta/dt
				dydt[3] = y[7]; //d\phi/dt
				const double r = y[1], sint = sin(y[2]), cost = cos(y[2]);
				const double r2m = r - 2. * m, r3m = r - 3. * m;
				const double r2mr_1 = 1. / (r2m * r);
				//d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
				dydt[4] = 2. * m * y[5] * r2mr_1 * y[4];
				//d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[5] = -m * r2m / gsl_pow_3(r) + 3. * m * r2mr_1 * gsl_pow_2(y[5]) + r2m * (gsl_pow_2(y[6]) + gsl_pow_2(sint * y[7]));
				//d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[6] = -2. * r3m * r2mr_1 * y[5] * y[6] + sint * cost * gsl_pow_2(y[7]);
				//d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[7] = -2. * (r3m * r2mr_1 * y[5] + cost / sint * y[6]) * y[7];
				return GSL_SUCCESS;
			}
			int functionHamiltonian(double t, const double y[], double dydt[], void *params) {
				const double r_1 = 1. / y[1], r_2 = gsl_pow_2(r_1), g00 = 1. - 2. * m * r_1, E = 1. - y[4], L2 = gsl_pow_2(y[7]);
				const double sint_1 = 1. / sin(y[2]), sint_2 = gsl_pow_2(sint_1);
				//[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
				dydt[0] = g00 / E;						 //d\tau/dt
				dydt[1] = g00 * y[5] * dydt[0];			 //dr/dt
				dydt[2] = y[6] * r_2 * dydt[0];			 //d\theta/dt
				dydt[3] = y[7] * sint_2 * r_2 * dydt[0]; //d\phi/dt
				dydt[4] = 0.;
				dydt[5] = (-m * (gsl_pow_2(y[5]) + gsl_pow_2(E) / gsl_pow_2(g00)) + (gsl_pow_2(y[6]) + L2 * sint_2) * r_1) * r_2 * dydt[0];
				dydt[6] = sint_2 * L2 * cos(y[2]) * sint_1 * r_2 * dydt[0];
				dydt[7] = 0.;
				return GSL_SUCCESS;
			}
			int functionRIAF(double t, const double y[], double dydt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianRIAF(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			double energy(const double y[]) {
				return (y[1] - 2. * m) / (y[1] * y[4]);
			}
			double energyHamiltonian(const double y[]) {
				return 1. - y[4];
			}
			double angularMomentum(const double y[]) {
				return gsl_pow_2(y[1]) * y[7] / y[4];
			}
			double angularMomentumHamiltonian(const double y[]) {
				return y[7];
			}
			double carter(const double y[], const double mu2) {
				return gsl_pow_4(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(y[7] * cos(y[2]) * sin(y[2]))) / gsl_pow_2(y[4]);
			}
			double carterHamiltonian(const double y[], const double mu2) {
				return gsl_pow_2(y[6]) + gsl_pow_2(y[7] / tan(y[2]));
			}
			int particleNormalization(double y[]) {
				const double g00 = 1. - 2. * m / y[1];
				if (g00 <= 0)
					return 1;
				y[4] = sqrt(g00 - (gsl_pow_2(y[5]) / g00 + gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
				return std::isnan(y[4]);
			}
			int lightNormalization(double y[], double e) {
				const double g00 = 1. - 2. * m / y[1];
				if (g00 <= 0)
					return 1;
				const double eff = g00 / sqrt(gsl_pow_2(y[5]) + g00 * (gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
				y[4] = e;
				y[5] *= eff;
				y[6] *= eff;
				y[7] *= eff;
				return 0;
			}
		} // namespace Schwarzschild
		namespace Kerr {
			double dot(const double g[], const double x[], const double y[], const size_t dimension) {
				const double r = x[1], r2 = gsl_pow_2(r);
				const double sint = sin(x[2]), sint2 = gsl_pow_2(sint), cost = cos(x[2]), cost2 = gsl_pow_2(cost);
				const double Delta = r2 - 2. * m * r + a2, rho2 = r2 + a2 * cost2;
				const double mr_rho2 = 2. * m * r / rho2;
				if (dimension == 3)
					return (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + ((r2 + a2) * sint2 + mr_rho2 * a2 * gsl_pow_2(sint2)) * x[3] * y[3];
				return (mr_rho2 - 1.) * x[0] * y[0] - mr_rho2 * a * sint2 * (x[0] * y[3] + x[3] * y[0]) + (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + ((r2 + a2) * sint2 + mr_rho2 * a2 * gsl_pow_2(sint2)) * x[3] * y[3];
			}
			double ds2(const double x[], const double y[], const size_t dimension) {
				const double r = x[1], r2 = gsl_pow_2(r), d0 = x[0] - y[0], d3 = mod2Pi(x[3] - y[3]);
				const double sint = sin(x[2]), sint2 = gsl_pow_2(sint), cost = cos(x[2]), cost2 = gsl_pow_2(cost);
				const double Delta = r2 - 2. * m * r + a2, rho2 = r2 + a2 * cost2;
				const double mr_rho2 = 2. * m * r / rho2;
				if (dimension == 3)
					return (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + ((r2 + a2) * sint2 + 2. * m * r * a2 * gsl_pow_2(sint2) / rho2) * gsl_pow_2(d3);
				return (mr_rho2 - 1.) * gsl_pow_2(d0) - 2. * mr_rho2 * a * sint2 * d0 * d3 + (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + ((r2 + a2) * sint2 + mr_rho2 * a2 * gsl_pow_2(sint2)) * gsl_pow_2(d3);
			}
			int qdq2qp(double y[]) {
				const double dtdtau = 1. / y[4], r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2]));
				const double Delta = r2 - 2. * m * y[1] + a2, rho2 = r2 + a2 * gsl_pow_2(cos(y[2]));
				const double mr_rho2 = 2. * m * y[1] / rho2;
				y[4] = (y[4] - 1. + mr_rho2 * (1. - a * sint2 * y[7])) * dtdtau; // 1 + p_t
				y[5] *= rho2 / Delta * dtdtau;
				y[6] *= rho2 * dtdtau;
				y[7] = (-mr_rho2 * a + (a2 + r2 + mr_rho2 * a2 * sint2) * y[7]) * sint2 * dtdtau;
				return 0;
			}
			int qp2qdq(double y[]) {
				const double pt = y[4] - 1., r2 = gsl_pow_2(y[1]);
				const double Delta = r2 - 2. * m * y[1] + a2, rho2 = r2 + a2 * gsl_pow_2(cos(y[2]));
				const double rho_2 = 1. / rho2, mr_rho2 = 2. * m * y[1] * rho_2;
				y[4] = -Delta / ((Delta + mr_rho2 * (a2 + r2)) * pt + mr_rho2 * a * y[7]);
				y[5] *= Delta * y[4] * rho_2;
				y[6] *= y[4] * rho_2;
				y[7] = (-mr_rho2 * a * pt + (1. - mr_rho2) / gsl_pow_2(sin(y[2])) * y[7]) / Delta * y[4];
				return 0;
			}
			int function(double t, const double y[], double dydt[], void *params) {
				dydt[0] = y[4]; //d\tau/dt
				dydt[1] = y[5]; //dr/dt
				dydt[2] = y[6]; //d\theta/dt
				dydt[3] = y[7]; //d\phi/dt
				const double r = y[1], r2 = gsl_pow_2(r), r4 = gsl_pow_4(r), a2r2 = a2 + r2;
				const double sint = sin(y[2]), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint), cost = cos(y[2]), cost2 = gsl_pow_2(cost), cott = cost / sint;
				const double Delta = r2 - 2. * m * r + a2, Delta_1 = 1. / Delta;
				const double rho2 = r2 + a2 * cost2, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2), r2a2cost2 = r2 - a2 * cost2;
				const double dydt4 = 2. * m * rho_4 * (Delta_1 * a2r2 * r2a2cost2 * y[5] - 2. * a2 * r * sint * cost * y[6] * (1. - a * sint2 * y[7]) - Delta_1 * a * (2. * r4 + r2 * rho2 + a2 * r2a2cost2) * sint2 * y[5] * y[7]);
				//d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
				dydt[4] = dydt4 * y[4];
				//d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[5] = (-Delta * m * r2a2cost2 - (r * (a2 * sint2 - m * r) + m * (a2 * cost2)) * Delta_1 * rho4 * gsl_pow_2(y[5]) + 2. * a2 * sint * cost * rho4 * y[5] * y[6] + r * Delta * rho4 * gsl_pow_2(y[6]) + 2. * Delta * m * a * r2a2cost2 * sint2 * y[7] - Delta * sint2 * (m * a2 * sint2 * r2a2cost2 - r * rho4) * gsl_pow_2(y[7])) * rho_6 + dydt4 * y[5];
				//d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[6] = (2. * m * a2 * r * sint * cost - a2 * sint * cost * rho4 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho4 * y[5] * y[6] - 4. * m * a * r * sint * cost * a2r2 * y[7] + sint * cost * (2. * m * a4 * r * sint4 + 4. * m * a2 * r * sint2 * rho2 + a2r2 * rho4) * gsl_pow_2(y[7])) * rho_6 + dydt4 * y[6];
				//d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[7] = (-2. * m * a * r2a2cost2 * Delta_1 * y[5] + 4. * m * a * r * cott * y[6] - 2. * Delta_1 * (r * rho4 - 2. * m * r2 * rho2 - r2a2cost2 * m * a2 * sint2) * y[5] * y[7] - 2. * cott * (rho4 + 2. * a2 * m * r * sint2) * y[6] * y[7]) * rho_4 + dydt4 * y[7];
				return GSL_SUCCESS;
			}
			int functionHamiltonian(double t, const double y[], double dydt[], void *params) {
				const double r = y[1], r2 = gsl_pow_2(y[1]), a2r2 = a2 + r2, pr2 = gsl_pow_2(y[5]), ptheta2 = gsl_pow_2(y[6]);
				const double E = 1. - y[4], E2 = gsl_pow_2(E), deltaE2 = (2. - y[4]) * y[4], L2 = gsl_pow_2(y[7]);
				const double sint = sin(y[2]), sint2 = gsl_pow_2(sint), sint_2 = 1. / sint2, sint_4 = gsl_pow_2(sint_2), cost = cos(y[2]), cost2 = gsl_pow_2(cost);
				const double Delta = a2r2 - 2. * m * r, Delta_1 = 1. / Delta, Delta_2 = gsl_pow_2(Delta_1);
				const double rho2 = r2 + a2 * cost2, rho_2 = 1. / rho2, rho_4 = gsl_pow_2(rho_2);
				const double Q = ptheta2 + cost2 * (a2 * deltaE2 + L2 * sint_2);
				const double R = -a2r2 * r2 * deltaE2 + E2 * (2. * m * r * a2) + 2. * m * r * r2 - Delta * Q - (r2 - 2. * m * r) * L2 - 4. * m * r * a * E * y[7]; // R = gsl_pow_2(E * a2r2 - a * y[7]) - Delta * (r2 + gsl_pow_2(y[7] - a * E) + Q);
				//[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
				dydt[0] = rho2 * Delta / (a2r2 * (E * a2r2 - a * y[7]) + a * Delta * (y[7] - a * E * sint2));			  //d\tau/dt
				dydt[1] = Delta * rho_2 * y[5] * dydt[0];																  //dr/dt
				dydt[2] = rho_2 * y[6] * dydt[0];																		  //d\theta/dt
				dydt[3] = (2. * E * a * m * r - y[7] * (a2 - Delta * (1. + cost2 * sint_2))) * Delta_1 * rho_2 * dydt[0]; //d\phi/dt
				dydt[4] = 0.;
				dydt[5] = ((m * a2 * cost2 + r * a2 * sint2 - m * r2) * pr2 + ((-2. * deltaE2 * r2 - a2 * deltaE2 + 3. * m * r - L2 - Q) * r + m * (a2 * E2 + L2 - 2. * a * E * y[7] + Q)) * rho2 * Delta_1 - ((r - m) * rho2 + Delta * r) * R * Delta_2) * rho_4 * dydt[0];
				dydt[6] = (-(Delta * pr2 - R * Delta_1) * a2 * rho_2 + a2 * deltaE2 + L2 * sint_2 + cost2 * sint_4 * L2) * sint * cost * rho_2 * dydt[0];
				dydt[7] = 0.;
				return GSL_SUCCESS;
			}
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			double energy(const double y[]) {
				return (1. - 2. * m * y[1] / (gsl_pow_2(y[1]) + a2 * gsl_pow_2(cos(y[2]))) * (1. - a * gsl_pow_2(sin(y[2])) * y[7])) / y[4];
			}
			double energyHamiltonian(const double y[]) {
				return 1. - y[4];
			}
			double angularMomentum(const double y[]) {
				const double r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2]));
				const double mr_rho2 = 2. * m * y[1] / (r2 + a2 * gsl_pow_2(cos(y[2])));
				return (-mr_rho2 * a + (a2 + r2 + mr_rho2 * a2 * sint2) * y[7]) * sint2 / y[4];
			}
			double angularMomentumHamiltonian(const double y[]) {
				return y[7];
			}
			double carter(const double y[], const double mu2) {
				const double r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2])), cost2 = gsl_pow_2(cos(y[2]));
				const double rho2 = r2 + a2 * cost2;
				const double mr_rho2 = 2. * m * y[1] / rho2;
				return mu2 * cost2 * a2 + (gsl_pow_2(rho2 * y[6]) + cost2 * (gsl_pow_2(-mr_rho2 * a + (a2 + r2 + mr_rho2 * a2 * sint2) * y[7]) * sint2 - a2 * gsl_pow_2(mr_rho2 * (1. - a * sint2 * y[7]) - 1.))) / gsl_pow_2(y[4]);
			}
			double carterHamiltonian(const double y[], const double mu2) {
				return gsl_pow_2(y[6]) + gsl_pow_2(cos(y[2])) * (a2 * (mu2 - gsl_pow_2(1. - y[4])) + gsl_pow_2(y[7] / sin(y[2])));
			}
			int particleNormalization(double y[]) {
				const double r = y[1], r2 = gsl_pow_2(r), a2r2 = a2 + r2;
				const double sint = sin(y[2]), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint);
				const double Delta = r2 - 2. * m * r + a2, rho2 = r2 + a2 * gsl_pow_2(cos(y[2]));
				const double mr_rho2 = 2. * m * r / rho2;
				y[7] += 2. * m * a * r / (gsl_pow_2(a2r2) - a2 * Delta * sint2);
				y[4] = sqrt(1. - mr_rho2 + 2. * mr_rho2 * a * sint2 * y[7] - (rho2 / (r2 - 2. * m * r + a2) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + (a2r2 * sint2 + mr_rho2 * a2 * sint4) * gsl_pow_2(y[7])));
				return std::isnan(y[4]);
			}
			int lightNormalization(double y[], double e) {
				const double r = y[1], r2 = gsl_pow_2(r);
				const double sint = sin(y[2]), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint);
				const double rho2 = r2 + a2 * gsl_pow_2(cos(y[2]));
				const double mr_rho2 = 2. * m * r / rho2;
				const double effa = rho2 / (r2 - 2. * m * r + a2) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2 + r2) * sint2 + mr_rho2 * a2 * sint4) * gsl_pow_2(y[7]);
				const double effb = -2. * mr_rho2 * a * sint2 * y[7];
				const double effc = mr_rho2 - 1.;
				const double eff = 0.5 * (-effb + sqrt(gsl_pow_2(effb) - 4. * effa * effc)) / effa;
				y[4] = e;
				y[5] *= eff;
				y[6] *= eff;
				y[7] *= eff;
				return 0;
			}
		} // namespace Kerr
		namespace KerrTaubNUT {
			double dot(const double g[], const double x[], const double y[], const size_t dimension) {
				const double r = g[1], sint = sin(g[2]), cost = cos(g[2]);
				const double r2 = gsl_pow_2(r), sint2 = gsl_pow_2(sint);
				const double Delta = r2 - 2. * m * r - l2 + a2, rho2 = r2 + gsl_pow_2(l + a * cost), chi = a * sint2 - 2. * l * cost;
				const double rho_2 = 1. / rho2;
				if (dimension == 3)
					return (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + (gsl_pow_2(rho2 + a * chi) * sint2 - chi * chi * Delta) * rho_2 * x[3] * y[3];
				return (a2 * sint2 - Delta) * rho_2 * x[0] * y[0] - 2. * ((m * r + l2) * a * sint2 + Delta * l * cost) * rho_2 * (x[0] * y[3] + x[3] * y[0]) + (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + (gsl_pow_2(rho2 + a * chi) * sint2 - chi * chi * Delta) * rho_2 * x[3] * y[3];
			}
			double ds2(const double x[], const double y[], const size_t dimension) {
				const double r = x[1], sint = sin(x[2]), cost = cos(x[2]), d0 = x[0] - y[0], d3 = mod2Pi(x[3] - y[3]);
				const double r2 = gsl_pow_2(r), sint2 = gsl_pow_2(sint);
				const double Delta = r2 - 2. * m * r - l2 + a2, rho2 = r2 + gsl_pow_2(l + a * cost), chi = a * sint2 - 2. * l * cost;
				const double rho_2 = 1. / rho2;
				if (dimension == 3)
					return (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + (gsl_pow_2(rho2 + a * chi) * sint2 - gsl_pow_2(chi) * Delta) * rho_2 * gsl_pow_2(d3);
				return (a2 * sint2 - Delta) * rho_2 * gsl_pow_2(d0) - 4. * ((m * r + l2) * a * sint2 + Delta * l * cost) * rho_2 * d0 * d3 + (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + (gsl_pow_2(rho2 + a * chi) * sint2 - gsl_pow_2(chi) * Delta) * rho_2 * gsl_pow_2(d3);
			}
			int qdq2qp(double y[]) {
				const double dtdtau = 1. / y[4], r = y[1], r2 = gsl_pow_2(r);
				const double sint2 = gsl_pow_2(sin(y[2])), cost = sign(y[2]) * cos(y[2]);
				const double Delta = r2 - 2. * m * r - l2 + a2;
				const double rho2 = r2 + gsl_pow_2(l + a * cost), rho_2 = 1. / rho2;
				y[4] = (y[4] * rho2 - Delta + a2 * sint2 - 2. * ((m * r + l2) * a * sint2 + Delta * l * cost) * y[7]) * rho_2 * dtdtau; // 1 + p_t
				y[5] *= rho2 / Delta * dtdtau;
				y[6] *= rho2 * dtdtau;
				y[7] = (-2. * ((m * r + l2) * a * sint2 + Delta * l * cost) + (gsl_pow_2(r2 + l2 + a2) * sint2 - gsl_pow_2(a * sint2 - 2. * l * cost) * Delta) * y[7]) * rho_2 * dtdtau;
				return 0;
			}
			int qp2qdq(double y[]) {
				const double pt = y[4] - 1., r = y[1], r2 = gsl_pow_2(r);
				const double sint2 = gsl_pow_2(sin(y[2])), cost = sign(y[2]) * cos(y[2]);
				const double Delta = r2 - 2. * m * r - l2 + a2;
				const double rho2 = r2 + gsl_pow_2(l + a * cost), rho_2 = 1. / rho2;
				y[4] = Delta * rho2 * sint2 / ((Delta * gsl_pow_2(a * sint2 - 2. * l * cost) - gsl_pow_2(r2 + l2 + a2) * sint2) * pt - 2. * ((m * r + l2) * a * sint2 + Delta * l * cost) * y[7]);
				y[5] *= Delta * rho_2 * y[4];
				y[6] *= rho_2 * y[4];
				y[7] = (-2. * ((m * r + l2) * a * sint2 + Delta * l * cost) * pt + (Delta - a2 * sint2) * y[7]) / (Delta * rho2 * sint2) * y[4];
				return 0;
			}
			int function(double t, const double y[], double dydt[], void *params) {
				dydt[0] = y[4]; //d\tau/dt
				dydt[1] = y[5]; //dr/dt
				dydt[2] = y[6]; //d\theta/dt
				dydt[3] = y[7]; //d\phi/dt
				const double r = y[1], r2 = gsl_pow_2(r);
				const double sint = abs(sin(y[2])), sint_1 = 1. / sint, sint2 = gsl_pow_2(sint), cost = sign(y[2]) * cos(y[2]), cost2 = gsl_pow_2(cost);
				const double Delta = r2 - 2. * m * r - l2 + a2, Delta_1 = 1. / Delta;
				const double lacost = l + a * cost, lacost2 = gsl_pow_2(lacost);
				const double rho2 = r2 + lacost2, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2);
				const double chi = a * sint2 - 2. * l * cost, rho2achi = r2 + l2 + a2;
				const double dydt4 = 2. * rho_4 * (Delta_1 * rho2achi * ((m * r + l2 - a2 * cost2 + lacost2) * r - m * lacost2) * y[5] * (1. - chi * y[7]) - sint_1 * chi * ((2. * m * r + l * lacost) * lacost - r2 * l) * y[6] * (1. - chi * y[7]) - 2. * Delta_1 * rho2 * r * ((m * r + l2) * a * sint2 + Delta * l * cost) * y[5] * y[7] - sint_1 * rho4 * l * (1. + cost2) * y[6] * y[7]);
				//d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
				dydt[4] = dydt4 * y[4];
				//d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[5] = -Delta * (r * (m * r + 2. * l * lacost) - lacost2 * m) * rho_6 + (r * (m * r + l2 - a2 + lacost2) - m * lacost2) * Delta_1 * rho_2 * gsl_pow_2(y[5]) + 2. * a * sint * lacost * rho_2 * y[5] * y[6] + r * Delta * rho_2 * gsl_pow_2(y[6]) + Delta * chi * ((m * r + l2 - a2 * cost2 + lacost2) * r - m * lacost2) * rho_6 * y[7] * (2. - chi * y[7]) + Delta * r * sint2 * rho_2 * gsl_pow_2(y[7]) + dydt4 * y[5];
				//d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[6] = ((2. * m * r + l * lacost) * lacost - r2 * l) * sint * rho_6 * (a - 2. * rho2achi * y[7]) - a * sint * lacost * rho_2 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho_2 * y[5] * y[6] - (rho2 * (Delta * chi * lacost - gsl_pow_2(rho2achi) * cost) + (Delta * chi - a * rho2achi * sint2) * rho2achi * lacost) * sint * rho_6 * gsl_pow_2(y[7]) + dydt4 * y[6];
				//d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[7] = -2. * a * (2. * (l2 + l * a * cost) * r + (r2 - lacost2) * m) * Delta_1 * rho_4 * y[5] * (1. - chi * y[7]) + 2. * ((2. * m * r + l * lacost) * lacost - r2 * l) * rho_4 * sint_1 * y[6] * (1. - chi * y[7]) - 2. * (1. - a2 * sint2 * Delta_1) * r * rho_2 * y[5] * y[7] - 2. * cost * sint_1 * y[6] * y[7] + dydt4 * y[7];
				return GSL_SUCCESS;
			}
			int functionTau(double t, const double y[], double dydt[], void *params) {
				dydt[0] = y[4]; //dt/d\tau
				dydt[1] = y[5]; //dr/d\tau
				dydt[2] = y[6]; //d\theta/d\tau
				dydt[3] = y[7]; //d\phi/d\tau
				const double r = y[1], r2 = gsl_pow_2(r);
				const double sint = abs(sin(y[2])), sint_1 = 1. / sint, sint2 = gsl_pow_2(sint), cost = sign(y[2]) * cos(y[2]), cost2 = gsl_pow_2(cost);
				const double Delta = r2 - 2. * m * r - l2 + a2, Delta_1 = 1. / Delta;
				const double lacost = l + a * cost, lacost2 = gsl_pow_2(lacost);
				const double rho2 = r2 + lacost2, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2);
				const double chi = a * sint2 - 2. * l * cost, rho2achi = r2 + l2 + a2;
				//d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
				dydt[4] = -2. * rho_4 * (Delta_1 * rho2achi * ((m * r + l2 - a2 * cost2 + lacost2) * r - m * lacost2) * y[5] * (y[4] - chi * y[7]) - sint_1 * chi * ((2. * m * r + l * lacost) * lacost - r2 * l) * y[6] * (y[4] - chi * y[7]) - 2. * Delta_1 * rho2 * r * ((m * r + l2) * a * sint2 + Delta * l * cost) * y[5] * y[7] - sint_1 * rho4 * l * (1. + cost2) * y[6] * y[7]);
				//d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[5] = -Delta * (r * (m * r + 2. * l * lacost) - lacost2 * m) * rho_6 * gsl_pow_2(y[4]) + (r * (m * r + l2 - a2 + lacost2) - m * lacost2) * Delta_1 * rho_2 * gsl_pow_2(y[5]) + 2. * a * sint * lacost * rho_2 * y[5] * y[6] + r * Delta * rho_2 * gsl_pow_2(y[6]) + Delta * chi * ((m * r + l2 - a2 * cost2 + lacost2) * r - m * lacost2) * rho_6 * y[7] * (2. * y[4] - chi * y[7]) + Delta * r * sint2 * rho_2 * gsl_pow_2(y[7]);
				//d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[6] = ((2. * m * r + l * lacost) * lacost - r2 * l) * sint * rho_6 * y[4] * (a * y[4] - 2. * rho2achi * y[7]) - a * sint * lacost * rho_2 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho_2 * y[5] * y[6] - (rho2 * (Delta * chi * lacost - gsl_pow_2(rho2achi) * cost) + (Delta * chi - a * rho2achi * sint2) * rho2achi * lacost) * sint * rho_6 * gsl_pow_2(y[7]);
				//d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
				dydt[7] = -2. * a * (2. * (l2 + l * a * cost) * r + (r2 - lacost2) * m) * Delta_1 * rho_4 * y[5] * (y[4] - chi * y[7]) + 2. * ((2. * m * r + l * lacost) * lacost - r2 * l) * rho_4 * sint_1 * y[6] * (y[4] - chi * y[7]) - 2. * (1. - a2 * sint2 * Delta_1) * r * rho_2 * y[5] * y[7] - 2. * cost * sint_1 * y[6] * y[7];
				return GSL_SUCCESS;
			}
			int functionHamiltonian(double t, const double y[], double dydt[], void *params) { //TODO:
				return GSL_SUCCESS;
			}
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			double energy(const double y[]) {
				const double r = y[1], r2 = gsl_pow_2(r);
				const double sint2 = gsl_pow_2(sin(y[2])), cost = sign(y[2]) * cos(y[2]);
				const double Delta = r2 - 2. * m * r - l2 + a2;
				return (Delta - a2 * sint2 + 2. * ((m * r + l2) * a * sint2 + Delta * l * cost) * y[7]) / ((r2 + gsl_pow_2(l + a * cost)) * y[4]);
			}
			double energyHamiltonian(const double y[]) {
				return 1. - y[4];
			}
			double angularMomentum(const double y[]) {
				const double r = y[1], r2 = gsl_pow_2(r);
				const double sint2 = gsl_pow_2(sin(y[2])), cost = sign(y[2]) * cos(y[2]);
				const double Delta = r2 - 2. * m * r - l2 + a2;
				return (-2. * ((m * r + l2) * a * sint2 + Delta * l * cost) + (gsl_pow_2(r2 + l2 + a2) * sint2 - gsl_pow_2(a * sint2 - 2. * l * cost) * Delta) * y[7]) / ((r2 + gsl_pow_2(l + a * cost)) * y[4]);
			}
			double angularMomentumHamiltonian(const double y[]) {
				return y[7];
			}
			double carter(const double y[], const double mu2) { //TODO:
				const double r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2])), cost2 = gsl_pow_2(cos(y[2]));
				const double rho2 = r2 + a2 * cost2;
				const double mr_rho2 = 2. * m * y[1] / rho2;
				return mu2 * cost2 * a2 + (gsl_pow_2(rho2 * y[6]) + cost2 * (gsl_pow_2(-mr_rho2 * a + (a2 + r2 + mr_rho2 * a2 * sint2) * y[7]) * sint2 - a2 * gsl_pow_2(mr_rho2 * (1. - a * sint2 * y[7]) - 1.))) / gsl_pow_2(y[4]);
			}
			double carterHamiltonian(const double y[], const double mu2) { //TODO:
				return gsl_pow_2(y[6]) + gsl_pow_2(cos(y[2])) * (a2 * (mu2 - gsl_pow_2(1. - y[4])) + gsl_pow_2(y[7] / sin(y[2])));
			}
			int particleNormalization(double y[]) {
				const double r = y[1], r2 = gsl_pow_2(r);
				const double sint2 = gsl_pow_2(sin(y[2])), cost = sign(y[2]) * cos(y[2]);
				const double Delta = r2 - 2. * m * r - l2 + a2;
				const double rho2 = r2 + gsl_pow_2(l + a * cost);
				y[7] += 2. * m * a * r / (gsl_pow_2(a2 + r2) - a2 * Delta * sint2);
				y[4] = sqrt(((Delta - a2 * sint2) + 2. * ((m * r + l2) * a * sint2 + Delta * l * cost) * y[7] - (gsl_pow_2(r2 + l2 + a2) * sint2 - gsl_pow_2(a * sint2 - 2. * l * cost) * Delta) * gsl_pow_2(y[7])) / rho2 - rho2 * (gsl_pow_2(y[5]) / Delta + gsl_pow_2(y[6])));
				return std::isnan(y[4]);
			}
			int lightNormalization(double y[], double e) {
				const double r = y[1], r2 = gsl_pow_2(r);
				const double sint = abs(sin(y[2])), sint2 = gsl_pow_2(sint), cost = sign(y[2]) * cos(y[2]);
				const double Delta = r2 - 2. * m * r - l2 + a2;
				const double lacost = l + a * cost;
				const double rho2 = r2 + gsl_pow_2(lacost), rho_2 = 1. / rho2;
				const double chi = a * sint2 - 2. * l * cost, rho2achi = r2 + l2 + a2;
				const double effa = rho2 / Delta * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + rho_2 * (gsl_pow_2(rho2achi) * sint2 - gsl_pow_2(chi) * Delta) * gsl_pow_2(y[7]);
				const double effb = 2. * rho_2 * (Delta * chi - a * rho2achi * sint2) * y[7];
				const double effc = -rho_2 * (Delta - a2 * sint2);
				const double eff = 0.5 * (-effb + sqrt(gsl_pow_2(effb) - 4. * effa * effc)) / effa;
				y[4] = e;
				y[5] *= eff;
				y[6] *= eff;
				y[7] *= eff;
				return 0;
			}
		} // namespace KerrTaubNUT
	}	  // namespace Metric
} // namespace SBody
