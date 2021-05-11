#include "Kerr.h"

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include "Constant.h"
#include "Metric.h"
#include "Utility.h"

namespace SBody {
	namespace Metric {
		namespace Kerr {
			const int dimension = 8;
			// from cartesian to spherical
			int c2s(const double x[], double r[]) {
				r[0] = x[0];
				r[4] = x[4];
				return Metric::c2s(x + 1, x + 5, r + 1, r + 5);
			}
			// from spherical to cartesian
			int s2c(const double r[], double x[]) {
				x[0] = r[0];
				x[4] = r[4];
				return Metric::s2c(r + 1, r + 5, x + 1, x + 5);
			}
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
			double carter(const double r[], void *params) {
				const double m = ((source *)params)->mass, a = ((source *)params)->spin;
				const double a2 = gsl_pow_2(a), r2 = gsl_pow_2(r[1]), sint2 = gsl_pow_2(sin(r[2])), cost2 = gsl_pow_2(cos(r[2]));
				const double rho2 = r2 + a2 * cost2;
				const double mr_rho2 = 2. * m * r[1] / rho2;
				return cost2 * a2 + (gsl_pow_2(rho2 * r[6]) + cost2 * (gsl_pow_2(-mr_rho2 * a + (a2 + r2 + mr_rho2 * a2 * sint2) * r[7]) * sint2 - a2 * gsl_pow_2(mr_rho2 * (1. - a * sint2 * r[7]) - 1.))) / gsl_pow_2(r[4]);
			}
			namespace particle {
				int normalization(double y[], void *params) {
					const double m = ((source *)params)->mass, a = ((source *)params)->spin, r = y[1], sint = sin(y[2]);
					const double a2 = gsl_pow_2(a), r2 = gsl_pow_2(r), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint);
					const double rho2 = r2 + a2 * gsl_pow_2(cos(y[2]));
					const double mr_rho2 = 2. * m * r / rho2;
					y[4] = sqrt(1. - mr_rho2 + 2. * mr_rho2 * a * sint2 * y[7] - (rho2 / (r2 - 2. * m * r + a2) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2 + r2) * sint2 + mr_rho2 * a2 * sint4) * gsl_pow_2(y[7])));
					return std::isnan(y[4]);
				}
			} // namespace particle
			namespace light {
				int normalization(double y[], void *params) {
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
			} // namespace light
		}	  // namespace Kerr
	}		  // namespace Metric
} // namespace SBody
