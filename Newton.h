#ifndef _NEWTON_H
#define _NEWTON_H

namespace SBody {
	namespace Metric {
		namespace Newton {
			extern const int dimension;
			int c2s(const double x[], double r[]);
			int s2c(const double r[], double x[]);
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double x[], void *params);
			double angularMomentum(const double x[], void *params);
		} // namespace Newton
	}	  // namespace Metric
} // namespace SBody

#endif
