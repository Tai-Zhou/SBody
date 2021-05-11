#ifndef _KERRH_H
#define _KERRH_H

namespace SBody {
	namespace Metric {
		namespace KerrH {
			extern const int dimension;
			int c2s(const double x[], double r[]);
			int s2c(const double r[], double x[]);
			int qdq2qp(const double r[], double u[], void *params);
			int qp2qdq(const double u[], double r[], void *params);
			int function(double t, const double y[], double dydt[], void *params);
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[], void *params);
			double angularMomentum(const double y[], void *params);
			double carter(const double r[], void *params);
			namespace particle {
				int normalization(double y[], void *params);
			} // namespace particle
			namespace light {
				int normalization(double y[], void *params);
			} // namespace light
		}	  // namespace KerrH
	}		  // namespace Metric
} // namespace SBody

#endif
