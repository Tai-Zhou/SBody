#ifndef _METRIC_H
#define _METRIC_H

#include <string>

namespace SBody {
	struct source {
		const double mass;
		const double spin;
		source(double mass, double spin = 0);
	};

	namespace Metric {
		extern std::string name;
		extern int (*function)(double, const double[], double[], void *);
		extern int (*jacobian)(double, const double[], double *, double[], void *);
		extern double (*energy)(const double[], void *);
		extern double (*angularMomentum)(const double[], void *);
		extern double (*carter)(const double[], void *params);

		// set function pointers above
		void setMetric(int NSK);

		// from cartesian to spherical
		int c2s(const double x[], const double v[], double r[], double w[]);
		int c2s(const double x[], double r[]);

		// from spherical to cartesian
		int s2c(const double r[], const double w[], double x[], double v[]);
		int s2c(const double r[], double x[]);

		namespace Newton {
			extern int PN;
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double x[], void *params);
			double angularMomentum(const double x[], void *params);
			double carter(const double y[], void *params);
		} // namespace Newton
		namespace Schwarzschild {
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[], void *params);
			double angularMomentum(const double y[], void *params);
			double carter(const double y[], void *params);
			int particleNormalization(double y[], void *params);
			int lightNormalization(double y[], void *params);
		} // namespace Schwarzschild
		namespace Kerr {
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[], void *params);
			double angularMomentum(const double y[], void *params);
			double carter(const double y[], void *params);
			int particleNormalization(double y[], void *params);
			int lightNormalization(double y[], void *params);
		} // namespace Kerr
		namespace KerrH {
			int qdq2qp(const double r[], double u[], void *params);
			int qp2qdq(const double u[], double r[], void *params);
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[], void *params);
			double angularMomentum(const double y[], void *params);
			double carter(const double y[], void *params);
			int particleNormalization(double y[], void *params);
			int lightNormalization(double y[], void *params);
		} // namespace KerrH
	}	  // namespace Metric
} // namespace SBody

#endif
