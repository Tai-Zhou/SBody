#ifndef _METRIC_H
#define _METRIC_H

#include <string>

namespace SBody {
	namespace Metric {
		extern double m;
		extern double a, a2, a4;
		extern std::string name;
		extern double (*ds2)(const double[], const double[]);
		extern int (*function)(double, const double[], double[], void *);
		extern int (*jacobian)(double, const double[], double *, double[], void *);
		extern double (*energy)(const double[]);
		extern double (*angularMomentum)(const double[]);
		extern double (*carter)(const double[]);
		extern int (*particleNormalization)(double[]);
		extern int (*lightNormalization)(double[]);

		// set function pointers above
		void setMetric(int NSK, double mass, double spin);

		// from cartesian to spherical
		int c2s(const double x[], const double v[], double r[], double w[]);
		int c2s(const double x[], double r[]);

		// from spherical to cartesian
		int s2c(const double r[], const double w[], double x[], double v[]);
		int s2c(const double r[], double x[]);

		namespace Newton {
			extern int PN;
			double ds2(const double x[], const double y[]);
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double x[]);
			double angularMomentum(const double x[]);
			double carter(const double y[]);
			int particleNormalization(double y[]);
			int lightNormalization(double y[]);
		} // namespace Newton
		namespace Schwarzschild {
			double ds2(const double x[], const double y[]);
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[]);
			double angularMomentum(const double y[]);
			double carter(const double y[]);
			int particleNormalization(double y[]);
			int lightNormalization(double y[]);
		} // namespace Schwarzschild
		namespace Kerr {
			double ds2(const double x[], const double y[]);
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[]);
			double angularMomentum(const double y[]);
			double carter(const double y[]);
			int particleNormalization(double y[]);
			int lightNormalization(double y[]);
		} // namespace Kerr
		namespace KerrH {
			int qdq2qp(double r[]);
			int qp2qdq(double r[]);
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[]);
			double angularMomentum(const double y[]);
			double carter(const double y[]);
		} // namespace KerrH
	}	  // namespace Metric
} // namespace SBody

#endif
