#ifndef _METRIC_H
#define _METRIC_H

namespace SBody {
	namespace Metric {
		// from cartesian to spherical
		int c2s(const double x[], const double v[], double r[], double w[]);
		int c2s(const double x[], double r[], const int dimension);
		// from spherical to cartesian
		int s2c(const double r[], const double w[], double x[], double v[]);
		int s2c(const double r[], double x[], const int dimension);
		namespace Newton {
			extern const int dimension;
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double x[], void *params);
			double angularMomentum(const double x[], void *params);
		} // namespace Newton
		namespace Postnewtonian {
			extern const int dimension;
			extern int PN;
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[], void *params);
			double angularMomentum(const double y[], void *params);
		} // namespace Postnewtonian
		namespace Schwarzschild {
			extern const int dimension;
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[], void *params);
			double angularMomentum(const double y[], void *params);
			int particleNormalization(double y[], void *params);
			int lightNormalization(double y[], void *params);
		} // namespace Schwarzschild
		namespace Kerr {
			extern const int dimension;
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[], void *params);
			double angularMomentum(const double y[], void *params);
			double carter(const double r[], void *params);
			int particleNormalization(double y[], void *params);
			int lightNormalization(double y[], void *params);
		} // namespace Kerr
		namespace KerrH {
			extern const int dimension;
			int qdq2qp(const double r[], double u[], void *params);
			int qp2qdq(const double u[], double r[], void *params);
			int function(double t, const double y[], double dydt[], void *params);
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[], void *params);
			double angularMomentum(const double y[], void *params);
			double carter(const double r[], void *params);
			int particleNormalization(double y[], void *params);
			int lightNormalization(double y[], void *params);
		} // namespace KerrH
	}	  // namespace Metric
} // namespace SBody

#endif
