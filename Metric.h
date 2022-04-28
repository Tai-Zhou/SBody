#ifndef _METRIC_H
#define _METRIC_H

#include <string>

namespace SBody {
	namespace Metric {
		extern double m;
		extern double a, a2, a4;
		extern double e, e2, e4;
		extern double l, l2, l4;
		extern std::string name;
		extern double (*dot)(const double[], const double[], const double[], const size_t);
		extern double (*ds2)(const double[], const double[], const size_t);
		extern int (*qdq2qp)(double[]);
		extern int (*qp2qdq)(double[]);
		extern int (*function)(double, const double[], double[], void *);
		extern int (*functionTau)(double, const double[], double[], void *);
		extern int (*functionHamiltonian)(double, const double[], double[], void *);
		extern int (*jacobian)(double, const double[], double *, double[], void *);
		extern int (*jacobianHamiltonian)(double, const double[], double *, double[], void *);
		extern double (*energy)(const double[]);
		extern double (*energyHamiltonian)(const double[]);
		extern double (*angularMomentum)(const double[]);
		extern double (*angularMomentumHamiltonian)(const double[]);
		extern double (*carter)(const double[], const double);
		extern double (*carterHamiltonian)(const double[], const double);
		extern int (*particleNormalization)(double[]);
		extern int (*lightNormalization)(double[], double);

		// set function pointers above
		void setMetric(int NSK, int PN, double mass, double spin, double charge, double NUT);

		// from cartesian to spherical
		int c2s(const double x[], const double v[], double r[], double w[]);
		int c2s(const double x[], double r[]);

		// from spherical to cartesian
		int s2c(const double r[], const double w[], double x[], double v[]);
		int s2c(const double r[], double x[]);

		namespace Newton {
			extern int PN;
			double dot(const double g[], const double x[], const double y[], const size_t dimension);
			double ds2(const double x[], const double y[], const size_t dimension);
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[]);
			double angularMomentum(const double y[]);
			double carter(const double y[], const double mu2);
			int particleNormalization(double y[]);
			int lightNormalization(double y[], double e);
		} // namespace Newton
		namespace Schwarzschild {
			double dot(const double g[], const double x[], const double y[], const size_t dimension);
			double ds2(const double x[], const double y[], const size_t dimension);
			int qdq2qp(double y[]);
			int qp2qdq(double y[]);
			int function(double t, const double y[], double dydt[], void *params);
			int functionTau(double t, const double y[], double dydt[], void *params);
			int functionHamiltonian(double t, const double y[], double dydt[], void *params);
			int functionRIAF(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianRIAF(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[]);
			double energyHamiltonian(const double y[]);
			double angularMomentum(const double y[]);
			double angularMomentumHamiltonian(const double y[]);
			double carter(const double y[], const double mu2);
			double carterHamiltonian(const double y[], const double mu2);
			int particleNormalization(double y[]);
			int lightNormalization(double y[], double e);
		} // namespace Schwarzschild
		namespace Kerr {
			double dot(const double g[], const double x[], const double y[], const size_t dimension);
			double ds2(const double x[], const double y[], const size_t dimension);
			int qdq2qp(double y[]);
			int qp2qdq(double y[]);
			int function(double t, const double y[], double dydt[], void *params);
			int functionTau(double t, const double y[], double dydt[], void *params);
			int functionHamiltonian(double t, const double y[], double dydt[], void *params);
			int functionHamiltonianTau(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[]);
			double energyHamiltonian(const double y[]);
			double angularMomentum(const double y[]);
			double angularMomentumHamiltonian(const double y[]);
			double carter(const double y[], const double mu2);
			double carterHamiltonian(const double y[], const double mu2);
			int particleNormalization(double y[]);
			int lightNormalization(double y[], double e);
		} // namespace Kerr
		namespace KerrTaubNUT {
			double dot(const double g[], const double x[], const double y[], const size_t dimension);
			double ds2(const double x[], const double y[], const size_t dimension);
			int qdq2qp(double y[]);
			int qp2qdq(double y[]);
			int function(double t, const double y[], double dydt[], void *params);
			int functionTau(double t, const double y[], double dydt[], void *params);
			int functionHamiltonian(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[]);
			double energyHamiltonian(const double y[]);
			double angularMomentum(const double y[]);
			double angularMomentumHamiltonian(const double y[]);
			double carter(const double y[], const double mu2);
			double carterHamiltonian(const double y[], const double mu2);
			int particleNormalization(double y[]);
			int lightNormalization(double y[], double e);
		} // namespace KerrTaubNUT
	}	  // namespace Metric
} // namespace SBody

#endif
