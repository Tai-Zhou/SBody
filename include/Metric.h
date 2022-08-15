/**
 * @file Metric.h
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_METRIC_H
#define SBODY_METRIC_H

#include <string>

#include "Utility.h"

namespace SBody {
	/**
	 * @brief
	 *
	 */
	class Metric {
	  protected:
		std::string name_;

	  public:
		Metric(std::string name);
		double (*dot)(const double[], const double[], const double[], const size_t);
		double (*ds2)(const double[], const double[], const size_t);
		int (*qdq2qp)(double[]);
		int (*qp2qdq)(double[]);
		int (*function)(double, const double[], double[], void *);
		int (*functionTau)(double, const double[], double[], void *);
		int (*functionHamiltonian)(double, const double[], double[], void *);
		int (*jacobian)(double, const double[], double *, double[], void *);
		int (*jacobianHamiltonian)(double, const double[], double *, double[], void *);
		double (*energy)(const double[]);
		double (*energyHamiltonian)(const double[]);
		double (*angularMomentum)(const double[]);
		double (*angularMomentumHamiltonian)(const double[]);
		double (*carter)(const double[], const double);
		double (*carterHamiltonian)(const double[], const double);
		int (*particleNormalization)(double[]);
		int (*lightNormalization)(double[], double);
		void SetMetric(int metric, int PN, double mass, double spin, double charge, double NUT);
		Integrator SetIntegrator(int coordinate, void *params = nullptr);
	};
	class Newton : public Metric {
	};
	class Schwarzschild : public Metric {
	};
	class Kerr : public Schwarzschild {
	  public:
		double a_, a2_, a4_;
	};
	class KerrTaubNUT : public Kerr {
	  public:
		double e_, e2_, e4_;
		double l_, l2_, l4_;
	};
	namespace metric {
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
		void setMetric(int metric, int PN, double mass, double spin, double charge, double NUT);

		// from cartesian to spherical
		int c2s(const double x[], const double v[], double r[], double w[]);
		int c2s(const double x[], double r[]);

		// from spherical to cartesian
		int s2c(const double r[], const double w[], double x[], double v[]);
		int s2c(const double r[], double x[]);

		namespace Newton {
			extern int PN;
			/**
			 * @brief This calculate the dot product of 2 vectors.
			 *
			 * @param g the tetrad position in the spacetime.
			 * @param x vector x.
			 * @param y vector y.
			 * @param dimension 3 if calculate space part only; 4 if the time added.
			 * @return double
			 */
			double dot(const double g[], const double x[], const double y[], const size_t dimension);
			double ds2(const double x[], const double y[], const size_t dimension);
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			double energy(const double y[]);
			double AngularMomentum(const double y[]);
			double carter(const double y[], const double mu2);
			int particleNormalization(double y[]);
			int lightNormalization(double y[], double e);
		} // namespace Newton
		namespace Schwarzschild {
			int gmunu(const double pos[], double gmunu[]);
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
			int LocalInertialFrame(const double y[], double coordinate[]);
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
	}	  // namespace metric
} // namespace SBody

#endif
