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

#include <gsl/gsl_matrix.h>

#include "Utility.h"

namespace SBody {
	/**
	 * @brief
	 *
	 */
	enum metric_mode { T,
					   TAU,
					   HAMILTONIAN };

	/**
	 * @brief
	 *
	 */
	class Metric {
	  protected:
		int mode_;
		std::string name_;

	  public:
		Metric(metric_mode mode, std::string name);
		std::string &Name();
		virtual int MetricTensor(const double position[], gsl_matrix *metric);
		virtual double DotProduct(const double position[], const double x[], const double y[], const size_t dimension);
		virtual double Distance(const double x[], const double y[], const size_t dimension);
		int LocalInertialFrame(const double position[], gsl_matrix *coordinate, const double timelike[] = nullptr);
		virtual int LagrangianToHamiltonian(double y[]);
		virtual int HamiltonianToLagrangian(double y[]);
		virtual double Energy(const double y[]);
		virtual double AngularMomentum(const double y[]);
		virtual double CarterConstant(const double y[], const double mu2);
		virtual int NormalizeTimelikeGeodesic(double y[]);
		virtual int NormalizeNullGeodesic(double y[], double frequency = 1.);
		virtual Integrator GetIntegrator(int coordinate);
	};
	class Newton : public Metric {
	  public:
		const int PN_;
		Newton(int PN, metric_mode mode, std::string name = "Newton");
		int MetricTensor(const double position[], gsl_matrix *metric);
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension);
		double Distance(const double x[], const double y[], const size_t dimension);
		int LagrangianToHamiltonian(double y[]);
		int HamiltonianToLagrangian(double y[]);
		double Energy(const double y[]);
		double AngularMomentum(const double y[]);
		double CarterConstant(const double y[], const double mu2);
		int NormalizeTimelikeGeodesic(double y[]);
		int NormalizeNullGeodesic(double y[], double frequency = 1.);
		Integrator GetIntegrator(int coordinate);
	};
	class Schwarzschild : public Metric {
	  public:
		Schwarzschild(metric_mode mode, std::string name = "Schwarzschild");
		int MetricTensor(const double position[], gsl_matrix *metric);
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension);
		double Distance(const double x[], const double y[], const size_t dimension);
		int LagrangianToHamiltonian(double y[]);
		int HamiltonianToLagrangian(double y[]);
		double Energy(const double y[]);
		double AngularMomentum(const double y[]);
		double CarterConstant(const double y[], const double mu2);
		int NormalizeTimelikeGeodesic(double y[]);
		int NormalizeNullGeodesic(double y[], double frequency = 1.);
		Integrator GetIntegrator(int coordinate);
	};
	class Kerr : public Schwarzschild {
	  public:
		const double a_, a2_, a4_;
		Kerr(double spin, metric_mode mode, std::string name = "Kerr");
		int MetricTensor(const double position[], gsl_matrix *metric);
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension);
		double Distance(const double x[], const double y[], const size_t dimension);
		int LagrangianToHamiltonian(double y[]);
		int HamiltonianToLagrangian(double y[]);
		double Energy(const double y[]);
		double AngularMomentum(const double y[]);
		double CarterConstant(const double y[], const double mu2);
		int NormalizeTimelikeGeodesic(double y[]);
		int NormalizeNullGeodesic(double y[], double frequency = 1.);
		Integrator GetIntegrator(int coordinate);
	};
	class KerrTaubNUT : public Kerr {
	  public:
		const double e_, e2_, e4_;
		const double l_, l2_, l4_;
		KerrTaubNUT(double spin, double charge, double NUT, metric_mode mode, std::string name = "Kerr-Taub-NUT");
		int MetricTensor(const double position[], gsl_matrix *metric);
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension);
		double Distance(const double x[], const double y[], const size_t dimension);
		int LagrangianToHamiltonian(double y[]);
		int HamiltonianToLagrangian(double y[]);
		double Energy(const double y[]);
		double AngularMomentum(const double y[]);
		double CarterConstant(const double y[], const double mu2);
		int NormalizeTimelikeGeodesic(double y[]);
		int NormalizeNullGeodesic(double y[], double frequency = 1.);
		Integrator GetIntegrator(int coordinate);
	};
	namespace metric {
		// from cartesian to spherical
		int c2s(const double x[], const double v[], double r[], double w[]);
		int c2s(const double x[], double r[]);

		// from spherical to cartesian
		int s2c(const double r[], const double w[], double x[], double v[]);
		int s2c(const double r[], double x[]);

		namespace Newton {
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		} // namespace Newton
		namespace Schwarzschild {
			int function(double t, const double y[], double dydt[], void *params);
			int functionTau(double t, const double y[], double dydt[], void *params);
			int functionHamiltonian(double t, const double y[], double dydt[], void *params);
			int functionRIAF(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianTau(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianRIAF(double t, const double y[], double *dfdy, double dfdt[], void *params);
		} // namespace Schwarzschild
		namespace Kerr {
			int function(double t, const double y[], double dydt[], void *params);
			int functionTau(double t, const double y[], double dydt[], void *params);
			int functionHamiltonian(double t, const double y[], double dydt[], void *params);
			int functionHamiltonianTau(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianTau(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		} // namespace Kerr
		namespace KerrTaubNUT {
			int function(double t, const double y[], double dydt[], void *params);
			int functionTau(double t, const double y[], double dydt[], void *params);
			int functionHamiltonian(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianTau(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		} // namespace KerrTaubNUT
	}	  // namespace metric
} // namespace SBody

#endif
