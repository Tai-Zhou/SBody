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
					   HAMILTONIAN,
					   RIAF,
					   HELICAL };

	/**
	 * @brief
	 *
	 */
	class Metric {
	  protected:
		int mode_;

	  public:
		Metric(metric_mode mode);
		virtual std::string Name() = 0;
		virtual int MetricTensor(const double position[], gsl_matrix *metric) = 0;
		virtual double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) = 0;
		virtual double Distance(const double x[], const double y[], const size_t dimension) = 0;
		int LocalInertialFrame(const double position[], gsl_matrix *coordinate, const double timelike[] = nullptr);
		virtual int LagrangianToHamiltonian(double y[]) = 0;
		virtual int HamiltonianToLagrangian(double y[]) = 0;
		virtual double Energy(const double y[]) = 0;
		virtual double AngularMomentum(const double y[]) = 0;
		virtual double CarterConstant(const double y[], const double mu2) = 0;
		virtual int NormalizeTimelikeGeodesic(double y[]) = 0;
		virtual int NormalizeNullGeodesic(double y[], double frequency = 1.) = 0;
		virtual Integrator GetIntegrator() = 0;
	};
	class Newton : public Metric {
	  private:
		const int PN_;

	  public:
		Newton(int PN, metric_mode mode);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double Distance(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		double Energy(const double y[]) override;
		double AngularMomentum(const double y[]) override;
		double CarterConstant(const double y[], const double mu2) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		Integrator GetIntegrator() override;
	};
	class PN1 : public Newton {
	  private:
		const double PN1_;

	  public:
		PN1(double fSP, metric_mode mode);
		Integrator GetIntegrator() override;
	};
	class Schwarzschild : public Metric {
	  public:
		Schwarzschild(metric_mode mode);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double Distance(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		double Energy(const double y[]) override;
		double AngularMomentum(const double y[]) override;
		double CarterConstant(const double y[], const double mu2) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		Integrator GetIntegrator() override;
	};
	class Kerr : public Schwarzschild {
	  public:
		const double a_, a2_, a4_;
		Kerr(double spin, metric_mode mode);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double Distance(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		double Energy(const double y[]) override;
		double AngularMomentum(const double y[]) override;
		double CarterConstant(const double y[], const double mu2) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		Integrator GetIntegrator() override;
	};
	class KerrTaubNUT : public Kerr {
	  public:
		const double e_, e2_, e4_;
		const double l_, l2_, l4_;
		KerrTaubNUT(double spin, double charge, double NUT, metric_mode mode);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double Distance(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		double Energy(const double y[]) override;
		double AngularMomentum(const double y[]) override;
		double CarterConstant(const double y[], const double mu2) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		Integrator GetIntegrator() override;
	};
	namespace metric {
		namespace Newton {
			int function(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
		} // namespace Newton
		namespace Schwarzschild {
			int function(double t, const double y[], double dydt[], void *params);
			int functionTau(double t, const double y[], double dydt[], void *params);
			int functionHamiltonian(double t, const double y[], double dydt[], void *params);
			int functionRIAF(double t, const double y[], double dydt[], void *params);
			int functionHelicalWithFixedRadialSpeed(double t, const double y[], double dydt[], void *params);
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianTau(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianRIAF(double t, const double y[], double *dfdy, double dfdt[], void *params);
			int jacobianHelicalWithFixedRadialSpeed(double t, const double y[], double *dfdy, double dfdt[], void *params);
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
