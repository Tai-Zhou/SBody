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

#include <fmt/format.h>
#include <gsl/gsl_matrix.h>

#include "Utility.h"

namespace SBody {
	enum time_system { T,
					   TAU };

	/**
	 * @brief
	 *
	 */
	enum coordinate_system { LAGRANGIAN,
							 HAMILTONIAN };

	enum motion_mode { GEODESIC,
					   CIRCULAR,
					   RIAF,
					   HELICAL };
	template <typename T>
	auto format_as(T e) {
		return fmt::underlying(e);
	}
	/**
	 * @brief
	 *
	 */
	class Metric {
	  public:
		virtual std::string Name() = 0;
		virtual int MetricTensor(const double position[], gsl_matrix *metric) = 0;
		virtual double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) = 0;
		virtual double DistanceSquare(const double x[], const double y[], const size_t dimension) = 0;
		int LocalInertialFrame(const double position[], gsl_matrix *coordinate, const double timelike[] = nullptr);
		virtual int LagrangianToHamiltonian(double y[]) = 0;
		virtual int HamiltonianToLagrangian(double y[]) = 0;
		virtual int FastTrace(const double observer_r, const double observer_theta, const double sin_theta_observer, const double cos_theta_observer, const double target_r, const double target_theta, const double target_phi, double &alpha, double &beta, double *photon) = 0;
		virtual double Energy(const double y[], time_system time, coordinate_system coordinate) = 0;
		virtual double AngularMomentum(const double y[], time_system time, coordinate_system coordinate) = 0;
		virtual double CarterConstant(const double y[], const double mu2, time_system time, coordinate_system coordinate) = 0;
		virtual double Redshift(const double y[], const double photon[], time_system time = T); // TODO: add coordinate
		virtual int NormalizeTimelikeGeodesic(double y[]) = 0;
		virtual int NormalizeNullGeodesic(double y[], double frequency = 1.) = 0;
		virtual std::unique_ptr<Integrator> GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) = 0;
	};
	class Newton : public Metric {
	  public:
		const int PN_;
		Newton(int PN);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		int FastTrace(const double observer_r, const double observer_theta, const double sin_theta_observer, const double cos_theta_observer, const double target_r, const double target_theta, const double target_phi, double &alpha, double &beta, double *photon) override;
		double Energy(const double y[], time_system time, coordinate_system coordinate) override;
		double AngularMomentum(const double y[], time_system time, coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, time_system time, coordinate_system coordinate) override;
		double Redshift(const double y[], const double photon[], time_system time = T) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class PN1 : public Newton {
	  public:
		const double PN1_;
		PN1(double fSP);
		std::unique_ptr<Integrator> GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class Schwarzschild : public Metric {
	  public:
		Schwarzschild();
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		int FastTrace(const double observer_r, const double observer_theta, const double sin_theta_observer, const double cos_theta_observer, const double target_r, const double target_theta, const double target_phi, double &alpha, double &beta, double *photon) override;
		double Energy(const double y[], time_system time, coordinate_system coordinate) override;
		double AngularMomentum(const double y[], time_system time, coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, time_system time, coordinate_system coordinate) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class ReissnerNordstrom : public Schwarzschild {
	  public:
		const double r_Q_, r_Q2_, r_Q4_;
		ReissnerNordstrom(double charge);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		int FastTrace(const double observer_r, const double observer_theta, const double sin_theta_observer, const double cos_theta_observer, const double target_r, const double target_theta, const double target_phi, double &alpha, double &beta, double *photon) override;
		double Energy(const double y[], time_system time, coordinate_system coordinate) override;
		double AngularMomentum(const double y[], time_system time, coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, time_system time, coordinate_system coordinate) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class Kerr : public Schwarzschild {
	  public:
		const double a_, a2_, a4_;
		Kerr(double spin);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		int FastTrace(const double observer_r, const double observer_theta, const double sin_theta_observer, const double cos_theta_observer, const double target_r, const double target_theta, const double target_phi, double &alpha, double &beta, double *photon) override;
		double Energy(const double y[], time_system time, coordinate_system coordinate) override;
		double AngularMomentum(const double y[], time_system time, coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, time_system time, coordinate_system coordinate) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class KerrNewman : public Kerr {
	  public:
		const double r_Q_, r_Q2_, r_Q4_;
		KerrNewman(double charge);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		int FastTrace(const double observer_r, const double observer_theta, const double sin_theta_observer, const double cos_theta_observer, const double target_r, const double target_theta, const double target_phi, double &alpha, double &beta, double *photon) override;
		double Energy(const double y[], time_system time, coordinate_system coordinate) override;
		double AngularMomentum(const double y[], time_system time, coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, time_system time, coordinate_system coordinate) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class KerrTaubNUT : public Kerr {
	  public:
		const double l_, l2_, l4_;
		KerrTaubNUT(double spin, double NUT);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		int FastTrace(const double observer_r, const double observer_theta, const double sin_theta_observer, const double cos_theta_observer, const double target_r, const double target_theta, const double target_phi, double &alpha, double &beta, double *photon) override;
		double Energy(const double y[], time_system time, coordinate_system coordinate) override;
		double AngularMomentum(const double y[], time_system time, coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, time_system time, coordinate_system coordinate) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class Hayward : public Kerr {
	  public:
		const double alpha_, beta_, g_;
		Hayward(double alpha, double beta, double charge);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		int FastTrace(const double observer_r, const double observer_theta, const double sin_theta_observer, const double cos_theta_observer, const double target_r, const double target_theta, const double target_phi, double &alpha, double &beta, double *photon) override;
		double Energy(const double y[], time_system time, coordinate_system coordinate) override;
		double AngularMomentum(const double y[], time_system time, coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, time_system time, coordinate_system coordinate) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	// Post-Newtonian
	int NewtonTLagrangianGeodesic(double t, const double y[], double dydt[], void *params);

	// Post-Newtonian-1
	int PN1TLagrangianGeodesic(double t, const double y[], double dydt[], void *params);

	// Schwarzschild
	int SchwarzschildTLagrangianGeodesic(double t, const double y[], double dydt[], void *params);
	int SchwarzschildTLagrangianCircular(double t, const double y[], double dydt[], void *params);
	int SchwarzschildTLagrangianRIAF(double t, const double y[], double dydt[], void *params);
	int SchwarzschildTLagrangianHelical(double t, const double y[], double dydt[], void *params);
	int SchwarzschildTHamiltonianGeodesic(double t, const double y[], double dydt[], void *params);
	int SchwarzschildTauLagrangianGeodesic(double t, const double y[], double dydt[], void *params);

	// Kerr
	int KerrTLagrangianGeodesic(double t, const double y[], double dydt[], void *params);
	int KerrTLagrangianHelical(double t, const double y[], double dydt[], void *params);
	int KerrTHamiltonianGeodesic(double t, const double y[], double dydt[], void *params);
	int KerrTauLagrangianGeodesic(double t, const double y[], double dydt[], void *params);
	int KerrTauHamiltonianGeodesic(double t, const double y[], double dydt[], void *params);

	// Kerr-Taub-NUT
	int KerrTaubNutTLagrangianGeodesic(double t, const double y[], double dydt[], void *params);
	int KerrTaubNutTLagrangianHelical(double t, const double y[], double dydt[], void *params);
	int KerrTaubNutTHamiltonianGeodesic(double t, const double y[], double dydt[], void *params);
	int KerrTaubNutTauLagrangianGeodesic(double t, const double y[], double dydt[], void *params);

	int Jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
} // namespace SBody

#endif
