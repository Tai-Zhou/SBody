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
	enum coordinate_system { BASE,
							 LAGRANGIAN,
							 HAMILTONIAN };

	enum motion_mode { GEODESIC,
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
		virtual int BaseToHamiltonian(double y[]) = 0;
		virtual int HamiltonianToBase(double y[]) = 0;
		virtual double Energy(const double y[], coordinate_system coordinate) = 0;
		virtual double AngularMomentum(const double y[], coordinate_system coordinate) = 0;
		virtual double CarterConstant(const double y[], const double mu2, coordinate_system coordinate) = 0;
		virtual double Redshift(const double y[], const double photon[], time_system time = T);
		virtual int NormalizeTimelikeGeodesic(double y[]) = 0;
		virtual int NormalizeNullGeodesic(double y[], double frequency = 1.) = 0;
		virtual Integrator GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) = 0;
	};
	class Newton : public Metric {
	  public:
		const int PN_;
		Newton(int PN);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int BaseToHamiltonian(double y[]) override;
		int HamiltonianToBase(double y[]) override;
		double Energy(const double y[], coordinate_system coordinate) override;
		double AngularMomentum(const double y[], coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, coordinate_system coordinate) override;
		double Redshift(const double y[], const double photon[], time_system time = T) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		Integrator GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class PN1 : public Newton {
	  public:
		const double PN1_;
		PN1(double fSP);
		Integrator GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class Schwarzschild : public Metric {
	  public:
		Schwarzschild();
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int BaseToHamiltonian(double y[]) override;
		int HamiltonianToBase(double y[]) override;
		double Energy(const double y[], coordinate_system coordinate) override;
		double AngularMomentum(const double y[], coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, coordinate_system coordinate) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		Integrator GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class Kerr : public Schwarzschild {
	  public:
		const double a_, a2_, a4_;
		Kerr(double spin);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int BaseToHamiltonian(double y[]) override;
		int HamiltonianToBase(double y[]) override;
		double Energy(const double y[], coordinate_system coordinate) override;
		double AngularMomentum(const double y[], coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, coordinate_system coordinate) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		Integrator GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	class KerrTaubNUT : public Kerr {
	  public:
		const double l_, l2_, l4_;
		KerrTaubNUT(double spin, double NUT);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int BaseToHamiltonian(double y[]) override;
		int HamiltonianToBase(double y[]) override;
		double Energy(const double y[], coordinate_system coordinate) override;
		double AngularMomentum(const double y[], coordinate_system coordinate) override;
		double CarterConstant(const double y[], const double mu2, coordinate_system coordinate) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		Integrator GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion = GEODESIC) override;
	};
	namespace metric {
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
