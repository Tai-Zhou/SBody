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

#include <memory>
#include <string>

#include <gsl/gsl_matrix.h>

#include "Utility.h"

namespace SBody {
	/**
	 * @brief The time system used in the record of the object.
	 *
	 * - T: the time step of the integrator is the observer's time \f$t\f$, and the first element of the object vector is its intrinsic time \f$\tau\f$.
	 * - TAU: the time step of the integrator is the intrinsic time of the object, and the first element is the observer's time.
	 *
	 */
	enum TimeSystem { T,
					  TAU };

	/**
	 * @brief The dynamical system to describe the motion of the object.
	 *
	 * - LAGRANGIAN: the velocity part of the object vector is \f$\mathrm d{x^\mu}/\mathrm dt\f$ or \f$\mathrm d{x^\mu}/\mathrm d\tau\f$.
	 * - HAMILTONIAN: the velocity is recorded as \f$p_\mu\f$.
	 *
	 */
	enum DynamicalSystem { LAGRANGIAN,
						   HAMILTONIAN };

	/**
	 * @brief The mode of the motion of the object.
	 *
	 * - GEODESIC: the motion is described by the geodesic equation.
	 * - CIRCULAR: the orbit has fixed radius.
	 * - HELICAL: the motion is described by the helical orbit equation.
	 *
	 */
	enum MotionMode { GEODESIC,
					  CIRCULAR,
					  HELICAL };

	/**
	 * @brief The base class of all kinds of metrics.
	 *
	 */
	class Metric {
	  public:
		/**
		 * @brief Return the name of the metric.
		 *
		 * @return name of the metric.
		 */
		virtual std::string Name() = 0;

		/**
		 * @brief Calculate the metric tensor at `position`, stored in `metric`.
		 *
		 * @param position 4 dimensional vector
		 * @param metric matrix with size 4×4
		 * @return status
		 */
		virtual int MetricTensor(const double position[], gsl_matrix *metric) = 0;

		/**
		 * @brief Dot product of vector `x` and `y` at `position`. \f$g_{\mu\nu}x^\mu y^\nu\f$
		 *
		 * @param position 4 dimensional vector, position to calcuate the dot product of `x` and `y`.
		 * @param x 4 dimensional vector
		 * @param y 4 dimensional vector
		 * @param dimension dimension of the vector, should be 3 or 4.
		 * @return result
		 */
		virtual double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) = 0;

		/**
		 * @brief Calculate the square of the distance between `x` and `y` at `x`. \f$g_{\mu\nu}(x^\mu-y^\mu)(x^\nu-y^\nu)\f$
		 *
		 * @param x 4 dimensional vector
		 * @param y 4 dimensional vector
		 * @param dimension dimension of the vector
		 * @return result
		 */
		virtual double DistanceSquare(const double x[], const double y[], const size_t dimension) = 0;

		/**
		 * @brief Calculate the local inertial frame coordinate of the object at `position`, stored in `coordinate`.
		 *
		 * @param position 8 dimensional vector
		 * @param time time systme of the object
		 * @param coordinate matrix with size 4×4
		 * @return status
		 */
		int LocalInertialFrame(const double position[], TimeSystem time, gsl_matrix *coordinate);

		/**
		 * @brief Convert the coordinate system from Lagrangian to Hamiltonian.
		 *
		 * @param y 8 dimensional vector
		 * @return status
		 */
		virtual int LagrangianToHamiltonian(double y[]) = 0;

		/**
		 * @brief Convert the coordinate system from Hamiltonian to Lagrangian.
		 *
		 * @param y 8 dimensional vector
		 * @return status
		 */
		virtual int HamiltonianToLagrangian(double y[]) = 0;

		int InitializePhoton(double photon[], double alpha, double beta, double r, double r2, double theta, double sin_theta);

		int AngularMomentumCarterConstantToAlphaBeta(double l, double q2, double cos_theta, double sin_theta, double &alpha, double &abs_beta);

		/**
		 * @brief Trace the photon from the observer to the target, using elliptic integrals.
		 *
		 * @param r_observer radius of the observer, \f$r_\text{obs}\f$.
		 * @param theta_observer theta of the observer, \f$\theta_\text{obs}\f$.
		 * @param sin_theta_observer sin(theta_observer), \f$\sin\theta_\text{obs}\f$.
		 * @param cos_theta_observer cos(theta_observer), \f$\cos\theta_\text{obs}\f$.
		 * @param r_object radius of the target, \f$r_\text{tar}\f$.
		 * @param theta_object theta of the target, \f$\theta_\text{tar}\f$.
		 * @param phi_object phi of the target, \f$\phi_\text{tar}\f$.
		 * @param alpha x position of the target in the observer's view.
		 * @param beta y position of the target in the observer's view.
		 * @param photon 9 dimensional vector, position and the velocity of the photon traced to the target. photon[8] is used to store the look back time.
		 * @return status
		 */
		virtual int FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_object, const double theta_object, const double phi_object, double &alpha, double &beta, double photon[]) = 0;

		virtual int FastShadow(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double alpha, const double beta, const double r_min);
		/**
		 * @brief Calculate the energy of the object.
		 *
		 * @param y 8 dimensional vector
		 * @param time time system of the object
		 * @param dynamics dynamical system of the object
		 * @return result
		 */
		virtual double Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) = 0;
		/**
		 * @brief Calculate the angular momentum of the object.
		 *
		 * @param y 8 dimensional vector
		 * @param time time system of the object
		 * @param dynamics dynamical system of the object
		 * @return result
		 */
		virtual double AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) = 0;

		/**
		 * @brief Calculate the Carter constant of the object.
		 *
		 * @param y 8 dimensional vector
		 * @param time time system of the object
		 * @param dynamics dynamical system of the object
		 * @return result
		 */
		virtual double CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) = 0;

		/**
		 * @brief Calculate the redshift of the object, \f$1+z\f$
		 *
		 * @param y 8 dimensional vector of the object
		 * @param photon photon traced to the object
		 * @param object_time time system of the object
		 * @param photon_time time system of the photon
		 * @return result
		 */
		virtual double Redshift(const double y[], const double photon[], TimeSystem object_time, TimeSystem photon_time);

		/**
		 * @brief Normalize the timelike geodesic.
		 *
		 * @param y 8 dimensional vector
		 * @return status
		 */
		virtual int NormalizeTimelikeGeodesic(double y[]) = 0;

		/**
		 * @brief Normalize the null geodesic.
		 *
		 * @param y 8 dimensional vector
		 * @return status
		 */
		virtual int NormalizeNullGeodesic(double y[], double frequency = 1.) = 0;

		/**
		 * @brief Get the integrator to calculate the motion of the object.
		 *
		 * @param time time system of the object
		 * @param dynamics dynamical system of the object
		 * @param motion motion mode of the object
		 * @return pointer to the integrator
		 */
		virtual std::unique_ptr<Integrator> GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) = 0;
	};

	/// Post-Newtonian
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
		int FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_object, const double theta_object, const double phi_object, double &alpha, double &beta, double photon[]) override;
		double Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) override;
		double Redshift(const double y[], const double photon[], TimeSystem object_time, TimeSystem photon_time) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override;
	};
	class PN1 : public Newton {
	  public:
		const double PN1_;
		PN1(double fSP);
		std::unique_ptr<Integrator> GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override;
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
		int FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_object, const double theta_object, const double phi_object, double &alpha, double &beta, double photon[]) override;
		int FastShadow(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double alpha, const double beta, const double r_min) override;
		double Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override;
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
		int FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_object, const double theta_object, const double phi_object, double &alpha, double &beta, double photon[]) override;
		double Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override;
		static int LagrangianGeodesic(double t, const double y[], double dydt[], void *params);
	};
	class Kerr : public Schwarzschild {
	  protected:
		static int DeltaUMuPhi(const gsl_vector *alpha_beta, void *params, gsl_vector *delta_u_mu_phi);
		static int UIntegral(double a, double a2, double u_plus_1, double u_minus_1, double u_plus, double u_minus, double l, double l2, double q2, double r_obs, double u_obs, double u_obj, double &I_u_0, double &I_u_1, double &I_u_plus_0, double &I_u_plus_1, double &I_u_minus_0, double &I_u_minus_1, double &I_u_2_0, double &I_u_2_1, double &I_u_4_0, double &I_u_4_1);
		static int Mu0Integral(double a, double q2, double M_plus, double M_minus, double delta_M, double delta_M_plus, int beta_sign, double mu_obs, double &mu_plus, double &mu_minus, double &A, double &k, double &n, int &alpha_1, double &I_mu_0, double &I_mu_full_turn, double &I_t_mu_0, double &I_t_mu_full_turn, double &I_phi_mu_0, double &I_phi_mu_full_turn);
		static int MuFIntegral(double a, double a2, double l, double q2, double M_plus, double M_minus, double delta_M, double delta_M_plus, double mu_plus, double mu_minus, double I_u, double I_mu_0, double I_mu_full_turn, double I_t_mu_0, double I_t_mu_full_turn, double I_phi_mu_0, double I_phi_mu_full_turn, double A, double k, double n, int alpha_1, int &alpha_2, double &mu_f, double &t_mu, double &phi_mu);

	  public:
		const double a_, a2_, a4_;
		const double sqrt_delta_a2_;
		const double u_plus_1, u_plus, u_minus_1, u_minus, u_r;
		Kerr(double spin);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		int FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_object, const double theta_object, const double phi_object, double &alpha, double &beta, double photon[]) override;
		int CalcThetaPhi(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, double alpha, double beta, const std::vector<double> &u, double theta_0[], double theta_1[], double phi_0[], double phi_1[]);
		int FastShadow(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double alpha, const double beta, const double r_min) override;
		double Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override;
	};

	struct KerrFastTraceParameters {
		Kerr *const kerr;
		const double r, r2;
		const double u_obs, u_obj;
		const double mu_obs, mu_obj;
		const double theta_obs;
		const double sin_theta_obs, sin_theta_obj;
		const double phi_obj;
		double E, L, Q;
		double t;
		double tau;
		int u_dir, mu_dir;
		KerrFastTraceParameters(Kerr *kerr, double r, double r2, double u_obs, double u_obj, double mu_obs, double mu_obj, double theta_obs, double sin_theta_obs, double sin_theta_obj, double phi_obj);
	};

	class KerrNewman : public Kerr {
	  public:
		const double r_Q_, r_Q2_, r_Q4_;
		KerrNewman(double spin, double charge);
		std::string Name() override;
		int MetricTensor(const double position[], gsl_matrix *metric) override;
		double DotProduct(const double position[], const double x[], const double y[], const size_t dimension) override;
		double DistanceSquare(const double x[], const double y[], const size_t dimension) override;
		int LagrangianToHamiltonian(double y[]) override;
		int HamiltonianToLagrangian(double y[]) override;
		int FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_object, const double theta_object, const double phi_object, double &alpha, double &beta, double photon[]) override;
		double Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override;
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
		int FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_object, const double theta_object, const double phi_object, double &alpha, double &beta, double photon[]) override;
		double Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override;
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
		int FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_object, const double theta_object, const double phi_object, double &alpha, double &beta, double photon[]) override;
		double Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) override;
		double CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) override;
		int NormalizeTimelikeGeodesic(double y[]) override;
		int NormalizeNullGeodesic(double y[], double frequency = 1.) override;
		std::unique_ptr<Integrator> GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override;
	};
	// Post-Newtonian
	int NewtonTLagrangianGeodesic(double t, const double y[], double dydt[], void *params);

	// PN1
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
	int KerrTLagrangianCircular(double t, const double y[], double dydt[], void *params);
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
