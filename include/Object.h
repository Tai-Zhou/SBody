/**
 * @file Object.h
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_OBJECT_H
#define SBODY_OBJECT_H

#include <memory>
#include <vector>

#include <gsl/gsl_matrix.h>

#include "Metric.h"

namespace SBody {
	/**
	 * @brief
	 *
	 */
	class Object {
	  protected:
		std::shared_ptr<Metric> metric_;

	  public:
		static std::vector<Object *> object_list_;
		/**
		 * @brief Construct a new Object object
		 *
		 * @param metric
		 */
		Object(std::shared_ptr<Metric> metric);
		/**
		 * @brief
		 *
		 * @param current
		 * @param last
		 * @return int
		 */
		virtual int Hit(const double current[], const double last[]) = 0;
		/**
		 * @brief
		 * \f[z=\frac{E_\mathrm{obj}}{E_\mathrm{obs}}-1\f]
		 * where
		 * \f[E_\mathrm{obj}=-g_{\mu\nu}u_\mathrm{ph}^\mu u_\mathrm{obj}^\nu\f]
		 * For the distant observer, the equation is the same. With the tetrad velocity \f$u_\mathrm{obs}=(1,0,0,0)\f$ and metric
		 * \f[
		 * 	g_\mathrm{obs}=
		 * 	\begin{pmatrix}
		 * 		-1 & 0 & 0 & 0 \\
		 * 		0  & 1 & 0 & 0 \\
		 * 		0  & 0 & 1 & 0 \\
		 * 		0  & 0 & 0 & 1
		 * 	\end{pmatrix},
		 * \f]
		 * the energy of the photon can be simplified to
		 * \f[
		 * 	E_\mathrm{obs}=u_\mathrm{ph}^0.
		 * \f]
		 * @param photon 8 dimensional information of photon
		 * @return double
		 */
		virtual double Redshift(const double photon[], TimeSystem photon_time) = 0;
	};

	/**
	 * @brief
	 *
	 */
	class Particle : public Object {
	  protected:
		TimeSystem time_;
		DynamicalSystem coordinate_;
		/// 8 dimensional information
		double position_[8];
		/// if the position fixed
		const bool fixed_;
		/// integrator
		std::unique_ptr<Integrator> integrator_;

	  public:
		/**
		 * @brief Construct a new star object
		 *
		 * @param metric
		 * @param radius radius
		 * @param fixed whether the position of the star is fixed
		 */
		Particle(std::shared_ptr<Metric> metric, TimeSystem time, DynamicalSystem coordinate, bool fixed = false);
		int Position(double *position);
		int InitializeKeplerian(double a, double e, double inclination, double periapsis, double ascending_node, double true_anomaly, double observer_inclination = 0., double observer_rotation = 0.);
		int InitializeGeodesic(double orbital_radius, double inclination, double periapsis, double ascending_node, double v_r, double v_phi, double observer_inclination = 0., double observer_rotation = 0.);
		int InitializeSchwarzschildKeplerianPericenter(double a, double e, double inclination, double periapsis, double ascending_node, double observer_inclination, double observer_rotation);
		int InitializeSchwarzschildKeplerianApocenter(double a, double e, double inclination, double periapsis, double ascending_node, double observer_inclination, double observer_rotation);
		int InitializeKeplerianHarmonic(double a, double e, double inclination, double periapsis, double ascending_node, double true_anomaly, double observer_inclination, double observer_rotation);
		int InitializeCircular(double r, double phi, double v_phi_ratio);
		int InitializeHelical(double r, double theta, double phi, double v_r, double v_phi);
		int IntegratorApply(double *t, double t1, double *h);
		int IntegratorApplyFixedStep(double *t, const double h);
		int IntegratorReset();

		/**
		 * @brief Check if the star hit by the photon
		 *
		 * @param current
		 * @param last
		 * @return int
		 */
		int Hit(const double current[], const double last[] = nullptr) override;

		/**
		 * @brief return frequency
		 *
		 * @param photon
		 * @return double
		 */
		double Redshift(const double photon[], TimeSystem photon_time) override;

		int MetricTensor(gsl_matrix *metric);
		double DotProduct(const double x[], const double y[], const size_t dimension);
		int LocalInertialFrame(gsl_matrix *coordinate);
		double Energy();
		double AngularMomentum();
		double CarterConstant();
	};

	class Star : public Particle {
	  protected:
		/// radius
		const double radius_;
		/// radius square
		const double radius_square_;

	  public:
		Star(std::shared_ptr<Metric> metric, TimeSystem time, DynamicalSystem coordinate, double radius = 0, bool fixed = false);
		int Hit(const double current[], const double last[]);
	};

	class HotSpot : public Particle {
	  protected:
		const double spectral_index_;
		const double luminosity_;
		const double luminosity_mu_;
		const double luminosity_sigma_;

	  public:
		HotSpot(std::shared_ptr<Metric> metric, TimeSystem time, DynamicalSystem coordinate, double spectral_index, double luminosity, double luminosity_mu, double luminosity_sigma, bool fixed = false);
		double Luminosity(double t);
		double SpectralDensity(double t, double redshift);
	};

	class Disk : public Object {
	  protected:
		const double inner_radius;
		const double outer_radius;
		const double inner_velocity;
		const double velocity_slope_index;

	  public:
		int Hit(const double current[], const double last[]);
	};
} // namespace SBody

#endif
