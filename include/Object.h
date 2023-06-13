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
		virtual double Redshift(const double photon[]) = 0;
	};
	/**
	 * @brief
	 *
	 */
	class Star : public Object {
	  protected:
		/// 8 dimensional information
		double position_[8];
		/// if the position fixed
		const bool fixed_;
		/// radius
		const double radius_;
		/// radius square
		const double radius_square_;
		/// integrator
		Integrator integrator_;

	  public:
		/**
		 * @brief Construct a new star object
		 *
		 * @param metric
		 * @param radius radius
		 * @param fixed whether the position of the star is fixed
		 */
		Star(std::shared_ptr<Metric> metric, double radius = 0, bool fixed = false);
		int Position(double *position);
		int InitializeKeplerian(double a, double e, double inclination, double periapsis, double ascending_node, double true_anomaly, double observer_inclination = 0., double observer_rotation = 0.);
		int InitializeGeodesic(double orbital_radius, double inclination, double periapsis, double ascending_node, double v_r, double v_phi, double observer_inclination = 0., double observer_rotation = 0.);
		int InitializeSchwarzschildKeplerian(double a, double e, double inclination, double periapsis, double ascending_node, double observer_inclination, double observer_rotation);
		int InitializeHelical(double r, double theta, double phi, double v_r, double v_phi);
		int IntegratorApply(double *t, double t1, double *h);
		int IntegratorApplyFixedStep(double *t, const double h);

		/**
		 * @brief Check if the star hit by the photon
		 *
		 * @param current
		 * @param last
		 * @return int
		 */
		int Hit(const double current[], const double last[] = nullptr);

		/**
		 * @brief return frequency
		 *
		 * @param photon
		 * @return double
		 */
		double Redshift(const double photon[]);

		/**
		 * @brief return frequency
		 *
		 * @param photon
		 * @return double
		 */
		double RedshiftTau(const double photon[]);
		int MetricTensor(gsl_matrix *metric);
		double DotProduct(const double x[], const double y[], const size_t dimension);
		int LocalInertialFrame(gsl_matrix *coordinate);
		double Energy();
		double AngularMomentum();
		double CarterConstant();
	};
	/**
	 * @brief
	 *
	 */
	class Disk : public Object {
	  protected:
		/// inner radius
		const double inner_radius_;
		/// outer radius
		const double outer_radius_;
		/// emissivity index
		const double emissivity_index_;

	  public:
		/**
		 * @brief Construct a new Disk object
		 *
		 * @param metric
		 * @param inner_radius
		 * @param outer_radius
		 */
		Disk(std::shared_ptr<Metric> metric, double inner_radius, double outer_radius, double emissivity_index = 0.);
		int Hit(const double current[], const double last[] = nullptr);
		double Redshift(const double ph[]);
	};
	/**
	 * @brief
	 *
	 */
	class ThickDisk : public Disk {
	  protected:
		const double half_angle_;

	  public:
		/**
		 * @brief Construct a new Thick Disk object
		 *
		 * @param metric
		 * @param inner_radius
		 * @param outer_radius
		 * @param half_angle
		 */
		ThickDisk(std::shared_ptr<Metric> metric, double inner_radius, double outer_radius, double half_angle);
		int Hit(const double current[], const double last[] = nullptr);
		double Redshift(const double ph[]);
	};
	/**
	 * @brief
	 *
	 */
	class Torus : public Object {
	  protected:
		const double major_radius_;
		const double minor_radius_;

	  public:
		/**
		 * @brief Construct a new Torus object
		 *
		 * @param metric
		 * @param major_radius
		 * @param minor_radius
		 */
		Torus(std::shared_ptr<Metric> metric, double major_radius, double minor_radius);
		int Hit(const double current[], const double last[] = nullptr);
		double Redshift(const double ph[]);
	};
} // namespace SBody

#endif
