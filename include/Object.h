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

#include <vector>

namespace SBody {
	/**
	 * @brief
	 *
	 */
	class Object {
	  public:
		static std::vector<Object *> object_list_;
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
		/// if the position of star fixed
		const bool fixed_;
		/// radius of star
		const double radius_;
		/// square of radius
		const double radius_square_;

	  public:
		/// 8 dimensional information of the star
		double pos[8]; // FIXME
		/**
		 * @brief Construct a new star object
		 *
		 * @param radius radius
		 * @param position initial position and velocity
		 * @param fixed whether the position of the star is fixed
		 */
		Star(double radius, const double position[], bool fixed = false);
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

	  public:
		/**
		 * @brief Construct a new Disk object
		 *
		 * @param inner_radius
		 * @param outer_radius
		 */
		Disk(double inner_radius, double outer_radius);
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
		 * @param inner_radius
		 * @param outer_radius
		 * @param half_angle
		 */
		ThickDisk(double inner_radius, double outer_radius, double half_angle);
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
		 * @param major_radius
		 * @param minor_radius
		 */
		Torus(double major_radius, double minor_radius);
		int Hit(const double current[], const double last[] = nullptr);
		double Redshift(const double ph[]);
	};
} // namespace SBody

#endif
