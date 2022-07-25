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
	namespace Object {
		class object {
		  public:
			/**
			 * @brief
			 *
			 * @param current
			 * @param last
			 * @return int
			 */
			virtual int hit(const double current[], const double last[]) = 0;
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
			virtual double frequency(const double photon[]) = 0;
		};
		extern std::vector<object *> objectList;
		class star : public object {
		  protected:
			///
			const int fixed;
			const double radius;
			/// square of radius
			const double r2;

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
			star(double radius, const double position[], int fixed = 0);
			/**
			 * @brief Check if the star hit by the photon
			 *
			 * @param current
			 * @param last
			 * @return int
			 */
			int hit(const double current[], const double last[] = nullptr);
			/**
			 * @brief return frequency
			 *
			 * @param photon
			 * @return double
			 */
			double frequency(const double photon[]);
			/**
			 * @brief return frequency
			 *
			 * @param photon
			 * @return double
			 */
			double frequencyTau(const double photon[]);
		};
		class disk : public object {
		  protected:
			const double innerRadius;
			const double outerRadius;

		  public:
			disk(double innerRadius, double outerRadius);
			int hit(const double current[], const double last[] = nullptr);
			double frequency(const double ph[]);
		};
		class thickDisk : public disk {
		  protected:
			const double halfAngle;

		  public:
			thickDisk(double innerRadius, double outerRadius, double halfAngle);
			int hit(const double current[], const double last[] = nullptr);
			double frequency(const double ph[]);
		};
		class torus : public object {
		  protected:
			const double majorRadius;
			const double minorRadius;

		  public:
			torus(double majorRadius, double minorRadius);
			int hit(const double current[], const double last[] = nullptr);
			double frequency(const double ph[]);
		};
	} // namespace Object
} // namespace SBody

#endif
