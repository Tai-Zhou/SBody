/**
 * @file Utility.h
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_UTILITY_H
#define SBODY_UTILITY_H

#include <string>

#include <gsl/gsl_odeiv2.h>

namespace SBody {
	/// Absolute accuracy
	extern double absAcc;

	/// Relative accuracy
	extern double relAcc;

	/// Epsilon
	constexpr double epsilon = 1e-10;

	/// \f$2\pi\f$
	constexpr double M_2PI = 6.28318530717958647692528676655900576;

	/// \f$\pi^2\f$
	constexpr double M_PI2 = 9.86960440108935861883449099987615111;

	/**
	 * @brief
	 *
	 */
	class integrator {
	  private:
		///
		const int coordinate;
		///
		const gsl_odeiv2_step_type *type;
		///
		gsl_odeiv2_control *control;
		///
		gsl_odeiv2_evolve *evolve;
		///
		gsl_odeiv2_step *step;
		///
		gsl_odeiv2_system system;

	  public:
		/**
		 * @brief Construct a new integrator object
		 *
		 * @param function
		 * @param jacobian
		 * @param coordinate
		 * @param params
		 * @param type
		 */
		integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), int coordinate, void *params = nullptr, const gsl_odeiv2_step_type *type = gsl_odeiv2_step_rk8pd);

		/// Destroy the integrator object
		~integrator();

		/**
		 * @brief
		 *
		 * @param t
		 * @param t1
		 * @param h
		 * @param y
		 * @return int
		 */
		int apply(double *t, double t1, double *h, double *y);

		/**
		 * @brief
		 *
		 * @param t
		 * @param h
		 * @param y
		 * @return int
		 */
		int apply_fixed(double *t, const double h, double *y);

		/**
		 * @brief
		 *
		 * @return int
		 */
		int reset();

		/**
		 * @brief
		 *
		 * @param y
		 * @return int
		 */
		int checkCoordinate(double *y);
	};

	/**
	 * @brief Dot product of vector x·y, or x·x if y == nullptr
	 *
	 * @param x
	 * @param y
	 * @param dimension
	 * @return double
	 */
	double dot(const double x[], const double y[] = nullptr, size_t dimension = 3);

	/**
	 * @brief Length of vector x, with 3 dimensions set by default
	 *
	 * @param x
	 * @param dimension
	 * @return double
	 */
	double norm(const double x[], size_t dimension = 3);

	/**
	 * @brief Cross product of vector x \times y, stored in z
	 *
	 * @param x
	 * @param y
	 * @param z
	 */
	void cross(const double x[], const double y[], double z[]);

	/**
	 * @brief return 1 if x, y have opposite signs
	 *
	 * @param x
	 * @param y
	 * @return int
	 */
	int oppositeSign(double x, double y);

	/**
	 * @brief return x in [0, 2*pi)
	 *
	 * @param x
	 * @return double
	 */
	double mod2Pi(double x);

	/**
	 * @brief
	 *
	 * @param x
	 * @return double
	 */
	double _0x(double x);

	/**
	 * @brief
	 *
	 * @param x
	 * @return double
	 */
	double _0x1(double x);
} // namespace SBody

#endif
