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
	class Integrator {
	  private:
		///
		const int coordinate_;
		///
		gsl_odeiv2_control *control_;
		///
		gsl_odeiv2_evolve *evolve_;
		///
		gsl_odeiv2_step *step_;
		///
		gsl_odeiv2_system system_;

	  public:
		/**
		 * @brief Construct a new Integrator object
		 *
		 * @param function
		 * @param jacobian
		 * @param coordinate
		 * @param params
		 * @param type
		 */
		Integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), int coordinate, void *params = nullptr, const gsl_odeiv2_step_type *type = gsl_odeiv2_step_rk8pd);

		/// Destroy the Integrator object
		~Integrator();

		/**
		 * @brief
		 *
		 * @param t
		 * @param t1
		 * @param h
		 * @param y
		 * @return int
		 */
		int Apply(double *t, double t1, double *h, double *y);

		/**
		 * @brief
		 *
		 * @param t
		 * @param h
		 * @param y
		 * @return int
		 */
		int ApplyFixedStep(double *t, const double h, double *y);

		/**
		 * @brief
		 *
		 * @return int
		 */
		int Reset();

		/**
		 * @brief
		 *
		 * @param y
		 * @return int
		 */
		int CheckCoordinate(double *y);
	};

	/**
	 * @brief Dot product of vector x·y, or x·x if y == nullptr
	 *
	 * @param x
	 * @param y
	 * @param dimension
	 * @return double
	 */
	double Dot(const double x[], const double y[] = nullptr, size_t dimension = 3);

	/**
	 * @brief Length of vector x, with 3 dimensions set by default
	 *
	 * @param x
	 * @param dimension
	 * @return double
	 */
	double Norm(const double x[], size_t dimension = 3);

	/**
	 * @brief Cross product of vector \f$x\timesy\f$, stored in z
	 *
	 * @param x
	 * @param y
	 * @param z
	 */
	void Cross(const double x[], const double y[], double z[]);

	/**
	 * @brief return 1 if x, y have opposite signs
	 *
	 * @param x
	 * @param y
	 * @return int
	 */
	int OppositeSign(double x, double y);

	/**
	 * @brief return x in [0, 2*pi)
	 *
	 * @param x
	 * @return double
	 */
	double ModBy2Pi(double x);

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
