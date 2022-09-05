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
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector.h>

namespace SBody {
	/// Global absolute accuracy.
	extern double absolute_accuracy;

	/// Global relative accuracy.
	extern double relative_accuracy;

	/// Sample number.
	constexpr int sample_number = 100;

	/// Epsilon, a small value. \f$\varepsilon\f$.
	constexpr double epsilon = 1e-10;

	/// \f$\sin\varepsilon\f$.
	constexpr double sin_epsilon = 1e-10;

	/// \f$\cos\varepsilon\f$.
	constexpr double cos_epsilon = 0.999999999999999999995;

	/// \f$2\pi\f$.
	constexpr double M_2PI = 6.28318530717958647692528676655900576;

	/// \f$\pi^2\f$.
	constexpr double M_PI2 = 9.86960440108935861883449099987615111;

	/**
	 * @brief A wrapper of the gsl_odeiv2_evolve
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
		 * @brief Construct a new Integrator object.
		 *
		 * @param function The function calculates \f[\frac{\mathrm dy_i(t)}{\mathrm dt}=f_i\left[t,y_1(t),\dots,y_n(t)\right].\f]
		 * @param jacobian The jacobian of the function \f[J_{ij}=\frac{\partial f_i\left[t,y(t)\right]}{\partial y_j}.\f]
		 * @param coordinate The coordinate of the system, `0` for Cartesian, `1` for spherical, and `2` for modified spherical. The modified spherical coordiante maps the polar angle from \f$\theta\in[0,\pi]\f$ to \f[
		 * \theta'\equiv\begin{cases}\theta & (\theta\leq\pi/2)\\
		 * \theta-\pi & (\theta>\pi/2)\end{cases}.
		 * \f]
		 * @param params The parameters passed to the function, like the PN parameter, or the spin of the black hole.
		 * @param type Type of the algorithms. Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method is set by default.
		 */
		Integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), int coordinate, void *params = nullptr, const gsl_odeiv2_step_type *type = gsl_odeiv2_step_rk8pd);

		/// Destructor
		~Integrator();

		/**
		 * @brief `gsl_odeiv2_evolve_apply`.
		 *
		 * @param t time.
		 * @param t1 maximum time not to be exceeded by the time step.
		 * @param h step size.
		 * @param y values of the integrated system.
		 * @return status
		 */
		int Apply(double *t, double t1, double *h, double *y);

		/**
		 * @brief `gsl_odeiv2_evolve_apply_fixed_step`.
		 *
		 * @param t time.
		 * @param h step size.
		 * @param y values of the integrated system.
		 * @return status
		 */
		int ApplyFixedStep(double *t, const double h, double *y);

		/**
		 * @brief Resets the evolution function and the stepping function.
		 *
		 * @return status
		 */
		int Reset();

		/**
		 * @brief
		 *
		 * @param y
		 * @return status
		 */
		int CheckCoordinate(double *y);
	};

	/**
	 * @brief A wrapper of the `gsl_vector`, `gsl_matrix`, and `gsl_permutation`.
	 *
	 */
	class GslBlock {
	  private:
		std::vector<gsl_vector *> vectors_;
		std::vector<gsl_matrix *> matrices_;
		std::vector<gsl_permutation *> permutations_;

	  public:
		/// Destructor
		~GslBlock();

		/**
		 * @brief `gsl_vector_alloc`, creates a vector of length `n`
		 *
		 * @param n length of the vector
		 * @return pointer
		 */
		gsl_vector *VectorAlloc(size_t n);

		/**
		 * @brief `gsl_vector_calloc`, creates a vector of length `n` and initializes all the elements of the vector to zero.
		 *
		 * @param n length of the vector
		 * @return pointer
		 */
		gsl_vector *VectorCalloc(size_t n);

		/**
		 * @brief `gsl_vector_alloc_from_block`, creates a vector with data from `block`.
		 *
		 * @param block the data block
		 * @param offset the offset to the data block
		 * @param n length of the vector
		 * @param stride step-size of the vector
		 * @return pointer
		 */
		gsl_vector *VectorAllocFromBlock(gsl_block *block, const size_t offset, const size_t n, const size_t stride = 1);

		/**
		 * @brief `gsl_matrix_alloc`, creates a matrix of size `n1` rows by `n2` columns.
		 *
		 * @param n1 rows of the matrix
		 * @param n2 columns of the matrix
		 * @return pointer
		 */
		gsl_matrix *MatrixAlloc(size_t n1, size_t n2);

		/**
		 * @brief `gsl_matrix_alloc`, creates a matrix of size `n1` rows by `n2` columns and initializes all the elements of the matrix to zero
		 *
		 * @param n1 rows of the matrix
		 * @param n2 columns of the matrix
		 * @return pointer
		 */
		gsl_matrix *MatrixCalloc(size_t n1, size_t n2);

		/**
		 * @brief `gsl_permutation_alloc`, allocates memory for a new permutation of size `n`
		 *
		 * @param n size of the permutation
		 * @return pointer
		 */
		gsl_permutation *PermutationAlloc(size_t n);

		/**
		 * @brief `gsl_permutation_alloc`, allocates memory for a new permutation of size `n` and initializes it to the identity
		 *
		 * @param n size of the permutation
		 * @return pointer
		 */
		gsl_permutation *PermutationCalloc(size_t n);
	};

	/**
	 * @brief Dot product of vector `x` and `y`, or `x` and `x` if `y == nullptr`. \f$\vec{x}\cdot\vec{y}\f$ or \f$\vec{x}\cdot\vec{x}\f$
	 *
	 * @param x vector
	 * @param y vector
	 * @param dimension dimension of the vector
	 * @return result
	 */
	double Dot(const double x[], const double y[] = nullptr, size_t dimension = 3);

	/**
	 * @brief Euclidean norm of `x`. \f$\|\vec{x}\|_2\f$.
	 *
	 * @param x vector
	 * @param dimension dimension of the vector
	 * @return result
	 */
	double Norm(const double x[], size_t dimension = 3);

	/**
	 * @brief Cross product of vector `x` and `y`, stored in `z`. \f$\vec{z}=\vec{x}\times\vec{y}\f$.
	 *
	 * @param x 3 dimensional vector
	 * @param y 3 dimensional vector
	 * @param z 3 dimensional vector
	 * @return status
	 */
	int Cross(const double x[], const double y[], double z[]);

	/**
	 * @brief Dot product of vector `x` and the cross product of vector `y` and vector `z`. \f$\vec{x}\cdot(\vec{y}\times\vec{z})\f$.
	 *
	 * @param x 3 dimensional vector
	 * @param y 3 dimensional vector
	 * @param z 3 dimensional vector
	 * @return result
	 */
	double DotCross(const double x[], const double y[], const double z[]);

	/**
	 * @brief Rotate vector `x` around the `axis` by `angle`.
	 *
	 * @param x 3 dimensional vector
	 * @param axis the subscript of the rotation axis, should be \f$\geq0\f$ and \f$\leq2\f$.
	 * @param angle in rad
	 * @return status
	 */
	int RotateAroundAxis(double x[], int axis, double angle);

	/**
	 * @brief
	 *
	 * @param x 4 or 8 dimensional vector
	 * @param dimension 3, 4, or 8
	 * @return status
	 */
	int CartesianToSpherical(double x[], size_t dimension = 8);

	/**
	 * @brief
	 *
	 * @param cartesian 3 dimensional vector
	 * @param spherical 3 dimensional vector
	 * @return status
	 */
	int CartesianToSpherical(const double cartesian[], double spherical[]);

	/**
	 * @brief
	 *
	 * @param cartesian_position 3 dimensional vector
	 * @param cartesian_velocity 3 dimensional vector
	 * @param spherical_position 3 dimensional vector
	 * @param spherical_velocity 3 dimensional vector
	 * @return status
	 */
	int CartesianToSpherical(const double cartesian_position[], const double cartesian_velocity[], double spherical_position[], double spherical_velocity[]);

	/**
	 * @brief
	 *
	 * @param x 4 or 8 dimensional vector
	 * @param dimension 3, 4, or 8.
	 * @return status
	 */
	int SphericalToCartesian(double x[], size_t dimension = 8);

	/**
	 * @brief
	 *
	 * @param spherical 3 dimensional vector
	 * @param cartesian 3 dimensional vector
	 * @return status
	 */
	int SphericalToCartesian(const double spherical[], double cartesian[]);

	/**
	 * @brief
	 *
	 * @param spherical_position 3 dimensional vector
	 * @param spherical_velocity 3 dimensional vector
	 * @param cartesian_position 3 dimensional vector or `nullptr`
	 * @param cartesian_velocity 3 dimensional vector
	 * @return status
	 */
	int SphericalToCartesian(const double spherical_position[], const double spherical_velocity[], double cartesian_position[], double cartesian_velocity[]);

	/**
	 * @brief return `1` if `x`, `y` have opposite signs, else `0`.
	 *
	 * @param x number
	 * @param y number
	 * @return result
	 */
	int OppositeSign(double x, double y);

	/**
	 * @brief return `x` in \f$[0, 2\pi)\f$.
	 *
	 * @param x
	 * @return result
	 */
	double ModBy2Pi(double x);

	/**
	 * @brief
	 *
	 * @param x
	 * @return result
	 */
	double _0x(double x);

	/**
	 * @brief
	 *
	 * @param x
	 * @return result
	 */
	double _0x1(double x);
} // namespace SBody

#endif
