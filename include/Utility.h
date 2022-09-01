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
	/// Global absolute accuracy
	extern double absolute_accuracy;

	/// Global relative accuracy
	extern double relative_accuracy;

	/// Sample number
	constexpr int sample_number = 100;

	/// Epsilon, \f$\epsilon\f$
	constexpr double epsilon = 1e-10;

	/// \f$\sin\epsilon\f$
	extern const double sin_epsilon;

	/// \f$\cos\epsilon\f$
	extern const double cos_epsilon;

	/// \f$2\pi\f$
	constexpr double M_2PI = 6.28318530717958647692528676655900576;

	/// \f$\pi^2\f$
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
		 * @brief Construct a new Integrator object
		 *
		 * @param function The function calculates \f$dy/dt\f$
		 * @param jacobian The jacobian of the function \f$J_{ij}=\partial f_i/\partial y_j\f$
		 * @param coordinate The coordinate of the system, `0` for Cartesian, `1` for spherical, and `2` for modified spherical
		 * @param params The parameters passed to the function
		 * @param type Type of the algorithms
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
	 * @brief A wrapper of the `gsl_vector`, `gsl_matrix`, and `gsl_permutation`
	 *
	 */
	class GslBlock {
	  private:
		std::vector<gsl_vector *> vectors_;
		std::vector<gsl_matrix *> matrices_;
		std::vector<gsl_permutation *> permutations_;

	  public:
		/// Destroy the GslBlock object, free the memory
		~GslBlock();

		/**
		 * @brief :func:`gsl_vector_alloc`, creates a vector of length `n`
		 *
		 * @param n length of the vector
		 * @return gsl_vector*
		 */
		gsl_vector *VectorAlloc(size_t n);

		/**
		 * @brief :func:`gsl_vector_calloc`, creates a vector of length `n` and initializes all the elements of the vector to zero
		 *
		 * @param n length of the vector
		 * @return gsl_vector*
		 */
		gsl_vector *VectorCalloc(size_t n);

		/**
		 * @brief :func:`gsl_matrix_alloc`, creates a matrix of size `n1` rows by `n2` columns
		 *
		 * @param n1 rows of the matrix
		 * @param n2 columns of the matrix
		 * @return gsl_matrix*
		 */
		gsl_matrix *MatrixAlloc(size_t n1, size_t n2);

		/**
		 * @brief :func:`gsl_matrix_alloc`, creates a matrix of size `n1` rows by `n2` columns and initializes all the elements of the matrix to zero
		 *
		 * @param n1 rows of the matrix
		 * @param n2 columns of the matrix
		 * @return gsl_matrix*
		 */
		gsl_matrix *MatrixCalloc(size_t n1, size_t n2);

		/**
		 * @brief :func:`gsl_permutation_alloc`, allocates memory for a new permutation of size `n`
		 *
		 * @param n size of the permutation
		 * @return gsl_permutation*
		 */
		gsl_permutation *PermutationAlloc(size_t n);

		/**
		 * @brief :func:`gsl_permutation_alloc`, allocates memory for a new permutation of size `n` and initializes it to the identity
		 *
		 * @param n size of the permutation
		 * @return gsl_permutation*
		 */
		gsl_permutation *PermutationCalloc(size_t n);
	};

	/**
	 * @brief Dot product of vector \f$x\cdot y\f$, or \f$x\cdot x\f$ if :code:`y == nullptr`
	 *
	 * @param x
	 * @param y
	 * @param dimension
	 * @return double
	 */
	double Dot(const double x[], const double y[] = nullptr, size_t dimension = 3);

	/**
	 * @brief Euclidean norm of `x`, \f$\|x\|_2\f$
	 *
	 * @param x vector
	 * @param dimension the dimension of the vector
	 * @return double
	 */
	double Norm(const double x[], size_t dimension = 3);

	/**
	 * @brief Cross product of vector \f$x\times y\f$, stored in z
	 *
	 * @param x 3 dimensional vector
	 * @param y 3 dimensional vector
	 * @param z 3 dimensional vector
	 */
	int Cross(const double x[], const double y[], double z[]);

	/**
	 * @brief Rotate vector `x` around the `axis` by `angle`
	 *
	 * @param x 3 dimensional vector
	 * @param axis the subscript of the rotation axis
	 * @param angle in unit rad
	 * @return int
	 */
	int RotateAroundAxis(double x[], int axis, double angle);

	/**
	 * @brief
	 *
	 * @param x 4 or 8 dimensional vector
	 * @param calculate_velocity if `true`, the velocity (`x[5]` - `x[8]`) is also calculated
	 * @return int
	 */
	int CartesianToSpherical(double x[], bool calculate_velocity = true);

	/**
	 * @brief
	 *
	 * @param cartesian 3 dimensional vector
	 * @param spherical 3 dimensional vector
	 * @return int
	 */
	int CartesianToSpherical(const double cartesian[], double spherical[]);

	/**
	 * @brief
	 *
	 * @param cartesian_position 3 dimensional vector
	 * @param cartesian_velocity 3 dimensional vector
	 * @param spherical_position 3 dimensional vector
	 * @param spherical_velocity 3 dimensional vector
	 * @return int
	 */
	int CartesianToSpherical(const double cartesian_position[], const double cartesian_velocity[], double spherical_position[], double spherical_velocity[]);

	/**
	 * @brief
	 *
	 * @param x 4 or 8 dimensional vector
	 * @param calculate_velocity if `true`, the velocity (`x[5]` - `x[8]`) is also calculated
	 * @return int
	 */
	int SphericalToCartesian(double x[], bool calculate_velocity = true);

	/**
	 * @brief
	 *
	 * @param spherical 3 dimensional vector
	 * @param cartesian 3 dimensional vector
	 * @return int
	 */
	int SphericalToCartesian(const double spherical[], double cartesian[]);

	/**
	 * @brief
	 *
	 * @param spherical_position 3 dimensional vector
	 * @param spherical_velocity 3 dimensional vector
	 * @param cartesian_position 3 dimensional vector or `nullptr`
	 * @param cartesian_velocity 3 dimensional vector
	 * @return int
	 */
	int SphericalToCartesian(const double spherical_position[], const double spherical_velocity[], double cartesian_position[], double cartesian_velocity[]);

	/**
	 * @brief return `1` if `x`, `y` have opposite signs, else `0`.
	 *
	 * @param x
	 * @param y
	 * @return int
	 */
	int OppositeSign(double x, double y);

	/**
	 * @brief return `x` in \f$[0, 2\pi)\f$.
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
