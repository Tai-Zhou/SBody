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
#include <gsl/gsl_mode.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_roots.h>
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

	/// Area of a circle with radius of epsilon. \f$\pi\varepsilon^2\f$
	constexpr double epsilon_circle_area = M_PI * epsilon * epsilon;

	/// \f$2\pi\f$.
	constexpr double M_2PI = 6.28318530717958647692528676655900576;

	/// \f$\pi^2\f$.
	constexpr double M_PI2 = 9.86960440108935861883449099987615111;

	/// \f$\sqrt27\f$
	constexpr double M_SQRT27 = 5.19615242270663188058233902451761710;

	enum Axis { X,
				Y,
				Z };

	/**
	 * @brief A wrapper of the gsl_odeiv2_evolve
	 *
	 */
	class Integrator {
	  private:
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
		 * @param params The parameters passed to the function, like the PN parameter, or the spin of the black hole.
		 * @param type Type of the algorithms. Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method is set by default.
		 */
		Integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), void *params = nullptr, const gsl_odeiv2_step_type *type = gsl_odeiv2_step_rk8pd);

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
		 * @brief `gsl_odeiv2_evolve_apply`.
		 *
		 * @param t time.
		 * @param t1 maximum time not to be exceeded by the time step.
		 * @param h step size.
		 * @param y values of the integrated system.
		 * @return status
		 */
		int ApplyStep(double *t, double t1, double *h, double *y);

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

	class Solver {
	  public:
		virtual int Iterate() = 0;
		virtual int Solve() = 0;
		virtual double Root() = 0;
	};

	/**
	 * @brief
	 *
	 */
	class FunctionSolver : public Solver {
	  private:
		gsl_root_fsolver *solver_;

	  public:
		FunctionSolver(const gsl_root_fsolver_type *type = gsl_root_fsolver_brent);
		FunctionSolver(gsl_function *function, double lower, double upper, const gsl_root_fsolver_type *type = gsl_root_fsolver_brent);
		~FunctionSolver();
		int Set(gsl_function *function, double lower, double upper);
		int Iterate() override;
		int Solve() override;
		double Root() override;
		double Lower();
		double Upper();
	};

	class DerivativeSolver : public Solver {
	  private:
		gsl_root_fdfsolver *solver_;

	  public:
		DerivativeSolver(const gsl_root_fdfsolver_type *type = gsl_root_fdfsolver_steffenson);
		DerivativeSolver(gsl_function_fdf *function, double root, const gsl_root_fdfsolver_type *type = gsl_root_fdfsolver_steffenson);
		~DerivativeSolver();
		int Set(gsl_function_fdf *function, double root);
		int Iterate() override;
		int Solve() override;
		double Root() override;
	};

	class MultiSolver {
	  public:
		virtual int Iterate() = 0;
		virtual int Solve() = 0;
		virtual gsl_vector *Root() = 0;
		virtual gsl_vector *Value() = 0;
		virtual gsl_vector *StepSize() = 0;
	};

	class MultiFunctionSolver : public MultiSolver {
	  private:
		gsl_multiroot_fsolver *solver_;

	  public:
		MultiFunctionSolver(size_t n, const gsl_multiroot_fsolver_type *type = gsl_multiroot_fsolver_broyden);
		MultiFunctionSolver(gsl_multiroot_function *function, const gsl_vector *x, size_t n, const gsl_multiroot_fsolver_type *type = gsl_multiroot_fsolver_broyden);
		~MultiFunctionSolver();
		int Set(gsl_multiroot_function *function, const gsl_vector *x);
		int Iterate() override;
		int Solve() override;
		gsl_vector *Root() override;
		gsl_vector *Value() override;
		gsl_vector *StepSize() override;
	};

	class MultiDerivativeSolver : public MultiSolver {
	  private:
		gsl_multiroot_fdfsolver *solver_;

	  public:
		MultiDerivativeSolver(size_t n, const gsl_multiroot_fdfsolver_type *type = gsl_multiroot_fdfsolver_gnewton);
		MultiDerivativeSolver(gsl_multiroot_function_fdf *function, const gsl_vector *x, size_t n, const gsl_multiroot_fdfsolver_type *type = gsl_multiroot_fdfsolver_gnewton);
		~MultiDerivativeSolver();
		int Set(gsl_multiroot_function_fdf *function, const gsl_vector *x);
		int Iterate() override;
		int Solve() override;
		gsl_vector *Root() override;
		gsl_vector *Value() override;
		gsl_vector *StepSize() override;
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
		 * @brief `gsl_vector_alloc_from_block`, creates a vector as a slice of an existing `block`.
		 *
		 * @param block target block
		 * @param offset offset to the data block
		 * @param n length of the vector
		 * @param stride step-size of the vector
		 * @return pointer
		 */
		gsl_vector *VectorAllocFromBlock(gsl_block *block, const size_t offset, const size_t n, const size_t stride = 1);

		/**
		 * @brief `gsl_vector_alloc_row_from_matrix`, allocates a new `gsl_vector` which points to the `i`-th row of the `matrix`.
		 *
		 * @param matrix target matrix
		 * @param i row index
		 * @return pointer
		 */
		gsl_vector *VectorAllocRowFromMatrix(gsl_matrix *matrix, const size_t i);

		/**
		 * @brief `gsl_vector_alloc_col_from_matrix`, allocates a new `gsl_vector` which points to the `j`-th column of the `matrix`.
		 *
		 * @param matrix target matrix
		 * @param j column index
		 * @return pointer
		 */
		gsl_vector *VectorAllocColFromMatrix(gsl_matrix *matrix, const size_t j);

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
	 * @brief Dot product of vector `x` and `y`. \f$\vec{x}\cdot\vec{y}\f$
	 *
	 * @param x vector
	 * @param y vector
	 * @param dimension dimension of the vector
	 * @return result
	 */
	double Dot(const double x[], const double y[], size_t dimension = 3);

	/**
	 * @brief Dot product of vector `x`. \f$\vec{x}\cdot\vec{x}\f$
	 *
	 * @param x vector
	 * @param dimension dimension of the vector
	 * @return result
	 */
	double Dot(const double x[], size_t dimension = 3);

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

	double TriangleArea(const double a, const double b, const double c);
	double TriangleArea(const double x[], const double y[], const double z[]);

	/**
	 * @brief Rotate vector `x` around the `axis` by `angle`.
	 *
	 * @param x 3 dimensional vector
	 * @param axis the rotation axis.
	 * @param angle in rad
	 * @return status
	 */
	int RotateAroundAxis(double x[], Axis axis, double angle);

	/**
	 * @brief
	 *
	 * @param cartesian 4 or 8 dimensional vector
	 * @param spherical 4 or 8 dimensional vector
	 * @param dimension 4 or 8.
	 * @return status
	 */
	int CartesianToSpherical(const double cartesian[], double spherical[], size_t dimension = 8);

	/**
	 * @brief
	 *
	 * @param x 4 or 8 dimensional vector
	 * @param dimension 4 or 8
	 * @return status
	 */
	int CartesianToSpherical(double x[], size_t dimension = 8);

	/**
	 * @brief
	 *
	 * @param spherical 4 or 8 dimensional vector
	 * @param cartesian 4 or 8 dimensional vector
	 * @param dimension 4 or 8.
	 * @return status
	 */
	int SphericalToCartesian(const double spherical[], double cartesian[], size_t dimension = 8);

	/**
	 * @brief
	 *
	 * @param x 4 or 8 dimensional vector
	 * @param dimension 4 or 8.
	 * @return status
	 */
	int SphericalToCartesian(double x[], size_t dimension = 8);

	/**
	 * @brief Return `1` if `x`, `y` have opposite signs, else `0`.
	 *
	 * @param x number
	 * @param y number
	 * @return result
	 */
	int OppositeSign(double x, double y);

	/**
	 * @brief The modified spherical coordiante maps the polar angle from \f$\theta\in[0,\pi]\f$ to \f[
	 * \theta'\equiv\begin{cases}\theta & (\theta\leq\pi/2)\\
	 * \theta-\pi & (\theta>\pi/2)\end{cases}.
	 * \f]
	 * @param theta_0
	 * @param y
	 */
	void MapTheta(const double theta_0, double *y);

	/**
	 * @brief Return `phi` in \f$[0, 2\pi)\f$.
	 *
	 * @param phi
	 */
	void ModBy2Pi(double &phi);

	/**
	 * @brief Similar to `ModBy2Pi`, but return `phi` in \f$(-\pi, \pi]\f$.
	 *
	 * @param phi
	 * @return result
	 */
	double PhiDifference(double phi);

	/**
	 * @brief \f[\int_y^x(a_5+b_5t)^{p_5/2}\prod_{i=1}^4(a_i+b_it)^{-1/2}dt\f].
	 *
	 * @return result
	 */
	double EllipticIntegral(int p5, double y, double x, double a5, double b5, double a1, double b1, double a2, double b2, double a3, double b3, double a4 = 1., double b4 = 0.);

	/**
	 * @brief \f[\int_y^x(a_5+b_5t)^{p_5/2}(f+gt+ht^2)^{-1/2}\prod_{i=1,4}(a_i+b_it)^{-1/2}dt\f].
	 *
	 * @return result
	 */
	double EllipticIntegral2Complex(int p5, double y, double x, double a5, double b5, double f, double g, double h, double a1, double b1, double a4 = 1., double b4 = 0.);

	/**
	 * @brief \f[\int_y^x(a_5+b_5t)^{p_5/2}\prod_{i=1}^2(f_i+g_it+h_it^2)^{-1/2}dt\f].
	 *
	 * @return result
	 */
	double EllipticIntegral4Complex(int p5, double y, double x, double a5, double b5, double f1, double g1, double h1, double f2, double g2, double h2);

	double Carlson_RC(double x, double y, gsl_mode_t mode = GSL_PREC_DOUBLE);
	double Carlson_RJ(double x, double y, double z, double p, gsl_mode_t mode = GSL_PREC_DOUBLE);
} // namespace SBody

#endif
