/**
 * @file Utility.h
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 * The MultiFunctionSolver::Hybrid functions are based on code from `gsl/multiroots/hybrid.c`.
 *
 * The MultiFunctionSolver::Dnewton functions are based on code from `gsl/multiroots/dnewton.c`.
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * ----------------------------------------------------------------------------------------------------
 *
 * The function PolySolveQuartic is based on `quartic_roots` from `boost/math/tools/quartic_roots.hpp`.
 *
 * Boost Software License - Version 1.0 - August 17th, 2003
 *
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef SBODY_UTILITY_H
#define SBODY_UTILITY_H

#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_multimin.h>
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

	/// \f$\pi/3\f$.
	constexpr double M_PI_3 = 1.047197551196597746154214461093167628;

	/// \f$3\pi/4\f$.
	constexpr double M_3PI_4 = 2.35619449019234492884698253745962716;

	/// \f$2\pi\f$.
	constexpr double M_2PI = 6.28318530717958647692528676655900576;

	/// \f$\pi^2\f$.
	constexpr double M_PI2 = 9.86960440108935861883449099987615111;

	/// \f$\sqrt{27}\f$
	constexpr double M_SQRT27 = 5.19615242270663188058233902451761710;

	/// Epsilon, a small value. \f$\varepsilon\f$.
	constexpr double EPSILON = 1e-10;

	/// \f$\sin\varepsilon\f$.
	constexpr double SIN_EPSILON = 1e-10;

	/// \f$\cos\varepsilon\f$.
	constexpr double COS_EPSILON = 0.999999999999999999995;

	/// Area of a circle with radius of GSL_SQRT_DBL_EPSILON. \f$\pi\varepsilon^2\f$
	constexpr double EPSILON_CIRCLE_AREA = M_PI * GSL_DBL_EPSILON;

	constexpr double GSL_ROOT3_2_DBL_EPSILON = 3.666852862501036033408990023698041e-11;

	/// Sample number.
	constexpr int SAMPLE_NUMBER = 100;

	/// The angle corresponding to the sample number.
	constexpr double ANGLE_INTERVAL = M_2PI / SAMPLE_NUMBER;

	/// Area of the regular polygon with SAMPLE_NUMBER edges.
	constexpr double EPSILON_POLYGON_AREA = 0.06279051952931337 / ANGLE_INTERVAL * EPSILON_CIRCLE_AREA;

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
		~FunctionSolver();
		int Set(gsl_function *function, double lower, double upper);
		int Iterate() override;
		int Solve();
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
		int Solve(double epsabs, int max_iteration = 256);
		double Root() override;
	};

	class MultiSolver {
	  public:
		virtual int Iterate() = 0;
		virtual int Solve(double epsabs, int max_iteration) = 0;
		virtual int Solve(double epsabs, double epsrel, int max_iteration) = 0;
		virtual gsl_vector *Root() = 0;
		virtual gsl_vector *Value() = 0;
		virtual gsl_vector *StepSize() = 0;
	};

	struct HybridState {
		gsl_vector *x_trial;
		gsl_vector *f_trial;
		gsl_matrix *jacobian;
		gsl_matrix *lu;
		gsl_permutation *permutation;
		gsl_vector *newton;
		gsl_vector *gradient;
		double trust_radius;
		double epsilon_coefficient;
		bool directional;
	};

	struct HybridAdditionState {
		gsl_vector *x_trial;
		gsl_vector *x_trial2;
		gsl_vector *f_trial;
		gsl_vector *f_trial2;
		gsl_matrix *jacobian;
		gsl_matrix *lu;
		gsl_permutation *permutation;
		gsl_vector *newton;
		gsl_vector *gradient;
		double trust_radius;
		double epsilon_coefficient;
		double gradient_coefficient;
		bool directional;
	};

	struct DNewtonState {
		gsl_vector *x_trial;
		gsl_vector *f_trial;
		gsl_matrix *jacobian;
		gsl_matrix *lu;
		gsl_permutation *permutation;
		double trust_radius;
		double epsilon_coefficient;
		bool directional;
	};
	struct DNewtonRotationTranslationState {
		gsl_vector *x_trial;
		gsl_vector *dx_translation;
		gsl_vector *f_trial;
		gsl_matrix *jacobian;
		gsl_matrix *lu;
		gsl_permutation *permutation;
		double theta_obs;
		double sin_theta_obs;
		double cos_theta_obs;
		double r_obj;
		double sin_theta_obj;
		double cos_theta_obj;
		double phi_obj;
		double projected_x;
		double projected_y;
		double iota_obj;
		double trust_radius;
		double epsilon_coefficient;
		bool directional;
		bool trace_to_plane;
	};
	struct D2NewtonState {
		gsl_vector *x_trial;
		gsl_vector *f_trial;
		gsl_vector *jacobian;
		gsl_matrix *hessian;
		gsl_matrix *lu;
		gsl_permutation *permutation;
		double trust_radius;
	};
	struct ConjugateGradientState {
		double iteration_coefficient;
		gsl_vector *x_trial;
		gsl_vector *f_trial;
		gsl_vector *last_dx;
		double gradient_norm;
		gsl_matrix *jacobian;
		gsl_matrix *lu;
		gsl_permutation *permutation;
		double trust_radius;
	};
	struct TriangleState {
		gsl_vector *x_trial;
		gsl_vector *f_trial;
		gsl_vector *x_a;
		gsl_vector *x_b;
		gsl_vector *f_a;
		gsl_vector *f_b;
	};
	struct DirectionState {
		int directional_num;
		double trust_radius;
		double delta_angle;
		double central_angle;
		gsl_vector *x_trial, *x_rec;
		gsl_vector *f_trial, *f_rec;
	};

	extern const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_hybrid;
	extern const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_hybrid_addition;
	extern const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_dnewton;
	extern const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_dnewton_rotation;
	extern const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_dnewton_translation;
	extern const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_d2newton;
	extern const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_gradient;
	extern const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_conjugate;
	extern const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_triangle;
	extern const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_direction;

	class MultiFunctionSolver : public MultiSolver {
	  private:
		gsl_multiroot_fsolver *solver_;
		static int ScaleX(gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f);
		static int OneSidedJacobian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *f, double epsrel, gsl_matrix *jacobian);
		static int TwoSidedJacobian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *f, double epsrel, gsl_matrix *jacobian);
		static int OneSidedDirectionalJacobian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *direction, const gsl_vector *f, double epsrel, gsl_matrix *jacobian);
		static int TwoSidedDirectionalJacobian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *direction, const gsl_vector *f, double epsrel, gsl_matrix *jacobian);
		static int Hessian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *f, double epsrel, gsl_vector *jacobian, gsl_matrix *hessian);
		static void ComputeDiag(const gsl_matrix *J, gsl_vector *diag);
		static void UpdateDiag(const gsl_matrix *J, gsl_vector *diag);
		static double ScaledEnorm(const gsl_vector *d, const gsl_vector *f);
		static int Dogleg(const gsl_matrix *r, const gsl_vector *qtf, const gsl_vector *diag, double delta, gsl_vector *newton, gsl_vector *gradient, gsl_vector *p);
		static int Dogleg(double trust_radius, double newton_norm, double newton_norm2, double gradient_norm, double gradient_norm2, double newton_gradient_dot, const gsl_vector *newton, const gsl_vector *gradient, gsl_vector *dx);

	  public:
		MultiFunctionSolver(size_t n, const gsl_multiroot_fsolver_type *type = gsl_multiroot_fsolver_dnewton);
		~MultiFunctionSolver();
		int Set(gsl_multiroot_function *function, const gsl_vector *x);
		int Set(gsl_multiroot_function *function, const gsl_vector *x, double theta_obs, double sin_theta_obs, double cos_theta_obs, double r_obj, double sin_theta_obj, double cos_theta_obj, double phi_obj, double sin_phi_obj, double cos_phi_obj, bool trace_to_plane);
		int Iterate() override;
		int Solve(double epsabs, int max_iteration = 256) override;
		int Solve(double epsabs, double epsrel, int max_iteration = 256) override;
		gsl_vector *Root() override;
		gsl_vector *Value() override;
		gsl_vector *StepSize() override;
		static int HybridAlloc(void *vstate, size_t n);
		static int HybridSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int HybridIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static void HybridFree(void *vstate);
		static int HybridAdditionAlloc(void *vstate, size_t n);
		static int HybridAdditionSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int HybridAdditionIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static void HybridAdditionFree(void *vstate);
		static int DNewtonAlloc(void *vstate, size_t n);
		static int DNewtonSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int DNewtonIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static void DNewtonFree(void *vstate);
		static int DNewtonRotationTranslationAlloc(void *vstate, size_t n);
		static int DNewtonRotationTranslationSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int DNewtonRotationIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int DNewtonTranslationIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static void DNewtonRotationTranslationFree(void *vstate);
		static int D2NewtonAlloc(void *vstate, size_t n);
		static int D2NewtonSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int D2NewtonIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static void D2NewtonFree(void *vstate);
		static int GradientIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int ConjugateGradientAlloc(void *vstate, size_t n);
		static int ConjugateGradientSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int ConjugateGradientIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static void ConjugateGradientFree(void *vstate);
		static int TriangleAlloc(void *vstate, size_t n);
		static int TriangleSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int TriangleRotationIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int TriangleLongestIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static void TriangleFree(void *vstate);
		static int DirectionAlloc(void *vstate, size_t n);
		static int DirectionSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static int DirectionIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx);
		static void DirectionFree(void *vstate);
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
		int Solve(double epsabs, double epsrel, int max_iteration = 128) override;
		gsl_vector *Root() override;
		gsl_vector *Value() override;
		gsl_vector *StepSize() override;
	};

	class MultiMinimizer {
	  public:
		virtual int Iterate() = 0;
		virtual int Solve(double epsabs) = 0;
		virtual gsl_vector *Root() = 0;
		virtual double Value() = 0;
	};

	class MultiFunctionMinimizer : public MultiMinimizer {
	  private:
		gsl_multimin_fminimizer *solver_;

	  public:
		MultiFunctionMinimizer(size_t n, const gsl_multimin_fminimizer_type *type = gsl_multimin_fminimizer_nmsimplex2rand);
		~MultiFunctionMinimizer();
		int Set(gsl_multimin_function *function, const gsl_vector *x, const gsl_vector *step_size);
		int Iterate() override;
		int Solve(double epsabs) override;
		gsl_vector *Root() override;
		double Value() override;
		double StepSize();
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

	int CoordinateOrthogonalization(const gsl_vector *x, gsl_matrix *coordinate);

	/**
	 * @brief Square root of `x`. Return `0` if `x` is negative.
	 *
	 * @param x number
	 * @return result
	 */
	double SquareRoot(double x);

	/**
	 * @brief Square of `x` with the sign of `x`.
	 *
	 * @param x number
	 * @return result
	 */
	double SignSquare(double x);

	/**
	 * @brief Square root of `abs(x)` with the sign of `x`.
	 *
	 * @param x number
	 * @return result
	 */
	double SignSquareRoot(double x);

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

	/**
	 * @brief Calculate the area of a triangle
	 *
	 * @param a side length
	 * @param b side length
	 * @param c side length
	 * @return result
	 */
	double TriangleArea(double a, double b, double c);

	/**
	 * @brief Calculate the area of a triangle
	 *
	 * @param a side length
	 * @param b side length
	 * @param c side length
	 * @return result
	 */
	double TriangleArea(const double x[], const double y[], const double z[]);

	/**
	 * @brief Calculate whether point \f$(x, y)\f$ locates inside the triangle \f$\Delta abc\f$.
	 *
	 * @param xa x coordiante of the vertex a
	 * @param ya y coordiante of the vertex a
	 * @param xb x coordiante of the vertex b
	 * @param yb y coordiante of the vertex b
	 * @param xc x coordiante of the vertex c
	 * @param yc y coordiante of the vertex c
	 * @param x x coordiante of the point
	 * @param y y coordiante of the point
	 * @return result
	 */
	bool PointInTriangle(double xa, double ya, double xb, double yb, double xc, double yc, double x = 0.0, double y = 0.0);

	/**
	 * @brief Calculate whether point \f$p\f$ locates inside the triangle \f$\Delta abc\f$.
	 *
	 * @param a vertex a
	 * @param b vertex b
	 * @param c vertex c
	 * @param p point p
	 * @return result
	 */
	bool PointInTriangle(const gsl_vector *a, const gsl_vector *b, const gsl_vector *c, const gsl_vector *p = nullptr);

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

	double SphericalAngle(double cos_theta_x, double cos_theta_y, double delta_theta_xy, double delta_phi_xy);

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
	 * @param theta_0 theta in the last integration step
	 * @param y 8 dimensional vector
	 * @return status
	 */
	int MapTheta(const double theta_0, double y[]);

	/**
	 * @brief Return `phi` in \f$[0, 2\pi)\f$.
	 *
	 * @param phi \f$\phi\f$
	 * @return result
	 */
	double ModBy2Pi(double phi);

	/**
	 * @brief Similar to `ModBy2Pi`, but return `phi` in \f$[-\pi, \pi)\f$.
	 *
	 * @param phi \f$\phi\f$
	 * @return result
	 */
	double PhiDifference(double phi);

	int AMinusPlusB(double a, double b, double &a_minus_b, double &a_plus_b);

	/**
	 * @brief Linear interpolation of points (`x0`, `y0`) and (`x1`, `y1`) at `x`. \f[y=\frac{y_0(x_1-x)+y_1(x-x_0)}{x_1-x_0}\f]
	 *
	 * @param x position to evaluate \f$y\f$
	 * @param x0 \f$x_0\f$
	 * @param x1 \f$x_1\f$
	 * @param y0 \f$y_0\f$
	 * @param y1 \f$y_1\f$
	 * @return result
	 */
	double LinearInterpolation(double x, double x0, double x1, double y0, double y1);

	/**
	 * @brief Linear interpolation of vectors (`x0`, `y0`) and (`x1`, `y1`) at `x`, stored in `y`. \f[y=\frac{y_0(x_1-x)+y_1(x-x_0)}{x_1-x_0}\f]
	 *
	 * @param x position to evaluate \f$y\f$
	 * @param x0 \f$x_0\f$
	 * @param x1 \f$x_1\f$
	 * @param y0 \f$y_0\f$
	 * @param y1 \f$y_1\f$
	 * @param y \f$y\f$
	 * @param size size of the vectors
	 * @return status
	 */
	int LinearInterpolation(double x, double x0, double x1, const double y0[], const double y1[], double y[], size_t size);

	/**
	 * @brief Linear interpolation of two spherical positions `y0` and `y1` at `t`, stored in `y`. \f[y=\frac{\text{SphericalToCartesian}(y_0)(t_1-t)+\text{SphericalToCartesian}(y_1)(t-t_0)}{t_1-t_0}\f]
	 *
	 * @param t time to evaluate \f$y\f$
	 * @param t0 \f$t_0\f$
	 * @param t1 \f$t_1\f$
	 * @param y0 8 dimensional vector, spherical, \f$y_0\f$
	 * @param y1 8 dimensional vector, spherical, \f$y_1\f$
	 * @param y 8 dimensional vector, cartesian, \f$y\f$
	 * @return status
	 */
	int InterpolateSphericalPositionToCartesian(double t, double t0, double t1, const double y0[], const double y1[], double y[]);

	/**
	 * @brief Calculate the bolometric flux. The result is \f[F_\text{obs}=\frac{L_\text{src}}{(1+z)^2}\frac{\Omega_\text{src}}{\Omega_\text{obs}}.\f] The additional \f$(1+z)^{-1}\f$ comes from the expansion of the frequency that \f$(1+z)^{-1}=\frac{\nu_\text{obs}}{\nu_\text{src}}\f$.
	 *
	 * @param luminosity intrinsic luminosity of the source, \f$L_\text{src}\f$
	 * @param magnification magnification of the source, including the relativistic redshift and beaming, and gravitational lensing. \f$\frac{\Omega_\text{src}}{(1+z)\Omega_\text{obs}}\f$
	 * @param redshift relativistic redshift of the source, \f$1+z\f$
	 * @return result
	 */
	double Flux(double luminosity, double magnification, double redshift);

	/**
	 * @brief Calculate the flux density. The result is \f[F_{\nu,\text{obs}}=\frac{L_{\nu,\text{src}}}{1+z}\frac{\Omega_\text{src}}{\Omega_\text{obs}}.\f]
	 *
	 * @param spectral_density intrinsic spectral density of the source, \f$L_{\nu,\text{src}}\f$
	 * @param magnification magnification of the source, including the relativistic redshift and beaming, and gravitational lensing. \f$\frac{\Omega_\text{src}}{(1+z)\Omega_\text{obs}}\f$
	 * @return result
	 */
	double FluxDensity(double spectral_density, double magnification);

	int PolishQuadraticRoot(double a, double b, double roots[], int root_num);
	int PolishCubicRoot(double a, double b, double c, double roots[], int root_num);
	int PolishQuarticRoot(double a, double b, double c, double d, double roots[], int root_num);

	int PolySolveQuarticWithZero(double a, double b, double c, double offset, double roots[]);

	/**
	 * @brief Solve for real roots of the quartic equation \f$x^4+ax^3+bx^2+cx=0\f$. The roots are returned in `x0`, `x1`, `x2`, and `x3` and satisfied \f$x_0\leq x_1\leq x_2\leq x_3\f$. The function is based on `quartic_roots` from `boost/math/tools/quartic_roots.hpp`.
	 *
	 * @param a \f$a\f$
	 * @param b \f$b\f$
	 * @param c \f$c\f$
	 * @param d \f$d\f$
	 * @param roots
	 *
	 * @return number of real roots
	 */
	int PolySolveQuartic(double a, double b, double c, double d, double roots[]);

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
