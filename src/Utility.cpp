/**
 * @file Utility.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Utility.hpp"

#include <cmath>
#include <vector>

#include <boost/algorithm/algorithm.hpp>
#include <fmt/core.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_vector.h>

#include "IO.hpp"

using namespace std;

namespace SBody {
	// DerivativeSolver
	DerivativeSolver::DerivativeSolver(const gsl_root_fdfsolver_type *type) {
		solver_ = gsl_root_fdfsolver_alloc(type);
	}
	DerivativeSolver::DerivativeSolver(gsl_function_fdf *function, long double root, const gsl_root_fdfsolver_type *type) {
		solver_ = gsl_root_fdfsolver_alloc(type);
		gsl_root_fdfsolver_set(solver_, function, root);
	}
	DerivativeSolver::~DerivativeSolver() {
		gsl_root_fdfsolver_free(solver_);
	}
	int DerivativeSolver::Set(gsl_function_fdf *function, long double root) {
		return gsl_root_fdfsolver_set(solver_, function, root);
	}
	int DerivativeSolver::Iterate() {
		return gsl_root_fdfsolver_iterate(solver_);
	}
	int DerivativeSolver::Solve(long double epsabs, int max_iteration) {
		for (; max_iteration > 0 && gsl_root_test_residual(GSL_FN_FDF_EVAL_F(solver_->fdf, solver_->root), absolute_accuracy) != GSL_SUCCESS; --max_iteration)
			if (int status = gsl_root_fdfsolver_iterate(solver_); status != GSL_SUCCESS)
				return status;
		if (max_iteration <= 0)
			return GSL_EMAXITER;
		return GSL_SUCCESS;
	}
	long double DerivativeSolver::Root() {
		return gsl_root_fdfsolver_root(solver_);
	}

	// MultiFunctionSolver
	MultiFunctionSolver::MultiFunctionSolver(size_t n, const gsl_multiroot_fsolver_type *type) {
		solver_ = gsl_multiroot_fsolver_alloc(type, n);
	}
	MultiFunctionSolver::~MultiFunctionSolver() {
		gsl_multiroot_fsolver_free(solver_);
	}
	int MultiFunctionSolver::Set(gsl_multiroot_function *function, const gsl_vector *x) {
		if (solver_->type == gsl_multiroot_fsolver_sbody_dnewton_rotation || solver_->type == gsl_multiroot_fsolver_sbody_dnewton_translation)
			return GSL_EINVAL;
		return gsl_multiroot_fsolver_set(solver_, function, x);
	}
	int MultiFunctionSolver::Set(gsl_multiroot_function *function, const gsl_vector *x, long double theta_obs, long double sin_theta_obs, long double cos_theta_obs, long double r_obj, long double sin_theta_obj, long double cos_theta_obj, long double phi_obj, long double sin_phi_obj, long double cos_phi_obj, bool trace_to_plane) {
		if (int status = gsl_multiroot_fsolver_set(solver_, function, x); status != GSL_SUCCESS)
			return status;
		if (solver_->type == gsl_multiroot_fsolver_sbody_dnewton_rotation || solver_->type == gsl_multiroot_fsolver_sbody_dnewton_translation) {
			auto state = static_cast<DNewtonRotationTranslationState *>(solver_->state);
			state->theta_obs = theta_obs;
			state->sin_theta_obs = sin_theta_obs;
			state->cos_theta_obs = cos_theta_obs;
			state->r_obj = r_obj;
			state->sin_theta_obj = sin_theta_obj;
			state->cos_theta_obj = cos_theta_obj;
			state->phi_obj = phi_obj;
			state->projected_x = sin_theta_obj * sin_phi_obj;
			state->projected_y = cos_theta_obj * sin_theta_obs - sin_theta_obj * cos_phi_obj * cos_theta_obs;
			state->iota_obj = atan2(state->projected_y, state->projected_x);
			state->trace_to_plane = trace_to_plane;
			return GSL_SUCCESS;
		}
		return GSL_EINVAL;
	}
	int MultiFunctionSolver::Iterate() {
		return gsl_multiroot_fsolver_iterate(solver_);
	}
	int MultiFunctionSolver::Solve(long double epsabs, int max_iteration) {
		for (; max_iteration > 0 && gsl_multiroot_test_residual(gsl_multiroot_fsolver_f(solver_), epsabs) != GSL_SUCCESS; --max_iteration) {
			if (int status = gsl_multiroot_fsolver_iterate(solver_); status != GSL_SUCCESS)
				return status;
			if (any_of(solver_->x->data, solver_->x->data + solver_->x->size, [](long double x) { return !isfinite(x); }))
				return GSL_EFAILED;
		}
		if (max_iteration <= 0)
			return GSL_EMAXITER;
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::Solve(long double epsabs, long double epsrel, int max_iteration) {
		for (; max_iteration > 0 && (gsl_multiroot_test_delta(gsl_multiroot_fsolver_dx(solver_), gsl_multiroot_fsolver_root(solver_), epsabs, epsrel) != GSL_SUCCESS || gsl_multiroot_test_residual(gsl_multiroot_fsolver_f(solver_), epsabs) != GSL_SUCCESS); --max_iteration) {
			if (int status = gsl_multiroot_fsolver_iterate(solver_); status != GSL_SUCCESS)
				return status;
			if (any_of(solver_->x->data, solver_->x->data + solver_->x->size, [](long double x) { return !isfinite(x); }))
				return GSL_EFAILED;
		}
		if (max_iteration <= 0)
			return GSL_EMAXITER;
		return GSL_SUCCESS;
	}
	gsl_vector *MultiFunctionSolver::Root() {
		return gsl_multiroot_fsolver_root(solver_);
	}
	gsl_vector *MultiFunctionSolver::Value() {
		return gsl_multiroot_fsolver_f(solver_);
	}
	gsl_vector *MultiFunctionSolver::StepSize() {
		return gsl_multiroot_fsolver_dx(solver_);
	}
	int MultiFunctionSolver::ScaleX(gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f) {
		for (int status = GSL_MULTIROOT_FN_EVAL(function, x, f); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x, f)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_scale(x, gsl_vector_get(f, 0));
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::OneSidedJacobian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *f, long double epsrel, gsl_matrix *jacobian) {
		const size_t n = x->size;
		const size_t m = f->size;
		const size_t n1 = jacobian->size1;
		const size_t n2 = jacobian->size2;
		if (m != n1 || n != n2)
			GSL_ERROR("function and jacobian are not conformant", GSL_EBADLEN);
		gsl_vector *x1 = gsl_vector_alloc(n);
		if (x1 == nullptr)
			GSL_ERROR("failed to allocate space for x1 workspace", GSL_ENOMEM);
		gsl_vector_memcpy(x1, x); /* copy x into x1 */
		for (size_t j = 0; j < n; ++j) {
			long double x_j = gsl_vector_get(x, j);
			long double dx = abs(x_j) >= 1. ? -epsrel * x_j : -epsrel * GSL_SIGN(x_j);
			gsl_vector_set(x1, j, x_j + dx);
			gsl_vector_view col_j = gsl_matrix_column(jacobian, j);
			if (int status = GSL_MULTIROOT_FN_EVAL(F, x1, &col_j.vector); status != GSL_SUCCESS) {
				gsl_vector_free(x1);
				return GSL_EBADFUNC;
			}
			gsl_blas_daxpy(-1., f, &col_j.vector);
			gsl_blas_dscal(1. / dx, &col_j.vector);
			gsl_vector_set(x1, j, x_j);
		}
		gsl_vector_free(x1);
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::TwoSidedJacobian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *f, long double epsrel, gsl_matrix *jacobian) {
		const size_t n = x->size;
		const size_t m = f->size;
		const size_t n1 = jacobian->size1;
		const size_t n2 = jacobian->size2;
		int status = GSL_SUCCESS;
		if (m != n1 || n != n2)
			GSL_ERROR("function and jacobian are not conformant", GSL_EBADLEN);
		gsl_vector *x1 = gsl_vector_alloc(n);
		if (x1 == nullptr)
			GSL_ERROR("failed to allocate space for x1 workspace", GSL_ENOMEM);
		gsl_vector *f_trial = gsl_vector_alloc(n);
		if (f_trial == nullptr) {
			gsl_vector_free(x1);
			GSL_ERROR("failed to allocate space for f_trial workspace", GSL_ENOMEM);
		}
		gsl_vector_memcpy(x1, x); /* copy x into x1 */
		for (size_t j = 0; j < n; ++j) {
			long double x_j = gsl_vector_get(x, j);
			gsl_vector_view col_j = gsl_matrix_column(jacobian, j);
			long double dx = abs(x_j) >= 1. ? -epsrel * x_j : -epsrel * GSL_SIGN(x_j);
			for (status = GSL_CONTINUE; status != GSL_SUCCESS; dx *= 0.5) {
				if (abs(dx) < abs(x_j) * GSL_DBL_EPSILON) {
					status = GSL_EBADFUNC;
					break;
				}
				gsl_vector_set(x1, j, x_j + dx);
				if (status = GSL_MULTIROOT_FN_EVAL(F, x1, &col_j.vector); status != GSL_SUCCESS)
					if (status != GSL_EDOM) {
						status = GSL_EBADFUNC;
						break;
					}
			}
			gsl_vector_set(x1, j, x_j - dx);
			if (status = GSL_MULTIROOT_FN_EVAL(F, x1, f_trial); status != GSL_SUCCESS) {
				status = GSL_SUCCESS;
				gsl_blas_daxpy(-1., f, &col_j.vector);
				gsl_blas_dscal(1. / dx, &col_j.vector);
			} else {
				gsl_blas_daxpy(-1., f_trial, &col_j.vector);
				gsl_blas_dscal(0.5 / dx, &col_j.vector);
			}
			gsl_vector_set(x1, j, x_j);
		}
		gsl_vector_free(x1);
		gsl_vector_free(f_trial);
		return status;
	}
	int MultiFunctionSolver::OneSidedDirectionalJacobian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *direction, const gsl_vector *f, long double epsrel, gsl_matrix *jacobian) {
		const size_t n = x->size;
		const size_t m = f->size;
		const size_t n1 = jacobian->size1;
		const size_t n2 = jacobian->size2;
		int status = GSL_SUCCESS;
		if (m != n1 || n != n2)
			GSL_ERROR("function and jacobian are not conformant", GSL_EBADLEN);
		gsl_vector *x1 = gsl_vector_alloc(n);
		if (x1 == nullptr)
			GSL_ERROR("failed to allocate space for x1 workspace", GSL_ENOMEM);
		gsl_vector *f_trial = gsl_vector_alloc(n);
		if (f_trial == nullptr) {
			gsl_vector_free(x1);
			GSL_ERROR("failed to allocate space for f_trial workspace", GSL_ENOMEM);
		}
		gsl_matrix *coordinate = gsl_matrix_alloc(n, n);
		if (coordinate == nullptr) {
			gsl_vector_free(x1);
			gsl_vector_free(f_trial);
			GSL_ERROR("failed to allocate space for coordinate workspace", GSL_ENOMEM);
		}
		if (status = CoordinateOrthogonalization(direction, coordinate); status != GSL_SUCCESS)
			return status;
		gsl_matrix_set_zero(jacobian);
		long double x_norm = gsl_blas_dnrm2(x), dx_limit = x_norm * GSL_DBL_EPSILON;
		for (size_t j = 0; j < n; ++j) {
			gsl_vector_view coordinate_j = gsl_matrix_column(coordinate, j);
			double dot_product;
			gsl_blas_ddot(x, &coordinate_j.vector, &dot_product);
			long double dx = -epsrel * max(1.l, x_norm) * GSL_SIGN(dot_product);
			gsl_blas_daxpy(dx, &coordinate_j.vector, x1);
			for (status = GSL_CONTINUE; status != GSL_SUCCESS; dx *= 0.5) {
				if (abs(dx) < dx_limit) {
					status = GSL_EBADFUNC;
					break;
				}
				gsl_vector_memcpy(x1, x);
				gsl_blas_daxpy(dx, &coordinate_j.vector, x1);
				if (status = GSL_MULTIROOT_FN_EVAL(F, x1, f_trial); status != GSL_SUCCESS)
					if (status != GSL_EDOM) {
						status = GSL_EBADFUNC;
						break;
					}
			}
			gsl_blas_daxpy(-1., f, f_trial);
			gsl_blas_dscal(1. / dx, f_trial);
			gsl_blas_dger(1., f_trial, &coordinate_j.vector, jacobian);
		}
		gsl_vector_free(x1);
		gsl_vector_free(f_trial);
		gsl_matrix_free(coordinate);
		return status;
	}
	int MultiFunctionSolver::TwoSidedDirectionalJacobian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *direction, const gsl_vector *f, long double epsrel, gsl_matrix *jacobian) {
		const size_t n = x->size;
		const size_t m = f->size;
		const size_t n1 = jacobian->size1;
		const size_t n2 = jacobian->size2;
		int status = GSL_SUCCESS;
		if (m != n1 || n != n2)
			GSL_ERROR("function and jacobian are not conformant", GSL_EBADLEN);
		gsl_vector *x1 = gsl_vector_alloc(n);
		if (x1 == nullptr)
			GSL_ERROR("failed to allocate space for x1 workspace", GSL_ENOMEM);
		gsl_vector *f_trial_plus = gsl_vector_alloc(n);
		if (f_trial_plus == nullptr) {
			gsl_vector_free(x1);
			GSL_ERROR("failed to allocate space for f_trial_plus workspace", GSL_ENOMEM);
		}
		gsl_vector *f_trial_minus = gsl_vector_alloc(n);
		if (f_trial_minus == nullptr) {
			gsl_vector_free(x1);
			gsl_vector_free(f_trial_plus);
			GSL_ERROR("failed to allocate space for f_trial_minus workspace", GSL_ENOMEM);
		}
		gsl_matrix *coordinate = gsl_matrix_alloc(n, n);
		if (coordinate == nullptr) {
			gsl_vector_free(x1);
			gsl_vector_free(f_trial_plus);
			gsl_vector_free(f_trial_minus);
			GSL_ERROR("failed to allocate space for coordinate workspace", GSL_ENOMEM);
		}
		if (status = CoordinateOrthogonalization(direction, coordinate); status != GSL_SUCCESS)
			return status;
		gsl_matrix_set_zero(jacobian);
		long double x_norm = gsl_blas_dnrm2(x), dx_limit = x_norm * GSL_DBL_EPSILON;
		for (size_t j = 0; j < n; ++j) {
			gsl_vector_view coordinate_j = gsl_matrix_column(coordinate, j);
			double dot_product;
			gsl_blas_ddot(x, &coordinate_j.vector, &dot_product);
			long double dx = -epsrel * max(1.l, x_norm) * GSL_SIGN(dot_product);
			for (status = GSL_CONTINUE; status != GSL_SUCCESS; dx *= 0.5) {
				if (abs(dx) < dx_limit) {
					status = GSL_EBADFUNC;
					break;
				}
				gsl_vector_memcpy(x1, x);
				gsl_blas_daxpy(dx, &coordinate_j.vector, x1);
				if (status = GSL_MULTIROOT_FN_EVAL(F, x1, f_trial_plus); status != GSL_SUCCESS)
					if (status != GSL_EDOM) {
						status = GSL_EBADFUNC;
						break;
					}
			}
			gsl_blas_daxpy(-1., f, f_trial_plus);
			gsl_blas_dscal(1. / dx, f_trial_plus);
			gsl_vector_memcpy(x1, x);
			gsl_blas_daxpy(-dx, &coordinate_j.vector, x1);
			if (status = GSL_MULTIROOT_FN_EVAL(F, x1, f_trial_minus); status != GSL_SUCCESS) {
				status = GSL_SUCCESS;
				gsl_blas_dger(1., f_trial_plus, &coordinate_j.vector, jacobian);
			} else {
				gsl_blas_daxpy(-1., f, f_trial_minus);
				gsl_blas_daxpy(-1. / dx, f_trial_minus, f_trial_plus);
				gsl_blas_dger(0.5, f_trial_plus, &coordinate_j.vector, jacobian);
			}
		}
		gsl_vector_free(x1);
		gsl_vector_free(f_trial_plus);
		gsl_vector_free(f_trial_minus);
		gsl_matrix_free(coordinate);
		return status;
	}
	int MultiFunctionSolver::Hessian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *f, long double epsrel, gsl_vector *jacobian, gsl_matrix *hessian) {
		const size_t n = x->size;
		const size_t m = f->size;
		const size_t n1 = hessian->size1;
		const size_t n2 = hessian->size2;
		if (m != n1 || n != n2)
			GSL_ERROR("function and hessian are not conformant", GSL_EBADLEN);
		gsl_vector *x1 = gsl_vector_alloc(n);
		if (x1 == nullptr)
			GSL_ERROR("failed to allocate space for x1 workspace", GSL_ENOMEM);
		gsl_vector *f_trial = gsl_vector_alloc(n);
		if (f_trial == nullptr) {
			gsl_vector_free(x1);
			GSL_ERROR("failed to allocate space for f_trial workspace", GSL_ENOMEM);
		}
		const long double f_norm = gsl_blas_dnrm2(f);
		long double dx[n];
		// f_{0, 0} = f0
		// f_{-1,0} = f0 - g0 * dx0 + 0.5 * G00 * dx0^2
		// f_{1, 0} = f0 + g0 * dx0 + 0.5 * G00 * dx0^2
		// f_{0,-1} = f0 - g1 * dx1 + 0.5 * G11 * dx1^2
		// f_{0, 1} = f0 + g1 * dx1 + 0.5 * G11 * dx1^2
		// f_{-1,-1}= f0 - g0 * dx0 - g1 * dx1 + 0.5 * G00 * dx0^2 + 0.5 * G11 * dx1^2 + G01 * dx0 * dx1
		// f_{1, 1} = f0 + g0 * dx0 + g1 * dx1 + 0.5 * G00 * dx0^2 + 0.5 * G11 * dx1^2 + G01 * dx0 * dx1
		// f_{-1,1} = f0 - g0 * dx0 + g1 * dx1 + 0.5 * G00 * dx0^2 + 0.5 * G11 * dx1^2 - G01 * dx0 * dx1
		// f_{1,-1} = f0 + g0 * dx0 - g1 * dx1 + 0.5 * G00 * dx0^2 + 0.5 * G11 * dx1^2 - G01 * dx0 * dx1

		gsl_vector_memcpy(x1, x); // copy x into x1
		for (size_t i = 0; i < n; ++i) {
			long double x_i = gsl_vector_get(x, i);
			dx[i] = abs(x_i) >= 1. ? epsrel * x_i : epsrel * x_i;
			gsl_vector_set(x1, i, x_i + dx[i]);
			if (int f_stat = GSL_MULTIROOT_FN_EVAL(F, x1, f_trial); f_stat != GSL_SUCCESS) {
				gsl_vector_free(x1);
				gsl_vector_free(f_trial);
				return GSL_EBADFUNC;
			}
			long double f_trial_norm = gsl_blas_dnrm2(f_trial);
			gsl_vector_set(x1, i, x_i - dx[i]);
			if (int f_stat = GSL_MULTIROOT_FN_EVAL(F, x1, f_trial); f_stat != GSL_SUCCESS) {
				gsl_vector_free(x1);
				gsl_vector_free(f_trial);
				return GSL_EBADFUNC;
			}
			gsl_vector_set(x1, i, x_i);
			gsl_vector_set(jacobian, i, (f_trial_norm - gsl_blas_dnrm2(f_trial)) / (2. * dx[i]));
			gsl_matrix_set(hessian, i, i, (f_trial_norm + gsl_blas_dnrm2(f_trial) - 2. * f_norm) / Power2(dx[i]));
		}
		for (size_t i = 1; i < n; ++i) {
			gsl_vector_memcpy(x1, x);
			long double x_i = gsl_vector_get(x, i);
			for (size_t j = 0; j < i; ++j) {
				long double df = 0., x_j = gsl_vector_get(x, j);
				gsl_vector_set(x1, i, x_i + dx[i]);
				gsl_vector_set(x1, j, x_j + dx[j]);
				if (int f_stat = GSL_MULTIROOT_FN_EVAL(F, x1, f_trial); f_stat != GSL_SUCCESS) {
					gsl_vector_free(x1);
					gsl_vector_free(f_trial);
					return GSL_EBADFUNC;
				}
				df += gsl_blas_dnrm2(f_trial);
				gsl_vector_set(x1, j, x_j - dx[j]);
				if (int f_stat = GSL_MULTIROOT_FN_EVAL(F, x1, f_trial); f_stat != GSL_SUCCESS) {
					gsl_vector_free(x1);
					gsl_vector_free(f_trial);
					return GSL_EBADFUNC;
				}
				df -= gsl_blas_dnrm2(f_trial);
				gsl_vector_set(x1, i, x_i - dx[i]);
				if (int f_stat = GSL_MULTIROOT_FN_EVAL(F, x1, f_trial); f_stat != GSL_SUCCESS) {
					gsl_vector_free(x1);
					gsl_vector_free(f_trial);
					return GSL_EBADFUNC;
				}
				df += gsl_blas_dnrm2(f_trial);
				gsl_vector_set(x1, j, x_j + dx[j]);
				if (int f_stat = GSL_MULTIROOT_FN_EVAL(F, x1, f_trial); f_stat != GSL_SUCCESS) {
					gsl_vector_free(x1);
					gsl_vector_free(f_trial);
					return GSL_EBADFUNC;
				}
				df -= gsl_blas_dnrm2(f_trial);
				df /= 4. * dx[i] * dx[j];
				gsl_matrix_set(hessian, i, j, df);
				gsl_matrix_set(hessian, j, i, df);
			}
		}
		gsl_vector_free(x1);
		gsl_vector_free(f_trial);
		return GSL_SUCCESS;
	}
	void MultiFunctionSolver::ComputeDiag(const gsl_matrix *J, gsl_vector *diag) {
		size_t j, n = diag->size;
		for (j = 0; j < n; j++) {
			auto column = gsl_matrix_const_column(J, j);
			long double norm = gsl_blas_dnrm2(&column.vector);
			if (norm == 0.)
				norm = 1.0;
			gsl_vector_set(diag, j, norm);
		}
	}
	void MultiFunctionSolver::UpdateDiag(const gsl_matrix *J, gsl_vector *diag) {
		size_t j, n = diag->size;
		for (j = 0; j < n; j++) {
			auto column = gsl_matrix_const_column(J, j);
			long double norm = gsl_blas_dnrm2(&column.vector);
			if (norm == 0.)
				norm = 1.0;
			if (norm > gsl_vector_get(diag, j))
				gsl_vector_set(diag, j, norm);
		}
	}
	long double MultiFunctionSolver::ScaledEnorm(const gsl_vector *d, const gsl_vector *f) {
		long double e2 = 0.;
		size_t i, n = f->size;
		for (i = 0; i < n; i++)
			e2 += gsl_pow_2(gsl_vector_get(f, i) * gsl_vector_get(d, i));
		return sqrt(e2);
	}
	int MultiFunctionSolver::Dogleg(const gsl_matrix *r, const gsl_vector *qtf, const gsl_vector *diag, long double delta, gsl_vector *newton, gsl_vector *gradient, gsl_vector *p) {
		long double qnorm, gnorm, sgnorm, bnorm, temp;
		gsl_linalg_R_solve(r, qtf, newton);
		gsl_vector_scale(newton, -1.);
		if (qnorm = ScaledEnorm(diag, newton); qnorm <= delta) {
			gsl_vector_memcpy(p, newton);
			return GSL_SUCCESS;
		}
		gsl_blas_dgemv(CblasTrans, -1., r, qtf, 0., gradient);
		gsl_vector_div(gradient, diag);
		if (gnorm = gsl_blas_dnrm2(gradient); gnorm == 0.) {
			gsl_vector_set_zero(p);
			gsl_blas_daxpy(delta / qnorm, newton, p);
			return GSL_SUCCESS;
		}
		// minimum_step(gnorm, diag, gradient);
		gsl_vector_scale(gradient, 1. / gnorm);
		gsl_vector_div(gradient, diag);
		// Use p as temporary space to compute Rg
		gsl_blas_dgemv(CblasNoTrans, 1., r, gradient, 0., p);

		temp = gsl_blas_dnrm2(p);
		sgnorm = (gnorm / temp) / temp;
		if (sgnorm > delta) {
			gsl_vector_set_zero(p);
			gsl_blas_daxpy(delta, gradient, p);
			return GSL_SUCCESS;
		}
		bnorm = gsl_blas_dnrm2(qtf);
		{
			long double bg = bnorm / gnorm;
			long double bq = bnorm / qnorm;
			long double dq = delta / qnorm;
			long double dq2 = dq * dq;
			long double sd = sgnorm / delta;
			long double sd2 = sd * sd;

			long double t1 = bg * bq * sd;
			long double u = t1 - dq;
			long double t2 = t1 - dq * sd2 + sqrt(u * u + (1 - dq2) * (1 - sd2));

			long double alpha = dq * (1 - sd2) / t2;
			long double beta = (1 - alpha) * sgnorm;
			gsl_vector_set_zero(p);
			gsl_blas_daxpy(alpha, newton, p);
			gsl_blas_daxpy(beta, gradient, p);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::Dogleg(long double trust_radius, long double newton_norm, long double newton_norm2, long double gradient_norm, long double gradient_norm2, long double newton_gradient_dot, const gsl_vector *newton, const gsl_vector *gradient, gsl_vector *dx) {
		if (newton_norm <= trust_radius) {
			gsl_vector_memcpy(dx, newton);
			return GSL_SUCCESS;
		}
		gsl_vector_set_zero(dx);
		if (gradient_norm >= trust_radius) {
			gsl_blas_daxpy(trust_radius / gradient_norm, gradient, dx);
			return GSL_SUCCESS;
		}
		double d_newton_gradient_dot;
		gsl_blas_ddot(newton, gradient, &d_newton_gradient_dot);
		newton_gradient_dot = d_newton_gradient_dot;
		const long double a = gradient_norm2 + newton_norm2 - 2. * newton_gradient_dot; // a > 0
		const long double b = 2. * (newton_gradient_dot - gradient_norm2);
		const long double c = gradient_norm2 - gsl_pow_2(trust_radius); // c < 0
		long double s_minus_plus[2];
		if (int root_num = PolySolveQuadratic(a, b, c, s_minus_plus); root_num != 2 || s_minus_plus[1] > 1.) {
			gsl_vector_memcpy(dx, gradient); // fallback to gradient
			return GSL_SUCCESS;
		}
		gsl_blas_daxpy(s_minus_plus[1], newton, dx);
		gsl_blas_daxpy(1. - s_minus_plus[1], gradient, dx);
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::HybridAlloc(void *vstate, size_t n) {
		auto state = static_cast<HybridState *>(vstate);
		if (state->x_trial = gsl_vector_calloc(n); state->x_trial == nullptr)
			GSL_ERROR("failed to allocate space for x_trial", GSL_ENOMEM);

		if (state->f_trial = gsl_vector_calloc(n); state->f_trial == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for f_trial", GSL_ENOMEM);
		}
		if (state->jacobian = gsl_matrix_calloc(n, n); state->jacobian == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for jacobian", GSL_ENOMEM);
		}
		if (state->lu = gsl_matrix_calloc(n, n); state->lu == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for lu", GSL_ENOMEM);
		}
		if (state->permutation = gsl_permutation_calloc(n); state->permutation == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for permutation", GSL_ENOMEM);
		}
		if (state->newton = gsl_vector_calloc(n); state->newton == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for newton", GSL_ENOMEM);
		}
		if (state->gradient = gsl_vector_calloc(n); state->gradient == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for gradient", GSL_ENOMEM);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::HybridSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<HybridState *>(vstate);
		for (int status = GSL_MULTIROOT_FN_EVAL(function, x, f); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x, f)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_scale(x, min(0.99999999l, static_cast<long double>(gsl_vector_get(f, 0))));
		}
		if (int status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON, state->jacobian); status != GSL_SUCCESS)
			return status;
		state->trust_radius = gsl_blas_dnrm2(x);
		state->epsilon_coefficient = 1.;
		state->directional = false;
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::HybridIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<HybridState *>(vstate);
		gsl_vector *x_trial = state->x_trial, *f_trial = state->f_trial, *newton = state->newton, *gradient = state->gradient;
		gsl_matrix *jacobian = state->jacobian;
		int signum, status;
		const long double f_norm = gsl_blas_dnrm2(f);
		const long double trust_radius_limit = gsl_blas_dnrm2(x) * GSL_DBL_EPSILON;
		long double dx_norm, newton_norm, newton_norm2, gradient_norm, gradient_norm2, newton_gradient_dot, jacobian_gradient_norm2, t;
		gsl_matrix_memcpy(state->lu, jacobian);
		if (status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
			return status;
		if (status = gsl_linalg_LU_solve(state->lu, state->permutation, f, dx); status != GSL_SUCCESS)
			return status;
		dx_norm = gsl_blas_dnrm2(dx);
		for (int iteration_count = 0; true;) {
			gsl_vector_memcpy(x_trial, x);
			if (dx_norm > state->trust_radius)
				gsl_blas_daxpy(-state->trust_radius / dx_norm, dx, x_trial);
			else
				gsl_blas_daxpy(-1., dx, x_trial);
			if (status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial); status != GSL_SUCCESS) {
				if (status != GSL_EDOM)
					return GSL_EBADFUNC;
				gsl_vector_scale(x_trial, min(0.99999999l, static_cast<long double>(gsl_vector_get(f_trial, 0))));
				gsl_vector_sub(x_trial, x);
				// here x_trial is the difference between x and x_trial
				state->trust_radius = min(0.5l * state->trust_radius, static_cast<long double>(gsl_blas_dnrm2(x_trial)));
				continue;
			}
			if (long double f_trial_norm = gsl_blas_dnrm2(f_trial); f_trial_norm < f_norm) {
				gsl_vector_memcpy(x, x_trial);
				gsl_vector_memcpy(f, f_trial);
				if (dx_norm > state->trust_radius)
					state->trust_radius *= 2.;
				else
					state->trust_radius = 2. * dx_norm;
				break;
			}
			if (state->epsilon_coefficient > GSL_ROOT5_DBL_EPSILON)
				state->epsilon_coefficient *= 0.5;
			if (++iteration_count == 2) {
				if (status = TwoSidedDirectionalJacobian(function, x, dx, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian); status != GSL_SUCCESS)
					return GSL_EBADFUNC;
				gsl_matrix_memcpy(state->lu, jacobian);
				if (status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
					return status;
				if (status = gsl_linalg_LU_solve(state->lu, state->permutation, f, newton); status != GSL_SUCCESS)
					return status;
				gsl_blas_dgemv(CblasTrans, 1., jacobian, f, 0., gradient);
				double d_newton_norm2, d_gradient_norm2, d_jacobian_gradient_norm2;
				gsl_blas_ddot(newton, newton, &d_newton_norm2);
				newton_norm2 = d_newton_norm2;
				gsl_blas_ddot(gradient, gradient, &d_gradient_norm2);
				gradient_norm2 = d_gradient_norm2;
				gsl_blas_dgemv(CblasTrans, 1., jacobian, gradient, 0., x_trial);
				gsl_blas_ddot(x_trial, x_trial, &d_jacobian_gradient_norm2);
				jacobian_gradient_norm2 = d_jacobian_gradient_norm2;
				t = gradient_norm2 / jacobian_gradient_norm2;
				gsl_vector_scale(gradient, t);
				gradient_norm2 *= gsl_pow_2(t);
				newton_norm = sqrt(newton_norm2);
				gradient_norm = sqrt(gradient_norm2);
				double d_newton_gradient_dot;
				gsl_blas_ddot(newton, gradient, &d_newton_gradient_dot);
				newton_gradient_dot = d_newton_gradient_dot;
				Dogleg(state->trust_radius, newton_norm, newton_norm2, gradient_norm, gradient_norm2, newton_gradient_dot, newton, gradient, dx);
				dx_norm = gsl_blas_dnrm2(dx);
				state->directional = true;
			} else {
				state->trust_radius *= 0.5;
				if (state->trust_radius < trust_radius_limit)
					return GSL_ENOPROG;
				if (iteration_count > 2) {
					Dogleg(state->trust_radius, newton_norm, newton_norm2, gradient_norm, gradient_norm2, newton_gradient_dot, newton, gradient, dx);
					dx_norm = gsl_blas_dnrm2(dx);
				}
			}
		}
		if (state->directional)
			status = TwoSidedDirectionalJacobian(function, x, dx, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		else
			status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		if (status != GSL_SUCCESS)
			return GSL_EBADFUNC;
		return GSL_SUCCESS;
	}
	void MultiFunctionSolver::HybridFree(void *vstate) {
		auto state = static_cast<HybridState *>(vstate);
		gsl_vector_free(state->x_trial);
		gsl_vector_free(state->f_trial);
		gsl_matrix_free(state->jacobian);
		gsl_matrix_free(state->lu);
		gsl_permutation_free(state->permutation);
		gsl_vector_free(state->newton);
		gsl_vector_free(state->gradient);
	}
	int MultiFunctionSolver::HybridAdditionAlloc(void *vstate, size_t n) {
		auto state = static_cast<HybridAdditionState *>(vstate);
		if (state->x_trial = gsl_vector_calloc(n); state->x_trial == nullptr)
			GSL_ERROR("failed to allocate space for x_trial", GSL_ENOMEM);
		if (state->x_trial2 = gsl_vector_calloc(n); state->x_trial2 == nullptr) {
			HybridAdditionFree(state);
			GSL_ERROR("failed to allocate space for x_trial2", GSL_ENOMEM);
		}
		if (state->f_trial = gsl_vector_calloc(n); state->f_trial == nullptr) {
			HybridAdditionFree(state);
			GSL_ERROR("failed to allocate space for f_trial", GSL_ENOMEM);
		}
		if (state->f_trial2 = gsl_vector_calloc(n); state->f_trial2 == nullptr) {
			HybridAdditionFree(state);
			GSL_ERROR("failed to allocate space for f_trial2", GSL_ENOMEM);
		}
		if (state->jacobian = gsl_matrix_calloc(n, n); state->jacobian == nullptr) {
			HybridAdditionFree(state);
			GSL_ERROR("failed to allocate space for jacobian", GSL_ENOMEM);
		}
		if (state->lu = gsl_matrix_calloc(n, n); state->lu == nullptr) {
			HybridAdditionFree(state);
			GSL_ERROR("failed to allocate space for lu", GSL_ENOMEM);
		}
		if (state->permutation = gsl_permutation_calloc(n); state->permutation == nullptr) {
			HybridAdditionFree(state);
			GSL_ERROR("failed to allocate space for permutation", GSL_ENOMEM);
		}
		if (state->newton = gsl_vector_calloc(n); state->newton == nullptr) {
			HybridAdditionFree(state);
			GSL_ERROR("failed to allocate space for newton", GSL_ENOMEM);
		}
		if (state->gradient = gsl_vector_calloc(n); state->gradient == nullptr) {
			HybridAdditionFree(state);
			GSL_ERROR("failed to allocate space for gradient", GSL_ENOMEM);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::HybridAdditionSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<HybridAdditionState *>(vstate);
		for (int status = GSL_MULTIROOT_FN_EVAL(function, x, f); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x, f)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_scale(x, min(0.99999999l, static_cast<long double>(gsl_vector_get(f, 0))));
		}
		if (int status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON, state->jacobian); status != GSL_SUCCESS)
			return status;
		state->trust_radius = gsl_blas_dnrm2(x);
		state->epsilon_coefficient = 1.;
		state->gradient_coefficient = 0.;
		state->directional = false;
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::HybridAdditionIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<HybridAdditionState *>(vstate);
		gsl_vector *x_trial = state->x_trial, *x_trial2 = state->x_trial2, *f_trial = state->f_trial, *f_trial2 = state->f_trial2, *newton = state->newton, *gradient = state->gradient;
		gsl_matrix *jacobian = state->jacobian;
		int signum, status;
		long double &trust_radius = state->trust_radius;
		gsl_matrix_memcpy(state->lu, jacobian);
		if (status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
			return status;
		if (status = gsl_linalg_LU_solve(state->lu, state->permutation, f, dx); status != GSL_SUCCESS)
			return status;
		const long double f_norm = gsl_blas_dnrm2(f);
		const long double recalc_radius_limit = gsl_blas_dnrm2(x) * GSL_ROOT4_DBL_EPSILON;
		const long double trust_radius_limit = gsl_blas_dnrm2(x) * GSL_DBL_EPSILON;
		long double newton_norm, gradient_norm = -1., dx_norm = gsl_blas_dnrm2(dx), f_trial_norm, f_trial2_norm, newton_gradient_dot, gradient_coefficient = 0., negative_direction_coefficient, positive_direction_coefficient;
		for (; true; trust_radius *= 0.5) {
			gsl_vector_memcpy(x_trial, x);
			if (dx_norm > trust_radius)
				gsl_blas_daxpy(-trust_radius / dx_norm, dx, x_trial);
			else
				gsl_blas_daxpy(-1., dx, x_trial);
			if (status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial); status != GSL_SUCCESS) {
				if (status != GSL_EDOM)
					return GSL_EBADFUNC;
				gsl_vector_scale(x_trial, min(0.99999999l, static_cast<long double>(gsl_vector_get(f_trial, 0))));
				gsl_vector_sub(x_trial, x);
				// here x_trial is the difference between x and x_trial
				if (long double trust_radius_limit = 2. * gsl_blas_dnrm2(x_trial); trust_radius_limit < trust_radius)
					trust_radius = trust_radius_limit; // there is a *= 0.5
				continue;
			}
			if (f_trial_norm = gsl_blas_dnrm2(f_trial); f_trial_norm < f_norm) {
				gsl_vector_memcpy(x, x_trial);
				gsl_vector_memcpy(f, f_trial);
				if (dx_norm > trust_radius)
					trust_radius *= 2.;
				else
					trust_radius = 2. * dx_norm;
				break;
			}
			if (state->epsilon_coefficient > GSL_ROOT5_DBL_EPSILON)
				state->epsilon_coefficient *= 0.5;
			if (trust_radius < trust_radius_limit)
				return GSL_ENOPROG;
			if (trust_radius < recalc_radius_limit) {
				if (!state->directional) {
					if (status = TwoSidedDirectionalJacobian(function, x, dx, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian); status != GSL_SUCCESS)
						return GSL_EBADFUNC;
					gsl_matrix_memcpy(state->lu, jacobian);
					if (status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
						return status;
					if (status = gsl_linalg_LU_solve(state->lu, state->permutation, f, dx); status != GSL_SUCCESS)
						return status;
					dx_norm = gsl_blas_dnrm2(dx);
					state->directional = true;
				}
				if (gradient_norm == -1.) {
					gsl_vector_memcpy(newton, dx);
					newton_norm = dx_norm;
					gsl_blas_dgemv(CblasTrans, 1., jacobian, f, 0., gradient);
					double d_newton_gradient_dot;
					gsl_blas_ddot(newton, gradient, &d_newton_gradient_dot);
					gsl_blas_daxpy(-d_newton_gradient_dot / newton_norm, newton, gradient); // gradient now is prep to newton.
					gradient_norm = gsl_blas_dnrm2(gradient);
				}
				if (gradient_coefficient == 0.) {
					if (newton_norm < trust_radius)
						positive_direction_coefficient = GSL_SQRT_DBL_EPSILON * newton_norm / gradient_norm;
					else
						positive_direction_coefficient = max(100. * trust_radius_limit, 2. * trust_radius_limit * dx_norm / trust_radius) / gradient_norm;
					negative_direction_coefficient = -positive_direction_coefficient;
				} else if (long double lower_limit = 10. * trust_radius_limit * dx_norm / trust_radius / gradient_norm; positive_direction_coefficient < lower_limit) {
					positive_direction_coefficient = lower_limit;
					negative_direction_coefficient = -positive_direction_coefficient;
				}
				gsl_vector_memcpy(dx, newton);
				gsl_blas_daxpy(gradient_coefficient + negative_direction_coefficient, gradient, dx);
				for (int iteration_limit = 64; iteration_limit > 0 && newton_norm >= gradient_norm * negative_direction_coefficient; --iteration_limit) {
					gsl_vector_memcpy(x_trial2, x);
					if (dx_norm = gsl_blas_dnrm2(dx); dx_norm > trust_radius)
						gsl_blas_daxpy(-trust_radius / dx_norm, dx, x_trial2);
					else
						gsl_blas_daxpy(-1., dx, x_trial2);
					if (status = GSL_MULTIROOT_FN_EVAL(function, x_trial2, f_trial2); status != GSL_SUCCESS)
						break;
					if (f_trial2_norm = gsl_blas_dnrm2(f_trial2); f_trial2_norm >= f_trial_norm)
						break;
					f_trial_norm = f_trial2_norm;
					gsl_vector_memcpy(x_trial, x_trial2);
					gsl_vector_memcpy(f_trial, f_trial2);
					gsl_blas_daxpy(negative_direction_coefficient, gradient, dx);
					if (f_norm < f_trial_norm)
						negative_direction_coefficient *= 2.;
					else if (gradient_coefficient == 0.)
						negative_direction_coefficient *= 1.5;
					else
						negative_direction_coefficient *= 0.9;
				}
				if (positive_direction_coefficient != -negative_direction_coefficient) {
					if (f_trial_norm < f_norm) {
						gsl_vector_memcpy(x, x_trial);
						gsl_vector_memcpy(f, f_trial);
						break;
					}
					negative_direction_coefficient *= 0.1;
					positive_direction_coefficient = -negative_direction_coefficient;
					gradient_coefficient += negative_direction_coefficient;
					gsl_blas_daxpy(-negative_direction_coefficient, gradient, dx);
					continue;
				}
				gsl_vector_memcpy(dx, newton);
				gsl_blas_daxpy(gradient_coefficient + positive_direction_coefficient, gradient, dx);
				for (int iteration_limit = 64; iteration_limit > 0 && newton_norm >= gradient_norm * positive_direction_coefficient; --iteration_limit) {
					gsl_vector_memcpy(x_trial2, x);
					if (dx_norm = gsl_blas_dnrm2(dx); dx_norm > trust_radius)
						gsl_blas_daxpy(-trust_radius / dx_norm, dx, x_trial2);
					else
						gsl_blas_daxpy(-1., dx, x_trial2);
					if (status = GSL_MULTIROOT_FN_EVAL(function, x_trial2, f_trial2); status != GSL_SUCCESS)
						break;
					if (f_trial2_norm = gsl_blas_dnrm2(f_trial2); f_trial2_norm >= f_trial_norm)
						break;
					f_trial_norm = f_trial2_norm;
					gsl_vector_memcpy(x_trial, x_trial2);
					gsl_vector_memcpy(f_trial, f_trial2);
					gsl_blas_daxpy(positive_direction_coefficient, gradient, dx);
					if (f_norm < f_trial_norm)
						positive_direction_coefficient *= 2.;
					else if (gradient_coefficient == 0.)
						positive_direction_coefficient *= 1.5;
					else
						positive_direction_coefficient *= 0.9;
				}
				if (positive_direction_coefficient != -negative_direction_coefficient) {
					if (f_trial_norm < f_norm) {
						gsl_vector_memcpy(x, x_trial);
						gsl_vector_memcpy(f, f_trial);
						break;
					}
					positive_direction_coefficient *= 0.1;
					negative_direction_coefficient = -positive_direction_coefficient;
					gradient_coefficient += positive_direction_coefficient;
					gsl_blas_daxpy(-positive_direction_coefficient, gradient, dx);
					continue;
				}
				gsl_vector_memcpy(dx, newton);
				gsl_blas_daxpy(gradient_coefficient, gradient, dx);
				dx_norm = gsl_blas_dnrm2(dx);
			}
		}
		if (state->directional)
			status = TwoSidedDirectionalJacobian(function, x, dx, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		else
			status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		if (status != GSL_SUCCESS)
			return GSL_EBADFUNC;
		return GSL_SUCCESS;
	}
	void MultiFunctionSolver::HybridAdditionFree(void *vstate) {
		auto state = static_cast<HybridAdditionState *>(vstate);
		gsl_vector_free(state->x_trial);
		gsl_vector_free(state->x_trial2);
		gsl_vector_free(state->f_trial);
		gsl_vector_free(state->f_trial2);
		gsl_matrix_free(state->jacobian);
		gsl_matrix_free(state->lu);
		gsl_permutation_free(state->permutation);
		gsl_vector_free(state->newton);
		gsl_vector_free(state->gradient);
	}
	int MultiFunctionSolver::DNewtonAlloc(void *vstate, size_t n) {
		auto state = static_cast<DNewtonState *>(vstate);
		if (state->x_trial = gsl_vector_calloc(n); state->x_trial == nullptr)
			GSL_ERROR("failed to allocate space for x_trial", GSL_ENOMEM);
		if (state->f_trial = gsl_vector_calloc(n); state->f_trial == nullptr) {
			DNewtonFree(vstate);
			GSL_ERROR("failed to allocate space for f_trial", GSL_ENOMEM);
		}
		if (state->lu = gsl_matrix_calloc(n, n); state->lu == nullptr) {
			DNewtonFree(vstate);
			GSL_ERROR("failed to allocate space for lu", GSL_ENOMEM);
		}
		if (state->permutation = gsl_permutation_calloc(n); state->permutation == nullptr) {
			DNewtonFree(vstate);
			GSL_ERROR("failed to allocate space for permutation", GSL_ENOMEM);
		}
		if (state->jacobian = gsl_matrix_calloc(n, n); state->jacobian == nullptr) {
			DNewtonFree(vstate);
			GSL_ERROR("failed to allocate space for d", GSL_ENOMEM);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::DNewtonSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<DNewtonState *>(vstate);
		for (int status = GSL_MULTIROOT_FN_EVAL(function, x, f); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x, f)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_scale(x, min(0.99999999l, static_cast<long double>(gsl_vector_get(f, 0))));
		}
		state->trust_radius = gsl_blas_dnrm2(x);
		state->epsilon_coefficient = 1.;
		state->directional = false;
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::DNewtonIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<DNewtonState *>(vstate);
		gsl_vector *x_trial = state->x_trial, *f_trial = state->f_trial;
		gsl_matrix *jacobian = state->jacobian;
		int signum, status;
		long double &trust_radius = state->trust_radius;
		if (state->directional)
			status = TwoSidedDirectionalJacobian(function, x, dx, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		else
			status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		if (status != GSL_SUCCESS)
			return GSL_EBADFUNC;
		gsl_matrix_memcpy(state->lu, jacobian);
		if (status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
			return status;
		if (status = gsl_linalg_LU_solve(state->lu, state->permutation, f, dx); status != GSL_SUCCESS)
			return status;
		const long double x_norm = gsl_blas_dnrm2(x), f_norm = gsl_blas_dnrm2(f);
		const long double trust_radius_limit = max(1.l, x_norm) * GSL_DBL_EPSILON;
		long double dx_norm = gsl_blas_dnrm2(dx);
		while (true) {
			gsl_vector_memcpy(x_trial, x);
			if (dx_norm > trust_radius)
				gsl_blas_daxpy(-trust_radius / dx_norm, dx, x_trial);
			else
				gsl_blas_daxpy(-1., dx, x_trial);
			if (status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial); status != GSL_SUCCESS) {
				if (status != GSL_EDOM)
					return GSL_EBADFUNC;
				gsl_vector_scale(x_trial, min(0.99999999l, static_cast<long double>(gsl_vector_get(f_trial, 0))));
				gsl_vector_sub(x_trial, x);
				// here x_trial is the difference between x and x_trial
				trust_radius = min(0.5l * trust_radius, static_cast<long double>(gsl_blas_dnrm2(x_trial)));
				continue;
			}
			if (long double f_trial_norm = gsl_blas_dnrm2(f_trial); f_trial_norm < f_norm) {
				gsl_vector_memcpy(x, x_trial);
				gsl_vector_memcpy(f, f_trial);
				if (dx_norm > trust_radius)
					trust_radius *= 2.;
				else
					trust_radius = 2. * dx_norm;
				return GSL_SUCCESS;
			}
			if (!state->directional) {
				state->directional = true;
				return GSL_SUCCESS;
			}
			if (state->epsilon_coefficient > GSL_ROOT5_DBL_EPSILON)
				state->epsilon_coefficient *= 0.5;
			if (trust_radius *= 0.5; trust_radius < trust_radius_limit)
				return GSL_ENOPROG;
		}
	}
	void MultiFunctionSolver::DNewtonFree(void *vstate) {
		auto state = static_cast<DNewtonState *>(vstate);
		gsl_vector_free(state->x_trial);
		gsl_vector_free(state->f_trial);
		gsl_matrix_free(state->jacobian);
		gsl_matrix_free(state->lu);
		gsl_permutation_free(state->permutation);
	}
	int MultiFunctionSolver::DNewtonRotationTranslationAlloc(void *vstate, size_t n) {
		auto state = static_cast<DNewtonRotationTranslationState *>(vstate);
		if (state->x_trial = gsl_vector_calloc(n); state->x_trial == nullptr)
			GSL_ERROR("failed to allocate space for x_trial", GSL_ENOMEM);
		if (state->dx_translation = gsl_vector_calloc(n); state->dx_translation == nullptr) {
			DNewtonRotationTranslationFree(vstate);
			GSL_ERROR("failed to allocate space for dx_translation", GSL_ENOMEM);
		}
		if (state->f_trial = gsl_vector_calloc(n); state->f_trial == nullptr) {
			DNewtonRotationTranslationFree(vstate);
			GSL_ERROR("failed to allocate space for f_trial", GSL_ENOMEM);
		}
		if (state->lu = gsl_matrix_calloc(n, n); state->lu == nullptr) {
			DNewtonRotationTranslationFree(vstate);
			GSL_ERROR("failed to allocate space for lu", GSL_ENOMEM);
		}
		if (state->permutation = gsl_permutation_calloc(n); state->permutation == nullptr) {
			DNewtonRotationTranslationFree(vstate);
			GSL_ERROR("failed to allocate space for permutation", GSL_ENOMEM);
		}
		if (state->jacobian = gsl_matrix_calloc(n, n); state->jacobian == nullptr) {
			DNewtonRotationTranslationFree(vstate);
			GSL_ERROR("failed to allocate space for d", GSL_ENOMEM);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::DNewtonRotationTranslationSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<DNewtonRotationTranslationState *>(vstate);
		if (int status = ScaleX(function, x, f); status != GSL_SUCCESS)
			return status;
		state->trust_radius = max(1., gsl_blas_dnrm2(x));
		state->epsilon_coefficient = 1.;
		state->directional = false;
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::DNewtonRotationIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<DNewtonRotationTranslationState *>(vstate);
		gsl_vector *x_trial = state->x_trial, *f_trial = state->f_trial;
		gsl_matrix *jacobian = state->jacobian;
		int signum, status;
		long double &trust_radius = state->trust_radius;
		if (state->directional)
			status = TwoSidedDirectionalJacobian(function, x, dx, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		else
			status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		if (status != GSL_SUCCESS)
			return GSL_EBADFUNC;
		gsl_matrix_memcpy(state->lu, jacobian);
		if (status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
			return status;
		if (status = gsl_linalg_LU_solve(state->lu, state->permutation, f, dx); status != GSL_SUCCESS)
			return status;
		const long double x_norm = gsl_blas_dnrm2(x), f_norm = gsl_blas_dnrm2(f);
		const long double rotation_radius_limit = max(1.l, x_norm) * GSL_ROOT3_DBL_EPSILON;
		const long double trust_radius_limit = max(1.l, x_norm) * GSL_DBL_EPSILON;
		long double dx_norm = gsl_blas_dnrm2(dx);
		while (true) {
			gsl_vector_memcpy(x_trial, x);
			if (dx_norm > trust_radius)
				gsl_blas_daxpy(-trust_radius / dx_norm, dx, x_trial);
			else
				gsl_blas_daxpy(-1., dx, x_trial);
			if (status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial); status != GSL_SUCCESS) {
				if (status != GSL_EDOM)
					return GSL_EBADFUNC;
				gsl_vector_scale(x_trial, gsl_vector_get(f_trial, 0));
				gsl_vector_sub(x_trial, x);
				// here x_trial is the difference between x and x_trial
				trust_radius = min(0.5l * trust_radius, static_cast<long double>(gsl_blas_dnrm2(x_trial)));
				continue;
			}
			if (long double f_trial_norm = gsl_blas_dnrm2(f_trial); f_trial_norm < f_norm) {
				gsl_vector_memcpy(x, x_trial);
				gsl_vector_memcpy(f, f_trial);
				if (dx_norm > trust_radius)
					trust_radius *= 2.;
				else
					trust_radius = 2. * dx_norm;
				return GSL_SUCCESS;
			}
			if (!state->directional) {
				state->directional = true;
				return GSL_SUCCESS;
			}
			if (state->epsilon_coefficient > GSL_ROOT5_DBL_EPSILON)
				state->epsilon_coefficient *= 0.5;
			if (trust_radius < rotation_radius_limit) {
				long double projected_x, projected_y;
				if (state->trace_to_plane) {
					projected_x = state->projected_x + gsl_vector_get(f, 0);
					projected_y = state->projected_y + gsl_vector_get(f, 1);
				} else {
					long double cos_theta_photon = state->cos_theta_obj + gsl_vector_get(f, 0);
					long double theta_photon = acos(cos_theta_photon);
					long double phi_photon = state->phi_obj - gsl_vector_get(f, 1);
					long double sin_theta_photon = sin(theta_photon);
					projected_x = sin_theta_photon * sin(phi_photon);
					projected_y = cos_theta_photon * state->sin_theta_obs - sin_theta_photon * cos(phi_photon) * state->cos_theta_obs;
				}
				long double iota_photon = atan2(projected_y, projected_x);
				long double cos_delta_iota = cos(state->iota_obj - iota_photon);
				long double sin_delta_iota = sin(state->iota_obj - iota_photon);
				gsl_vector_set(dx, 0, gsl_vector_get(x, 0) * (1. - cos_delta_iota) + gsl_vector_get(x, 1) * sin_delta_iota);
				gsl_vector_set(dx, 1, -gsl_vector_get(x, 0) * sin_delta_iota + gsl_vector_get(x, 1) * (1. - cos_delta_iota));
				dx_norm = gsl_blas_dnrm2(dx);
				gsl_vector_memcpy(x_trial, x);
				if (dx_norm * GSL_ROOT6_DBL_EPSILON > trust_radius)
					gsl_blas_daxpy(-GSL_ROOT6_DBL_EPSILON, dx, x_trial);
				else if (dx_norm > trust_radius)
					gsl_blas_daxpy(-trust_radius / dx_norm, dx, x_trial);
				else
					gsl_blas_daxpy(-1., dx, x_trial);
				if (status = ScaleX(function, x_trial, f_trial); status != GSL_SUCCESS)
					return status;
				gsl_vector_memcpy(x, x_trial);
				gsl_vector_memcpy(f, f_trial);
				if (long double f_trial_norm = gsl_blas_dnrm2(f_trial); f_trial_norm < f_norm)
					trust_radius *= 2.;
				else if (trust_radius *= 0.5; trust_radius < trust_radius_limit)
					return GSL_ENOPROG;
				return GSL_SUCCESS;
			}
			trust_radius *= 0.5;
		}
	}
	int MultiFunctionSolver::DNewtonTranslationIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<DNewtonRotationTranslationState *>(vstate);
		gsl_vector *x_trial = state->x_trial, *f_trial = state->f_trial, *dx_translation = state->dx_translation;
		gsl_matrix *jacobian = state->jacobian;
		int signum, status;
		long double &trust_radius = state->trust_radius;
		if (state->directional)
			status = TwoSidedDirectionalJacobian(function, x, dx, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		else
			status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		if (status != GSL_SUCCESS)
			return GSL_EBADFUNC;
		gsl_matrix_memcpy(state->lu, jacobian);
		if (status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
			return status;
		if (status = gsl_linalg_LU_solve(state->lu, state->permutation, f, dx); status != GSL_SUCCESS)
			return status;
		const long double f_norm = gsl_blas_dnrm2(f);
		const long double translation_radius_limit = gsl_blas_dnrm2(x) * GSL_ROOT3_DBL_EPSILON;
		const long double trust_radius_limit = gsl_blas_dnrm2(x) * GSL_DBL_EPSILON;
		long double dx_norm = gsl_blas_dnrm2(dx), dx_translation_norm;
		while (true) {
			gsl_vector_memcpy(x_trial, x);
			if (dx_norm > trust_radius)
				gsl_blas_daxpy(-trust_radius / dx_norm, dx, x_trial);
			else
				gsl_blas_daxpy(-1., dx, x_trial);
			if (status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial); status != GSL_SUCCESS) {
				if (status != GSL_EDOM)
					return GSL_EBADFUNC;
				gsl_vector_scale(x_trial, gsl_vector_get(f_trial, 0));
				gsl_vector_sub(x_trial, x);
				// here x_trial is the difference between x and x_trial
				trust_radius = min(0.5l * trust_radius, static_cast<long double>(gsl_blas_dnrm2(x_trial)));
				continue;
			}
			if (long double f_trial_norm = gsl_blas_dnrm2(f_trial); f_trial_norm < f_norm) {
				gsl_vector_memcpy(x, x_trial);
				gsl_vector_memcpy(f, f_trial);
				if (dx_norm > trust_radius)
					trust_radius *= 2.;
				else
					trust_radius = 2. * dx_norm;
				return GSL_SUCCESS;
			}
			if (!state->directional) {
				state->directional = true;
				return GSL_SUCCESS;
			}
			if (state->epsilon_coefficient > GSL_ROOT5_DBL_EPSILON)
				state->epsilon_coefficient *= 0.5;
			if (trust_radius < translation_radius_limit) {
				if (state->trace_to_plane) {
					gsl_vector_memcpy(dx_translation, f);
				} else {
					long double cos_theta_photon = state->cos_theta_obj + gsl_vector_get(f, 0);
					long double theta_photon = acos(cos_theta_photon);
					long double phi_photon = state->phi_obj - gsl_vector_get(f, 1);
					long double sin_theta_photon = sin(theta_photon);
					long double projected_x = sin_theta_photon * sin(phi_photon);
					long double projected_y = cos_theta_photon * state->sin_theta_obs - sin_theta_photon * cos(phi_photon) * state->cos_theta_obs;
					gsl_vector_set(dx_translation, 0, state->r_obj * (projected_x - state->projected_x));
					gsl_vector_set(dx_translation, 1, state->r_obj * (projected_y - state->projected_y));
				}
				dx_translation_norm = gsl_blas_dnrm2(dx_translation);
				gsl_vector_memcpy(x_trial, x);
				if (dx_translation_norm * GSL_ROOT6_DBL_EPSILON > trust_radius)
					gsl_blas_daxpy(-GSL_ROOT6_DBL_EPSILON, dx_translation, x_trial);
				else if (dx_translation_norm > trust_radius)
					gsl_blas_daxpy(-trust_radius / dx_translation_norm, dx_translation, x_trial);
				else
					gsl_blas_daxpy(-1., dx_translation, x_trial);
				if (status = ScaleX(function, x_trial, f_trial); status != GSL_SUCCESS)
					return status;
				gsl_vector_memcpy(x, x_trial);
				gsl_vector_memcpy(f, f_trial);
				if (long double f_trial_norm = gsl_blas_dnrm2(f_trial); f_trial_norm < f_norm)
					trust_radius *= 2.;
				else if (trust_radius *= 0.5; trust_radius < trust_radius_limit)
					return GSL_ENOPROG;
				return GSL_SUCCESS;
			}
			trust_radius *= 0.5;
		}
		return GSL_SUCCESS;
	}
	void MultiFunctionSolver::DNewtonRotationTranslationFree(void *vstate) {
		auto state = static_cast<DNewtonRotationTranslationState *>(vstate);
		gsl_vector_free(state->x_trial);
		gsl_vector_free(state->dx_translation);
		gsl_vector_free(state->f_trial);
		gsl_matrix_free(state->jacobian);
		gsl_matrix_free(state->lu);
		gsl_permutation_free(state->permutation);
	}
	int MultiFunctionSolver::D2NewtonAlloc(void *vstate, size_t n) {
		auto state = static_cast<D2NewtonState *>(vstate);
		if (state->x_trial = gsl_vector_calloc(n); state->x_trial == nullptr)
			GSL_ERROR("failed to allocate space for x_trial", GSL_ENOMEM);
		if (state->f_trial = gsl_vector_calloc(n); state->f_trial == nullptr) {
			DNewtonFree(vstate);
			GSL_ERROR("failed to allocate space for f_trial", GSL_ENOMEM);
		}
		if (state->jacobian = gsl_vector_calloc(n); state->jacobian == nullptr) {
			DNewtonFree(vstate);
			GSL_ERROR("failed to allocate space for jacobian", GSL_ENOMEM);
		}
		if (state->hessian = gsl_matrix_calloc(n, n); state->hessian == nullptr) {
			DNewtonFree(vstate);
			GSL_ERROR("failed to allocate space for hessian", GSL_ENOMEM);
		}
		if (state->lu = gsl_matrix_calloc(n, n); state->lu == nullptr) {
			DNewtonFree(vstate);
			GSL_ERROR("failed to allocate space for lu", GSL_ENOMEM);
		}
		if (state->permutation = gsl_permutation_calloc(n); state->permutation == nullptr) {
			DNewtonFree(vstate);
			GSL_ERROR("failed to allocate space for permutation", GSL_ENOMEM);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::D2NewtonSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<D2NewtonState *>(vstate);
		for (int status = GSL_MULTIROOT_FN_EVAL(function, x, f); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x, f)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_scale(x, gsl_vector_get(f, 0));
		}
		if (int status = OneSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON, state->hessian); status != GSL_SUCCESS)
			return status;
		state->trust_radius = gsl_blas_dnrm2(x);
		gsl_vector_set_zero(dx);
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::D2NewtonIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<D2NewtonState *>(vstate);
		int signum;
		gsl_matrix_memcpy(state->lu, state->hessian);
		if (int status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
			return status;
		// gsl_vector_memcpy(state->dx, dx);
		if (int status = gsl_linalg_LU_solve(state->lu, state->permutation, f, dx); status != GSL_SUCCESS)
			return status;
		const long double f_norm = gsl_blas_dnrm2(f);
		const long double trust_radius_limit = gsl_blas_dnrm2(x) * GSL_SQRT_DBL_EPSILON;
		long double dx_norm = gsl_blas_dnrm2(dx);
		int limit_iteration = 0;
		while (true) {
			gsl_vector_memcpy(state->x_trial, x);
			if (dx_norm > state->trust_radius)
				gsl_blas_daxpy(-state->trust_radius / dx_norm, dx, state->x_trial);
			else
				gsl_blas_daxpy(-1., dx, state->x_trial);
			if (int status = GSL_MULTIROOT_FN_EVAL(function, state->x_trial, state->f_trial); status != GSL_SUCCESS) {
				if (status != GSL_EDOM)
					return GSL_EBADFUNC;
				gsl_vector_scale(state->x_trial, min(0.99999999l, static_cast<long double>(gsl_vector_get(state->f_trial, 0))));
				gsl_vector_sub(state->x_trial, x);
				// here x_trial is the difference between x and x_trial
				if (state->trust_radius = min(0.5 * state->trust_radius, static_cast<long double>(gsl_blas_dnrm2(state->x_trial))); state->trust_radius < trust_radius_limit) {
					if (limit_iteration == 20)
						return GSL_ENOPROG;
					else if (limit_iteration == 0) {
						if (int status = Hessian(function, x, f, GSL_SQRT_DBL_EPSILON, state->jacobian, state->hessian); status != GSL_SUCCESS)
							return GSL_EBADFUNC;
						gsl_matrix_memcpy(state->lu, state->hessian);
						if (int status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
							return status;
						// gsl_vector_memcpy(state->dx, dx);
						if (int status = gsl_linalg_LU_solve(state->lu, state->permutation, state->jacobian, dx); status != GSL_SUCCESS)
							return status;
						dx_norm = gsl_blas_dnrm2(dx);
						state->trust_radius *= 2.;
					}
					++limit_iteration;
				}
				continue;
			}
			if (long double f_trial_norm = gsl_blas_dnrm2(state->f_trial); f_trial_norm < f_norm) {
				gsl_vector_memcpy(x, state->x_trial);
				gsl_vector_memcpy(f, state->f_trial);
				if (dx_norm > state->trust_radius)
					state->trust_radius *= 2.;
				else
					state->trust_radius = 2. * dx_norm;
				break;
			} else {
				if (state->trust_radius *= 0.5; state->trust_radius < trust_radius_limit) {
					if (limit_iteration == 20)
						return GSL_ENOPROG;
					else if (limit_iteration == 0) {
						if (int status = Hessian(function, x, f, GSL_SQRT_DBL_EPSILON, state->f_trial, state->hessian); status != GSL_SUCCESS)
							return GSL_EBADFUNC;
						gsl_matrix_memcpy(state->lu, state->hessian);
						if (int status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
							return status;
						// gsl_vector_memcpy(state->dx, dx);
						if (int status = gsl_linalg_LU_solve(state->lu, state->permutation, state->f_trial, state->x_trial); status != GSL_SUCCESS)
							return status;
						gsl_vector_add(dx, state->x_trial);
						dx_norm = gsl_blas_dnrm2(dx);
						state->trust_radius *= 2.;
					}
					++limit_iteration;
				}
			}
		}
		if (int status = OneSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON, state->hessian); status != GSL_SUCCESS)
			return GSL_EBADFUNC;
		return GSL_SUCCESS;
	}
	void MultiFunctionSolver::D2NewtonFree(void *vstate) {
		auto state = static_cast<D2NewtonState *>(vstate);
		gsl_vector_free(state->x_trial);
		gsl_vector_free(state->f_trial);
		gsl_vector_free(state->jacobian);
		gsl_matrix_free(state->hessian);
		gsl_matrix_free(state->lu);
		gsl_permutation_free(state->permutation);
	}

	int MultiFunctionSolver::GradientIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<DNewtonState *>(vstate);
		gsl_vector *x_trial = state->x_trial, *f_trial = state->f_trial;
		gsl_matrix *jacobian = state->jacobian;
		int signum, status;
		gsl_matrix_memcpy(state->lu, jacobian);
		if (status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
			return status;
		if (status = gsl_linalg_LU_solve(state->lu, state->permutation, f, dx); status != GSL_SUCCESS)
			return status;
		const long double f_norm = gsl_blas_dnrm2(f);
		const long double trust_radius_limit = gsl_blas_dnrm2(x) * GSL_DBL_EPSILON;
		long double dx_norm = gsl_blas_dnrm2(dx);
		int iteration_count = 0;
		while (true) {
			gsl_vector_memcpy(x_trial, x);
			if (dx_norm > state->trust_radius)
				gsl_blas_daxpy(-state->trust_radius / dx_norm, dx, x_trial);
			else
				gsl_blas_daxpy(-1., dx, x_trial);
			if (status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial); status != GSL_SUCCESS) {
				if (status != GSL_EDOM)
					return GSL_EBADFUNC;
				gsl_vector_scale(x_trial, min(0.99999999l, static_cast<long double>(gsl_vector_get(f_trial, 0))));
				gsl_vector_sub(x_trial, x);
				// here x_trial is the difference between x and x_trial
				state->trust_radius = min(0.5 * state->trust_radius, static_cast<long double>(gsl_blas_dnrm2(x_trial)));
				continue;
			}
			if (long double f_trial_norm = gsl_blas_dnrm2(f_trial); f_trial_norm < f_norm) {
				gsl_vector_memcpy(x, x_trial);
				gsl_vector_memcpy(f, f_trial);
				if (dx_norm > state->trust_radius)
					state->trust_radius *= 2.;
				else
					state->trust_radius = 2. * dx_norm;
				break;
			}
			++iteration_count;
			if (state->epsilon_coefficient > GSL_ROOT5_DBL_EPSILON)
				state->epsilon_coefficient *= 0.5;
			if (!state->directional && iteration_count == 2) {
				state->trust_radius *= 2.;
				if (status = TwoSidedDirectionalJacobian(function, x, dx, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian); status != GSL_SUCCESS)
					return GSL_EBADFUNC;
				gsl_blas_dgemv(CblasTrans, 1., jacobian, f, 0., dx);
				dx_norm = gsl_blas_dnrm2(dx);
				state->directional = true;
			} else {
				state->trust_radius *= 0.5;
				if (state->trust_radius < trust_radius_limit)
					return GSL_ENOPROG;
			}
		}
		if (state->directional)
			status = TwoSidedDirectionalJacobian(function, x, dx, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		else
			status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON * state->epsilon_coefficient, jacobian);
		if (status != GSL_SUCCESS)
			return GSL_EBADFUNC;
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::ConjugateGradientAlloc(void *vstate, size_t n) {
		auto state = static_cast<ConjugateGradientState *>(vstate);
		if (state->x_trial = gsl_vector_calloc(n); state->x_trial == nullptr)
			GSL_ERROR("failed to allocate space for x_trial", GSL_ENOMEM);
		if (state->f_trial = gsl_vector_calloc(n); state->f_trial == nullptr) {
			ConjugateGradientFree(vstate);
			GSL_ERROR("failed to allocate space for f_trial", GSL_ENOMEM);
		}
		if (state->last_dx = gsl_vector_calloc(n); state->last_dx == nullptr) {
			ConjugateGradientFree(vstate);
			GSL_ERROR("failed to allocate space for jacobian", GSL_ENOMEM);
		}
		if (state->jacobian = gsl_matrix_calloc(n, n); state->jacobian == nullptr) {
			ConjugateGradientFree(vstate);
			GSL_ERROR("failed to allocate space for jacobian", GSL_ENOMEM);
		}
		if (state->lu = gsl_matrix_calloc(n, n); state->lu == nullptr) {
			ConjugateGradientFree(vstate);
			GSL_ERROR("failed to allocate space for lu", GSL_ENOMEM);
		}
		if (state->permutation = gsl_permutation_calloc(n); state->permutation == nullptr) {
			ConjugateGradientFree(vstate);
			GSL_ERROR("failed to allocate space for permutation", GSL_ENOMEM);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::ConjugateGradientSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<ConjugateGradientState *>(vstate);
		state->iteration_coefficient = 1.0;
		for (int status = GSL_MULTIROOT_FN_EVAL(function, x, f); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x, f)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_scale(x, gsl_vector_get(f, 0));
		}
		if (int status = OneSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON, state->jacobian); status != GSL_SUCCESS)
			return status;
		state->gradient_norm = 1.;
		state->trust_radius = state->iteration_coefficient * gsl_blas_dnrm2(x);
		gsl_vector_set_zero(dx);
		gsl_vector_set_zero(state->last_dx);
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::ConjugateGradientIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<ConjugateGradientState *>(vstate);
		const long double f_norm = gsl_blas_dnrm2(f);
		const long double trust_radius_limit = gsl_blas_dnrm2(x) * GSL_SQRT_DBL_EPSILON;
		if (int status = gsl_blas_dgemv(CblasTrans, 1., state->jacobian, f, 0., state->x_trial); status != GSL_SUCCESS)
			return status;
		long double gradient_norm = gsl_blas_dnrm2(state->x_trial);
		long double step_beta = gradient_norm / state->gradient_norm;
		state->gradient_norm = gradient_norm;
		gsl_vector_memcpy(dx, state->x_trial);
		gsl_blas_daxpy(step_beta, state->last_dx, dx);
		long double dx_norm = gsl_blas_dnrm2(dx);
		int limit_iteration = 0;
		while (true) {
			gsl_vector_memcpy(state->x_trial, x);
			if (dx_norm > state->trust_radius)
				gsl_blas_daxpy(-state->trust_radius / dx_norm, dx, state->x_trial);
			else
				gsl_blas_daxpy(-1., dx, state->x_trial);
			if (int status = GSL_MULTIROOT_FN_EVAL(function, state->x_trial, state->f_trial); status != GSL_SUCCESS) {
				if (status != GSL_EDOM)
					return GSL_EBADFUNC;
				gsl_vector_scale(state->x_trial, min(0.99999999l, static_cast<long double>(gsl_vector_get(state->f_trial, 0))));
				gsl_vector_sub(state->x_trial, x);
				// here x_trial is the difference between x and x_trial
				if (state->trust_radius = min(0.5 * state->trust_radius, static_cast<long double>(gsl_blas_dnrm2(state->x_trial))); state->trust_radius < trust_radius_limit) {
					if (int status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON, state->jacobian); status != GSL_SUCCESS)
						return GSL_EBADFUNC;
					state->trust_radius *= 2.;
				}
				continue;
			}
			if (long double f_trial_norm = gsl_blas_dnrm2(state->f_trial); f_trial_norm < f_norm) {
				gsl_vector_memcpy(x, state->x_trial);
				gsl_vector_memcpy(f, state->f_trial);
				gsl_vector_memcpy(state->last_dx, dx);
				if (dx_norm > state->trust_radius)
					state->trust_radius *= 2.;
				else
					state->trust_radius = 2. * dx_norm;
				break;
			} else {
				if (state->trust_radius *= 0.5; state->trust_radius < trust_radius_limit) {
					if (limit_iteration == 20)
						return GSL_ENOPROG;
					else if (limit_iteration == 0) {
						if (int status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON, state->jacobian); status != GSL_SUCCESS)
							return GSL_EBADFUNC;
						if (int status = gsl_blas_dgemv(CblasTrans, 1., state->jacobian, f, 0., state->x_trial); status != GSL_SUCCESS)
							return status;
						gsl_vector_memcpy(dx, state->x_trial);
						gsl_blas_daxpy(step_beta, state->last_dx, dx);
						dx_norm = gsl_blas_dnrm2(dx);
						state->trust_radius *= 2.;
					}
					++limit_iteration;
				}
			}
		}
		if (int status = TwoSidedJacobian(function, x, f, GSL_SQRT_DBL_EPSILON, state->jacobian); status != GSL_SUCCESS)
			return GSL_EBADFUNC;
		return GSL_SUCCESS;
	}
	void MultiFunctionSolver::ConjugateGradientFree(void *vstate) {
		auto state = static_cast<ConjugateGradientState *>(vstate);
		gsl_vector_free(state->x_trial);
		gsl_vector_free(state->f_trial);
		gsl_vector_free(state->last_dx);
		gsl_matrix_free(state->jacobian);
		gsl_matrix_free(state->lu);
		gsl_permutation_free(state->permutation);
	}
	int MultiFunctionSolver::TriangleAlloc(void *vstate, size_t n) {
		auto state = static_cast<TriangleState *>(vstate);
		if (state->x_trial = gsl_vector_calloc(n); state->x_trial == nullptr)
			GSL_ERROR("failed to allocate space for x_trial", GSL_ENOMEM);
		if (state->f_trial = gsl_vector_calloc(n); state->f_trial == nullptr) {
			TriangleFree(vstate);
			GSL_ERROR("failed to allocate space for f_trial", GSL_ENOMEM);
		}
		if (state->x_a = gsl_vector_calloc(n); state->x_a == nullptr) {
			TriangleFree(vstate);
			GSL_ERROR("failed to allocate space for x_a", GSL_ENOMEM);
		}
		if (state->x_b = gsl_vector_calloc(n); state->x_b == nullptr) {
			TriangleFree(vstate);
			GSL_ERROR("failed to allocate space for x_b", GSL_ENOMEM);
		}
		if (state->f_a = gsl_vector_calloc(n); state->f_a == nullptr) {
			TriangleFree(vstate);
			GSL_ERROR("failed to allocate space for f_a", GSL_ENOMEM);
		}
		if (state->f_b = gsl_vector_calloc(n); state->f_b == nullptr) {
			TriangleFree(vstate);
			GSL_ERROR("failed to allocate space for f_b", GSL_ENOMEM);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::TriangleSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<TriangleState *>(vstate);
		gsl_vector *x_a = state->x_a, *x_b = state->x_b, *f_a = state->f_a, *f_b = state->f_b;
		int status;
		for (status = GSL_MULTIROOT_FN_EVAL(function, x, f); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x, f)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_scale(x, min(0.99999999l, static_cast<long double>(gsl_vector_get(f, 0))));
		}
		// long double dx = -gsl_vector_get(x, 0) * GSL_SQRT_DBL_EPSILON, dy = gsl_vector_get(x, 1) * GSL_SQRT_DBL_EPSILON;
		gsl_vector_memcpy(x_a, x);
		gsl_vector_memcpy(x_b, x);
		gsl_vector_set(x_a, 0, gsl_vector_get(x, 0) * (1. - 0.01));
		gsl_vector_set(x_b, 1, gsl_vector_get(x, 1) * (1. - 0.01));
		for (int iteration = 0; iteration < 10; ++iteration) {
			for (status = GSL_MULTIROOT_FN_EVAL(function, x_a, f_a); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x_a, f_a)) {
				if (status != GSL_EDOM)
					return GSL_EBADFUNC;
				gsl_vector_scale(x_a, min(0.99999999l, static_cast<long double>(gsl_vector_get(f_a, 0))));
			}
			for (status = GSL_MULTIROOT_FN_EVAL(function, x_b, f_b); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x_b, f_b)) {
				if (status != GSL_EDOM)
					return GSL_EBADFUNC;
				gsl_vector_scale(x_b, min(0.99999999l, static_cast<long double>(gsl_vector_get(f_b, 0))));
			}
			if (PointInTriangle(f, f_a, f_b))
				return GSL_SUCCESS;
			gsl_vector_sub(f_a, f);
			gsl_vector_sub(f_b, f);
			long double x_a_coeff = -3. * (gsl_vector_get(f, 0) * gsl_vector_get(f_b, 1) - gsl_vector_get(f, 1) * gsl_vector_get(f_b, 0)) / (gsl_vector_get(f_a, 0) * gsl_vector_get(f_b, 1) - gsl_vector_get(f_a, 1) * gsl_vector_get(f_b, 0));
			long double x_b_coeff = -3. * (gsl_vector_get(f, 0) * gsl_vector_get(f_a, 1) - gsl_vector_get(f, 1) * gsl_vector_get(f_a, 0)) / (gsl_vector_get(f_b, 0) * gsl_vector_get(f_a, 1) - gsl_vector_get(f_b, 1) * gsl_vector_get(f_a, 0));
			gsl_vector_sub(x_a, x);
			gsl_vector_scale(x_a, x_a_coeff);
			gsl_vector_add(x_a, x);
			gsl_vector_sub(x_b, x);
			gsl_vector_scale(x_b, x_b_coeff);
			gsl_vector_add(x_b, x);
		}
		return GSL_FAILURE;
	}
	int MultiFunctionSolver::TriangleRotationIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<TriangleState *>(vstate);
		gsl_vector *x_trial = state->x_trial, *f_trial = state->f_trial;
		// find the middle point in edge a-b
		gsl_vector_memcpy(x_trial, state->x_a);
		gsl_vector_add(x_trial, state->x_b);
		gsl_vector_scale(x_trial, 0.5);
		for (int status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_scale(x_trial, min(0.99999999l, static_cast<long double>(gsl_vector_get(f_trial, 0))));
		}
		gsl_vector_memcpy(dx, x_trial);
		gsl_vector_sub(dx, x);
		if (PointInTriangle(f_trial, f, state->f_a)) {
			gsl_vector_memcpy(state->x_b, state->x_a);
			gsl_vector_memcpy(state->f_b, state->f_a);
			gsl_vector_memcpy(state->x_a, x);
			gsl_vector_memcpy(state->f_a, f);
			gsl_vector_memcpy(x, x_trial);
			gsl_vector_memcpy(f, f_trial);
			return GSL_SUCCESS;
		}
		if (PointInTriangle(f_trial, state->f_b, f)) {
			gsl_vector_memcpy(state->x_a, x);
			gsl_vector_memcpy(state->f_a, f);
			gsl_vector_memcpy(x, state->x_b);
			gsl_vector_memcpy(f, state->f_b);
			gsl_vector_memcpy(state->x_b, x_trial);
			gsl_vector_memcpy(state->f_b, f_trial);
			return GSL_SUCCESS;
		}
		return GSL_FAILURE;
	}
	int MultiFunctionSolver::TriangleLongestIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<TriangleState *>(vstate);
		gsl_vector *x_trial = state->x_trial, *f_trial = state->f_trial;
		// find the longest edge
		gsl_vector_memcpy(x_trial, x);
		gsl_vector_sub(x_trial, state->x_a);
		const long double x_a_norm = gsl_blas_dnrm2(x_trial);
		gsl_vector_memcpy(x_trial, x);
		gsl_vector_sub(x_trial, state->x_b);
		const long double x_b_norm = gsl_blas_dnrm2(x_trial);
		gsl_vector_memcpy(x_trial, state->x_a);
		gsl_vector_sub(x_trial, state->x_b);
		const long double a_b_norm = gsl_blas_dnrm2(x_trial);
		if (x_a_norm >= x_b_norm && x_a_norm >= a_b_norm) {
			gsl_blas_dswap(x, state->x_b);
			gsl_blas_dswap(f, state->f_b);
		} else if (x_b_norm >= x_a_norm && x_b_norm >= a_b_norm) {
			gsl_blas_dswap(x, state->x_a);
			gsl_blas_dswap(f, state->f_a);
		}
		gsl_vector_memcpy(x_trial, state->x_a);
		gsl_vector_add(x_trial, state->x_b);
		gsl_vector_scale(x_trial, 0.5);
		for (int status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_scale(x_trial, min(0.99999999l, static_cast<long double>(gsl_vector_get(f_trial, 0))));
		}
		gsl_vector_memcpy(dx, x_trial);
		gsl_vector_sub(dx, x);
		if (PointInTriangle(f_trial, f, state->f_a)) {
			gsl_vector_memcpy(state->x_b, state->x_a);
			gsl_vector_memcpy(state->f_b, state->f_a);
			gsl_vector_memcpy(state->x_a, x);
			gsl_vector_memcpy(state->f_a, f);
			gsl_vector_memcpy(x, x_trial);
			gsl_vector_memcpy(f, f_trial);
			return GSL_SUCCESS;
		}
		if (PointInTriangle(f_trial, state->f_b, f)) {
			gsl_vector_memcpy(state->x_a, x);
			gsl_vector_memcpy(state->f_a, f);
			gsl_vector_memcpy(x, state->x_b);
			gsl_vector_memcpy(f, state->f_b);
			gsl_vector_memcpy(state->x_b, x_trial);
			gsl_vector_memcpy(state->f_b, f_trial);
			return GSL_SUCCESS;
		}
		return GSL_FAILURE;
	}
	void MultiFunctionSolver::TriangleFree(void *vstate) {
		auto state = static_cast<TriangleState *>(vstate);
		gsl_vector_free(state->x_trial);
		gsl_vector_free(state->f_trial);
		gsl_vector_free(state->x_a);
		gsl_vector_free(state->x_b);
		gsl_vector_free(state->f_a);
		gsl_vector_free(state->f_b);
	}
	int MultiFunctionSolver::DirectionAlloc(void *vstate, size_t n) {
		auto state = static_cast<DirectionState *>(vstate);
		if (state->x_trial = gsl_vector_calloc(n); state->x_trial == nullptr)
			GSL_ERROR("failed to allocate space for x_trial", GSL_ENOMEM);
		if (state->x_rec = gsl_vector_calloc(n); state->x_rec == nullptr) {
			TriangleFree(vstate);
			GSL_ERROR("failed to allocate space for x_rec", GSL_ENOMEM);
		}
		if (state->f_trial = gsl_vector_calloc(n); state->f_trial == nullptr) {
			TriangleFree(vstate);
			GSL_ERROR("failed to allocate space for f_trial", GSL_ENOMEM);
		}
		if (state->f_rec = gsl_vector_calloc(n); state->f_rec == nullptr) {
			TriangleFree(vstate);
			GSL_ERROR("failed to allocate space for f_rec", GSL_ENOMEM);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::DirectionSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<DirectionState *>(vstate);
		state->directional_num = 1;
		state->delta_angle = 2. * M_PI_3;
		state->central_angle = 0.;
		if (int status = ScaleX(function, x, f); status != GSL_SUCCESS)
			return status;
		state->trust_radius = max(1., gsl_blas_dnrm2(x) / gsl_blas_dnrm2(f)) * GSL_SQRT_DBL_EPSILON;
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::DirectionIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		auto state = static_cast<DirectionState *>(vstate);
		gsl_vector *x_trial = state->x_trial, *x_rec = state->x_rec, *f_trial = state->f_trial, *f_rec = state->f_rec;
		int status;
		const long double f_norm = gsl_blas_dnrm2(f), trust_radius_limit = gsl_blas_dnrm2(x) * GSL_DBL_EPSILON * 100;
		for (long double effective_radius = max(trust_radius_limit, state->trust_radius * f_norm); state->directional_num < 64; effective_radius = max(trust_radius_limit, state->trust_radius * f_norm)) {
			long double min_norm, delta_angle = state->delta_angle, next_angle = state->central_angle;
			long double delta_x = effective_radius * cos(state->central_angle), delta_y = effective_radius * sin(state->central_angle);
			gsl_vector_set(x_trial, 0, gsl_vector_get(x, 0) + delta_x);
			gsl_vector_set(x_trial, 1, gsl_vector_get(x, 1) + delta_y);
			if (status = ScaleX(function, x_trial, f_trial); status != GSL_SUCCESS)
				return status;
			if (min_norm = gsl_blas_dnrm2(f_trial); min_norm < f_norm) {
				next_angle = state->central_angle;
				gsl_vector_memcpy(x_rec, x_trial);
				gsl_vector_memcpy(f_rec, f_trial);
			}
			for (int i = state->directional_num; i > 0; --i) {
				for (int sign = 1; sign > -2; sign -= 2) {
					long double angle = state->central_angle + sign * delta_angle;
					delta_x = state->trust_radius * cos(angle);
					delta_y = state->trust_radius * sin(angle);
					gsl_vector_set(x_trial, 0, gsl_vector_get(x, 0) + delta_x);
					gsl_vector_set(x_trial, 1, gsl_vector_get(x, 1) + delta_y);
					if (status = ScaleX(function, x_trial, f_trial); status != GSL_SUCCESS)
						return status;
					if (long double f_dir_norm = gsl_blas_dnrm2(f_trial); f_dir_norm < min_norm) {
						min_norm = f_dir_norm;
						next_angle = angle;
						gsl_vector_memcpy(x_rec, x_trial);
						gsl_vector_memcpy(f_rec, f_trial);
					}
				}
				if (delta_angle *= 2.; delta_angle > M_PI) {
					delta_angle *= 0.5;
					break;
				}
			}
			if (min_norm >= f_norm) {
				state->trust_radius *= 0.7;
				if (next_angle != state->central_angle + M_PI)
					state->delta_angle = max(state->delta_angle, abs(next_angle - state->central_angle)) * 0.5;
				else
					state->delta_angle = M_PI_4 * pow(0.5, state->directional_num);
				state->directional_num += 2;
				state->central_angle = ModBy2Pi(next_angle);
				continue;
			}
			delta_x = 0.5 * effective_radius * cos(next_angle);
			delta_y = 0.5 * effective_radius * sin(next_angle);
			while (effective_radius > trust_radius_limit) {
				gsl_vector_set(x_trial, 0, gsl_vector_get(x, 0) + delta_x);
				gsl_vector_set(x_trial, 1, gsl_vector_get(x, 1) + delta_y);
				if (status = ScaleX(function, x_trial, f_trial); status != GSL_SUCCESS)
					return status;
				if (long double f_mid_norm = gsl_blas_dnrm2(f_trial); f_mid_norm < min_norm) {
					min_norm = f_mid_norm;
					gsl_vector_memcpy(x_rec, x_trial);
					gsl_vector_memcpy(f_rec, f_trial);
					delta_x *= 0.5;
					delta_y *= 0.5;
					state->trust_radius *= 0.5;
				} else
					break;
			}
			gsl_vector_memcpy(dx, x_rec);
			gsl_vector_sub(dx, x);
			gsl_vector_memcpy(x, x_rec);
			gsl_vector_memcpy(f, f_rec);
			if (state->directional_num > 1) {
				if (next_angle != state->central_angle + M_PI)
					state->delta_angle = max(state->delta_angle, abs(next_angle - state->central_angle));
				else
					state->delta_angle = M_PI_2;
				state->directional_num = 1;
			} else {
				if (next_angle == state->central_angle)
					state->delta_angle *= 0.5;
				state->trust_radius *= 1.5;
			}
			state->central_angle = ModBy2Pi(next_angle);
			return GSL_SUCCESS;
		}
		return GSL_FAILURE;
	}
	void MultiFunctionSolver::DirectionFree(void *vstate) {
		auto state = static_cast<DirectionState *>(vstate);
		gsl_vector_free(state->x_trial);
		gsl_vector_free(state->x_rec);
		gsl_vector_free(state->f_trial);
		gsl_vector_free(state->f_rec);
	}

	static const gsl_multiroot_fsolver_type HybridType{"sbody_hybrid",
													   sizeof(HybridState),
													   &MultiFunctionSolver::HybridAlloc,
													   &MultiFunctionSolver::HybridSet,
													   &MultiFunctionSolver::HybridIterate,
													   &MultiFunctionSolver::HybridFree};

	static const gsl_multiroot_fsolver_type HybridAdditionType{"sbody_hybrid_addition",
															   sizeof(HybridAdditionState),
															   &MultiFunctionSolver::HybridAdditionAlloc,
															   &MultiFunctionSolver::HybridAdditionSet,
															   &MultiFunctionSolver::HybridAdditionIterate,
															   &MultiFunctionSolver::HybridAdditionFree};

	static const gsl_multiroot_fsolver_type DNewtonType{"sbody_dnewton",
														sizeof(DNewtonState),
														&MultiFunctionSolver::DNewtonAlloc,
														&MultiFunctionSolver::DNewtonSet,
														&MultiFunctionSolver::DNewtonIterate,
														&MultiFunctionSolver::DNewtonFree};

	static const gsl_multiroot_fsolver_type DNewtonRotationType{"sbody_dnewton_rotation",
																sizeof(DNewtonRotationTranslationState),
																&MultiFunctionSolver::DNewtonRotationTranslationAlloc,
																&MultiFunctionSolver::DNewtonRotationTranslationSet,
																&MultiFunctionSolver::DNewtonRotationIterate,
																&MultiFunctionSolver::DNewtonRotationTranslationFree};

	static const gsl_multiroot_fsolver_type DNewtonTranslationType{"sbody_dnewton_translation",
																   sizeof(DNewtonRotationTranslationState),
																   &MultiFunctionSolver::DNewtonRotationTranslationAlloc,
																   &MultiFunctionSolver::DNewtonRotationTranslationSet,
																   &MultiFunctionSolver::DNewtonTranslationIterate,
																   &MultiFunctionSolver::DNewtonRotationTranslationFree};

	static const gsl_multiroot_fsolver_type D2NewtonType{"sbody_d2newton",
														 sizeof(D2NewtonState),
														 &MultiFunctionSolver::D2NewtonAlloc,
														 &MultiFunctionSolver::D2NewtonSet,
														 &MultiFunctionSolver::D2NewtonIterate,
														 &MultiFunctionSolver::D2NewtonFree};

	static const gsl_multiroot_fsolver_type GradientType{"sbody_gradient",
														 sizeof(DNewtonState),
														 &MultiFunctionSolver::DNewtonAlloc,
														 &MultiFunctionSolver::DNewtonSet,
														 &MultiFunctionSolver::GradientIterate,
														 &MultiFunctionSolver::DNewtonFree};

	static const gsl_multiroot_fsolver_type ConjugateGradientType{"sbody_conjugate",
																  sizeof(ConjugateGradientState),
																  &MultiFunctionSolver::ConjugateGradientAlloc,
																  &MultiFunctionSolver::ConjugateGradientSet,
																  &MultiFunctionSolver::ConjugateGradientIterate,
																  &MultiFunctionSolver::ConjugateGradientFree};

	static const gsl_multiroot_fsolver_type TriangleType{"sbody_triangle_rotation",
														 sizeof(TriangleState),
														 &MultiFunctionSolver::TriangleAlloc,
														 &MultiFunctionSolver::TriangleSet,
														 &MultiFunctionSolver::TriangleRotationIterate,
														 &MultiFunctionSolver::TriangleFree};

	static const gsl_multiroot_fsolver_type DirectionType{"sbody_direction",
														  sizeof(DirectionState),
														  &MultiFunctionSolver::DirectionAlloc,
														  &MultiFunctionSolver::DirectionSet,
														  &MultiFunctionSolver::DirectionIterate,
														  &MultiFunctionSolver::DirectionFree};

	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_hybrid = &HybridType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_hybrid_addition = &HybridAdditionType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_dnewton = &DNewtonType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_dnewton_rotation = &DNewtonRotationType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_dnewton_translation = &DNewtonTranslationType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_d2newton = &D2NewtonType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_gradient = &GradientType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_conjugate = &ConjugateGradientType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_triangle = &TriangleType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_direction = &DirectionType;

	// MultiDerivativeSolver
	MultiDerivativeSolver::MultiDerivativeSolver(size_t n, const gsl_multiroot_fdfsolver_type *type) {
		solver_ = gsl_multiroot_fdfsolver_alloc(type, n);
	}
	MultiDerivativeSolver::MultiDerivativeSolver(gsl_multiroot_function_fdf *function, const gsl_vector *x, size_t n, const gsl_multiroot_fdfsolver_type *type) {
		solver_ = gsl_multiroot_fdfsolver_alloc(type, n);
		gsl_multiroot_fdfsolver_set(solver_, function, x);
	}
	MultiDerivativeSolver::~MultiDerivativeSolver() {
		gsl_multiroot_fdfsolver_free(solver_);
	}
	int MultiDerivativeSolver::Set(gsl_multiroot_function_fdf *function, const gsl_vector *x) {
		return gsl_multiroot_fdfsolver_set(solver_, function, x);
	}
	int MultiDerivativeSolver::Iterate() {
		return gsl_multiroot_fdfsolver_iterate(solver_);
	}
	int MultiDerivativeSolver::Solve(long double epsabs, long double epsrel, int max_iteration) {
		for (; max_iteration > 0 && gsl_multiroot_test_delta(gsl_multiroot_fdfsolver_dx(solver_), gsl_multiroot_fdfsolver_root(solver_), epsabs, epsrel) != GSL_SUCCESS; --max_iteration)
			if (int status = gsl_multiroot_fdfsolver_iterate(solver_); status != GSL_SUCCESS)
				return status;
		if (max_iteration <= 0)
			return GSL_EMAXITER;
		return GSL_SUCCESS;
	}
	gsl_vector *MultiDerivativeSolver::Root() {
		return gsl_multiroot_fdfsolver_root(solver_);
	}
	gsl_vector *MultiDerivativeSolver::Value() {
		return gsl_multiroot_fdfsolver_f(solver_);
	}
	gsl_vector *MultiDerivativeSolver::StepSize() {
		return gsl_multiroot_fdfsolver_dx(solver_);
	}

	MultiFunctionMinimizer::MultiFunctionMinimizer(size_t n, const gsl_multimin_fminimizer_type *type) {
		solver_ = gsl_multimin_fminimizer_alloc(type, n);
	}
	MultiFunctionMinimizer::~MultiFunctionMinimizer() {
		gsl_multimin_fminimizer_free(solver_);
	}
	int MultiFunctionMinimizer::Set(gsl_multimin_function *function, const gsl_vector *x, const gsl_vector *step_size) {
		solver_->fval = 100.;
		return gsl_multimin_fminimizer_set(solver_, function, x, step_size);
	}
	int MultiFunctionMinimizer::Iterate() {
		return gsl_multimin_fminimizer_iterate(solver_);
	}
	int MultiFunctionMinimizer::Solve(long double epsabs) {
		while (gsl_multimin_test_size(gsl_multimin_fminimizer_size(solver_), epsabs) != GSL_SUCCESS || gsl_multimin_test_size(gsl_multimin_fminimizer_minimum(solver_), epsabs) != GSL_SUCCESS)
			if (int status = gsl_multimin_fminimizer_iterate(solver_); status != GSL_SUCCESS)
				return status;
		return GSL_SUCCESS;
	}

	gsl_vector *MultiFunctionMinimizer::Root() {
		return gsl_multimin_fminimizer_x(solver_);
	}
	long double MultiFunctionMinimizer::Value() {
		return gsl_multimin_fminimizer_minimum(solver_);
	}
	long double MultiFunctionMinimizer::StepSize() {
		return gsl_multimin_fminimizer_size(solver_);
	}
	// GslBlock
	GslBlock::~GslBlock() {
		for (auto vector : vectors_)
			gsl_vector_free(vector);
		for (auto matrix : matrices_)
			gsl_matrix_free(matrix);
		for (auto permutation : permutations_)
			gsl_permutation_free(permutation);
	}
	gsl_vector *GslBlock::VectorAlloc(size_t n) {
		vectors_.push_back(gsl_vector_alloc(n));
		return vectors_.back();
	}
	gsl_vector *GslBlock::VectorCalloc(size_t n) {
		vectors_.push_back(gsl_vector_calloc(n));
		return vectors_.back();
	}
	gsl_vector *GslBlock::VectorAllocFromBlock(gsl_block *block, const size_t offset, const size_t n, const size_t stride) {
		vectors_.push_back(gsl_vector_alloc_from_block(block, offset, n, stride));
		return vectors_.back();
	}
	gsl_vector *GslBlock::VectorAllocRowFromMatrix(gsl_matrix *matrix, const size_t i) {
		vectors_.push_back(gsl_vector_alloc_row_from_matrix(matrix, i));
		return vectors_.back();
	}
	gsl_vector *GslBlock::VectorAllocColFromMatrix(gsl_matrix *matrix, const size_t j) {
		vectors_.push_back(gsl_vector_alloc_col_from_matrix(matrix, j));
		return vectors_.back();
	}
	gsl_matrix *GslBlock::MatrixAlloc(size_t n1, size_t n2) {
		matrices_.push_back(gsl_matrix_alloc(n1, n2));
		return matrices_.back();
	}
	gsl_matrix *GslBlock::MatrixCalloc(size_t n1, size_t n2) {
		matrices_.push_back(gsl_matrix_calloc(n1, n2));
		return matrices_.back();
	}
	gsl_permutation *GslBlock::PermutationAlloc(size_t n) {
		permutations_.push_back(gsl_permutation_alloc(n));
		return permutations_.back();
	}
	gsl_permutation *GslBlock::PermutationCalloc(size_t n) {
		permutations_.push_back(gsl_permutation_calloc(n));
		return permutations_.back();
	}

	int CoordinateOrthogonalization(const gsl_vector *x, gsl_matrix *coordinate) {
		size_t n = x->size;
		if (n != coordinate->size1 || n != coordinate->size2) {
			PrintlnError("CoordinateOrthogonalization: Invalid size");
			return GSL_EBADLEN;
		}
		gsl_matrix_set_zero(coordinate);
		gsl_vector_view column_view[n];
		for (int i = 0; i < n; ++i)
			column_view[i] = gsl_matrix_column(coordinate, i);
		int coordinate_idx = 0;
		gsl_blas_daxpy(1. / gsl_blas_dnrm2(x), x, &column_view[0].vector);
		for (int i = 1; i < n; ++i) {
			gsl_vector_set_basis(&column_view[i].vector, coordinate_idx++);
			for (int j = 0; j < i; ++j) {
				double dot_product;
				gsl_blas_ddot(&column_view[i].vector, &column_view[j].vector, &dot_product);
				gsl_blas_daxpy(-dot_product, &column_view[j].vector, &column_view[i].vector);
			}
			long double base_norm = gsl_blas_dnrm2(&column_view[i].vector);
			if (base_norm > GSL_SQRT_DBL_EPSILON)
				gsl_vector_scale(&column_view[i].vector, 1. / gsl_blas_dnrm2(&column_view[i].vector));
			else if (coordinate_idx == n)
				return GSL_EFAILED;
			else
				--i;
		}
		return GSL_SUCCESS;
	}

	bool PointInTriangle(const gsl_vector *a, const gsl_vector *b, const gsl_vector *c, const gsl_vector *p) {
		if (p == nullptr)
			return PointInTriangle(gsl_vector_get(a, 0), gsl_vector_get(a, 1), gsl_vector_get(b, 0), gsl_vector_get(b, 1), gsl_vector_get(c, 0), gsl_vector_get(c, 1));
		return PointInTriangle(gsl_vector_get(a, 0), gsl_vector_get(a, 1), gsl_vector_get(b, 0), gsl_vector_get(b, 1), gsl_vector_get(c, 0), gsl_vector_get(c, 1), gsl_vector_get(p, 0), gsl_vector_get(p, 1));
	}
} // namespace SBody
