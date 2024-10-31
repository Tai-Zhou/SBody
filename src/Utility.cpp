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

#include "Utility.h"

#include <cmath>
#include <vector>

#include <fmt/core.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_ellint.h>

#include "IO.h"

using namespace std;

namespace SBody {
	double absolute_accuracy = 1e-15, relative_accuracy = 1e-15;

	// Integrator
	Integrator::Integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), void *params, const gsl_odeiv2_step_type *type) : control_(gsl_odeiv2_control_y_new(absolute_accuracy, relative_accuracy)), evolve_(gsl_odeiv2_evolve_alloc(8)), step_(gsl_odeiv2_step_alloc(type, 8)) {
		system_ = gsl_odeiv2_system{function, jacobian, 8UL, params};
	}
	Integrator::~Integrator() {
		gsl_odeiv2_control_free(control_);
		gsl_odeiv2_evolve_free(evolve_);
		gsl_odeiv2_step_free(step_);
	}
	int Integrator::Apply(double *t, double t1, double *h, double *y) {
		int status = 0;
		double theta_0 = y[2];
		if (*h > 0)
			while (status <= 0 && *t < t1)
				status = gsl_odeiv2_evolve_apply(evolve_, control_, step_, &system_, t, t1, h, y);
		else
			while (status <= 0 && *t > t1)
				status = gsl_odeiv2_evolve_apply(evolve_, control_, step_, &system_, t, t1, h, y);
		MapTheta(theta_0, y);
		y[3] = ModBy2Pi(y[3]);
		return status;
	}
	int Integrator::ApplyStep(double *t, double t1, double *h, double *y) {
		double theta_0 = y[2];
		int status = gsl_odeiv2_evolve_apply(evolve_, control_, step_, &system_, t, t1, h, y);
		MapTheta(theta_0, y);
		y[3] = ModBy2Pi(y[3]);
		return status;
	}
	int Integrator::ApplyFixedStep(double *t, const double h, double *y) {
		double theta_0 = y[2];
		int status = gsl_odeiv2_evolve_apply_fixed_step(evolve_, control_, step_, &system_, t, h, y);
		MapTheta(theta_0, y);
		y[3] = ModBy2Pi(y[3]);
		return status;
	}
	int Integrator::Reset() {
		if (int status = gsl_odeiv2_evolve_reset(evolve_); status != GSL_SUCCESS)
			return status;
		if (int status = gsl_odeiv2_step_reset(step_); status != GSL_SUCCESS)
			return status;
		return GSL_SUCCESS;
	}

	// FunctionSolver
	FunctionSolver::FunctionSolver(const gsl_root_fsolver_type *type) {
		solver_ = gsl_root_fsolver_alloc(type);
	}
	FunctionSolver::~FunctionSolver() {
		gsl_root_fsolver_free(solver_);
	}
	int FunctionSolver::Set(gsl_function *function, double lower, double upper) {
		return gsl_root_fsolver_set(solver_, function, lower, upper);
	}
	int FunctionSolver::Iterate() {
		return gsl_root_fsolver_iterate(solver_);
	}
	int FunctionSolver::Solve() {
		while (gsl_root_test_interval(gsl_root_fsolver_x_lower(solver_), gsl_root_fsolver_x_upper(solver_), absolute_accuracy, relative_accuracy) != GSL_SUCCESS)
			if (int status = gsl_root_fsolver_iterate(solver_); status != GSL_SUCCESS)
				return status;
		return GSL_SUCCESS;
	}
	double FunctionSolver::Root() {
		return gsl_root_fsolver_root(solver_);
	}
	double FunctionSolver::Lower() {
		return gsl_root_fsolver_x_lower(solver_);
	}
	double FunctionSolver::Upper() {
		return gsl_root_fsolver_x_upper(solver_);
	}

	// DerivativeSolver
	DerivativeSolver::DerivativeSolver(const gsl_root_fdfsolver_type *type) {
		solver_ = gsl_root_fdfsolver_alloc(type);
	}
	DerivativeSolver::DerivativeSolver(gsl_function_fdf *function, double root, const gsl_root_fdfsolver_type *type) {
		solver_ = gsl_root_fdfsolver_alloc(type);
		gsl_root_fdfsolver_set(solver_, function, root);
	}
	DerivativeSolver::~DerivativeSolver() {
		gsl_root_fdfsolver_free(solver_);
	}
	int DerivativeSolver::Set(gsl_function_fdf *function, double root) {
		return gsl_root_fdfsolver_set(solver_, function, root);
	}
	int DerivativeSolver::Iterate() {
		return gsl_root_fdfsolver_iterate(solver_);
	}
	int DerivativeSolver::Solve(double epsabs, int max_iteration) {
		for (; max_iteration > 0 && gsl_root_test_residual(GSL_FN_FDF_EVAL_F(solver_->fdf, solver_->root), absolute_accuracy) != GSL_SUCCESS; --max_iteration)
			if (int status = gsl_root_fdfsolver_iterate(solver_); status != GSL_SUCCESS)
				return status;
		if (max_iteration <= 0)
			return GSL_EMAXITER;
		return GSL_SUCCESS;
	}
	double DerivativeSolver::Root() {
		return gsl_root_fdfsolver_root(solver_);
	}

	// MultiFunctionSolver
	MultiFunctionSolver::MultiFunctionSolver(size_t n, const gsl_multiroot_fsolver_type *type) {
		solver_ = gsl_multiroot_fsolver_alloc(type, n);
	}
	MultiFunctionSolver::MultiFunctionSolver(gsl_multiroot_function *function, const gsl_vector *x, const gsl_multiroot_fsolver_type *type) {
		solver_ = gsl_multiroot_fsolver_alloc(type, function->n);
		gsl_multiroot_fsolver_set(solver_, function, x);
	}
	MultiFunctionSolver::~MultiFunctionSolver() {
		gsl_multiroot_fsolver_free(solver_);
	}
	int MultiFunctionSolver::Set(gsl_multiroot_function *function, const gsl_vector *x) {
		return gsl_multiroot_fsolver_set(solver_, function, x);
	}
	int MultiFunctionSolver::Iterate() {
		return gsl_multiroot_fsolver_iterate(solver_);
	}
	int MultiFunctionSolver::Solve(double epsabs, int max_iteration) {
		for (; max_iteration > 0 && gsl_multiroot_test_residual(gsl_multiroot_fsolver_f(solver_), epsabs) != GSL_SUCCESS; --max_iteration)
			if (int status = gsl_multiroot_fsolver_iterate(solver_); status != GSL_SUCCESS)
				return status;
		if (max_iteration <= 0)
			return GSL_EMAXITER;
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::Solve(double epsabs, double epsrel, int max_iteration) {
		for (; max_iteration > 0 && (gsl_multiroot_test_delta(gsl_multiroot_fsolver_dx(solver_), gsl_multiroot_fsolver_root(solver_), epsabs, epsrel) != GSL_SUCCESS || gsl_multiroot_test_residual(gsl_multiroot_fsolver_f(solver_), epsabs) != GSL_SUCCESS); --max_iteration)
			if (int status = gsl_multiroot_fsolver_iterate(solver_); status != GSL_SUCCESS)
				return status;
		if (max_iteration <= 0) {
			// for (int status = gsl_multiroot_fsolver_iterate(solver_); status == GSL_SUCCESS;)
			// 	status = gsl_multiroot_fsolver_iterate(solver_);
			return GSL_EMAXITER;
		}
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
	int MultiFunctionSolver::SetIterationCoefficient(double coefficient) {
		if (solver_->type == gsl_multiroot_fsolver_sbody_dnewton)
			static_cast<DnewtonState *>(solver_->state)->iteration_coefficient = coefficient;
		else if (solver_->type == gsl_multiroot_fsolver_sbody_hybrid || solver_->type == gsl_multiroot_fsolver_sbody_hybrids)
			static_cast<HybridState *>(solver_->state)->iteration_coefficient = coefficient;
		else
			return GSL_FAILURE;
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::Jacobian(gsl_multiroot_function *F, const gsl_vector *x, const gsl_vector *f, double epsrel, gsl_matrix *jacobian) {
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
			double x_j = gsl_vector_get(x, j);
			double dx = abs(x_j) >= 1 ? -epsrel * x_j : -epsrel * GSL_SIGN(x_j);
			gsl_vector_set(x1, j, x_j + dx);
			gsl_vector_view col_j = gsl_matrix_column(jacobian, j);
			if (int f_stat = GSL_MULTIROOT_FN_EVAL(F, x1, &col_j.vector); f_stat != GSL_SUCCESS) {
				gsl_vector_free(x1);
				return GSL_EBADFUNC;
			}
			gsl_blas_daxpy(-1., f, &col_j.vector);
			gsl_blas_dscal(1. / dx, &col_j.vector);
			if (gsl_vector_isnull(&col_j.vector)) {
				gsl_vector_free(x1);
				return GSL_ESING;
			}
			gsl_vector_set(x1, j, x_j);
		}
		gsl_vector_free(x1);
		return GSL_SUCCESS;
	}
	void MultiFunctionSolver::ComputeDiag(const gsl_matrix *J, gsl_vector *diag) {
		size_t j, n = diag->size;
		for (j = 0; j < n; j++) {
			auto column = gsl_matrix_const_column(J, j);
			double norm = gsl_blas_dnrm2(&column.vector);
			if (norm == 0.)
				norm = 1.0;
			gsl_vector_set(diag, j, norm);
		}
	}
	void MultiFunctionSolver::UpdateDiag(const gsl_matrix *J, gsl_vector *diag) {
		size_t j, n = diag->size;
		for (j = 0; j < n; j++) {
			auto column = gsl_matrix_const_column(J, j);
			double norm = gsl_blas_dnrm2(&column.vector);
			if (norm == 0.)
				norm = 1.0;
			if (norm > gsl_vector_get(diag, j))
				gsl_vector_set(diag, j, norm);
		}
	}
	double MultiFunctionSolver::ScaledEnorm(const gsl_vector *d, const gsl_vector *f) {
		double e2 = 0.;
		size_t i, n = f->size;
		for (i = 0; i < n; i++)
			e2 += gsl_pow_2(gsl_vector_get(f, i) * gsl_vector_get(d, i));
		return sqrt(e2);
	}
	int MultiFunctionSolver::Dogleg(const gsl_matrix *r, const gsl_vector *qtf, const gsl_vector *diag, double delta, gsl_vector *newton, gsl_vector *gradient, gsl_vector *p) {
		double qnorm, gnorm, sgnorm, bnorm, temp;
		gsl_linalg_R_solve(r, qtf, p);
		gsl_vector_scale(p, -1.);
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
			double bg = bnorm / gnorm;
			double bq = bnorm / qnorm;
			double dq = delta / qnorm;
			double dq2 = dq * dq;
			double sd = sgnorm / delta;
			double sd2 = sd * sd;

			double t1 = bg * bq * sd;
			double u = t1 - dq;
			double t2 = t1 - dq * sd2 + sqrt(u * u + (1 - dq2) * (1 - sd2));

			double alpha = dq * (1 - sd2) / t2;
			double beta = (1 - alpha) * sgnorm;
			gsl_vector_set_zero(p);
			gsl_blas_daxpy(alpha, newton, p);
			gsl_blas_daxpy(beta, gradient, p);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::HybridAlloc(void *vstate, size_t n) {
		HybridState *state = static_cast<HybridState *>(vstate);
		gsl_matrix *J, *q, *r;
		gsl_vector *tau, *diag, *qtf, *newton, *gradient, *x_trial, *f_trial, *df, *qtdf, *rdx, *w, *v;
		if (J = gsl_matrix_calloc(n, n); J == nullptr)
			GSL_ERROR("failed to allocate space for J", GSL_ENOMEM);
		state->J = J;
		if (q = gsl_matrix_calloc(n, n); q == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for q", GSL_ENOMEM);
		}
		state->q = q;
		if (r = gsl_matrix_calloc(n, n); r == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for r", GSL_ENOMEM);
		}
		state->r = r;
		if (tau = gsl_vector_calloc(n); tau == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for tau", GSL_ENOMEM);
		}
		state->tau = tau;
		if (diag = gsl_vector_calloc(n); diag == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for diag", GSL_ENOMEM);
		}
		state->diag = diag;
		if (qtf = gsl_vector_calloc(n); qtf == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for qtf", GSL_ENOMEM);
		}
		state->qtf = qtf;
		if (newton = gsl_vector_calloc(n); newton == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for newton", GSL_ENOMEM);
		}
		state->newton = newton;
		if (gradient = gsl_vector_calloc(n); gradient == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for gradient", GSL_ENOMEM);
		}
		state->gradient = gradient;
		if (x_trial = gsl_vector_calloc(n); x_trial == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for x_trial", GSL_ENOMEM);
		}
		state->x_trial = x_trial;
		if (f_trial = gsl_vector_calloc(n); f_trial == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for f_trial", GSL_ENOMEM);
		}
		state->f_trial = f_trial;
		if (df = gsl_vector_calloc(n); df == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for df", GSL_ENOMEM);
		}
		state->df = df;
		if (qtdf = gsl_vector_calloc(n); qtdf == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for qtdf", GSL_ENOMEM);
		}
		state->qtdf = qtdf;
		if (rdx = gsl_vector_calloc(n); rdx == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for rdx", GSL_ENOMEM);
		}
		state->rdx = rdx;
		if (w = gsl_vector_calloc(n); w == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for w", GSL_ENOMEM);
		}
		state->w = w;
		if (v = gsl_vector_calloc(n); v == nullptr) {
			HybridFree(state);
			GSL_ERROR("failed to allocate space for v", GSL_ENOMEM);
		}
		state->v = v;
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::HybridSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		return HybridSetCore(vstate, function, x, f, dx, false);
	}
	int MultiFunctionSolver::HybridScaleSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		return HybridSetCore(vstate, function, x, f, dx, true);
	}
	int MultiFunctionSolver::HybridSetCore(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx, bool scale) {
		HybridState *state = static_cast<HybridState *>(vstate);
		gsl_matrix *J = state->J;
		gsl_matrix *q = state->q;
		gsl_matrix *r = state->r;
		gsl_vector *tau = state->tau;
		gsl_vector *diag = state->diag;
		if (int status = GSL_MULTIROOT_FN_EVAL(function, x, f); status == GSL_EDOM)
			gsl_vector_scale(x, gsl_vector_get(f, 0));
		else if (status != GSL_SUCCESS)
			return status;
		if (int status = Jacobian(function, x, f, GSL_SQRT_DBL_EPSILON, J); status != GSL_SUCCESS)
			return status;
		state->iter = 1;
		state->fnorm = gsl_blas_dnrm2(f);
		state->ncfail = 0;
		state->ncsuc = 0;
		state->nslow1 = 0;
		state->nslow2 = 0;
		gsl_vector_set_zero(dx);
		/* Store column norms in diag */
		if (scale)
			ComputeDiag(J, diag);
		else
			gsl_vector_set_all(diag, 1.0);
		/* Set delta to factor |D x| or to factor if |D x| is zero */
		double Dx = ScaledEnorm(diag, x);
		state->delta = Dx > 0. ? 100. * Dx : 100.;
		/* Factorize J into QR decomposition */
		if (int status = gsl_linalg_QR_decomp(J, tau); status != GSL_SUCCESS)
			return status;
		return gsl_linalg_QR_unpack(J, tau, q, r);
	}
	int MultiFunctionSolver::HybridIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		return HybridIterateCore(vstate, function, x, f, dx, false);
	}
	int MultiFunctionSolver::HybridScaleIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		return HybridIterateCore(vstate, function, x, f, dx, true);
	}
	int MultiFunctionSolver::HybridIterateCore(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx, bool scale) {
		HybridState *state = static_cast<HybridState *>(vstate);

		const double fnorm = state->fnorm;

		gsl_matrix *J = state->J;
		gsl_matrix *q = state->q;
		gsl_matrix *r = state->r;
		gsl_vector *tau = state->tau;
		gsl_vector *diag = state->diag;
		gsl_vector *qtf = state->qtf;
		gsl_vector *x_trial = state->x_trial;
		gsl_vector *f_trial = state->f_trial;
		gsl_vector *df = state->df;
		gsl_vector *qtdf = state->qtdf;
		gsl_vector *rdx = state->rdx;
		gsl_vector *w = state->w;
		gsl_vector *v = state->v;

		double prered, actred;
		double pnorm, fnorm1, fnorm1p;
		double ratio;
		double p1 = 0.1, p5 = 0.5, p001 = 0.001, p0001 = 0.0001;

		/* Compute qtf = Q^T f */

		gsl_blas_dgemv(CblasTrans, 1., q, f, 0., qtf);

		/* Compute dogleg step */

		Dogleg(r, qtf, diag, state->delta, state->newton, state->gradient, dx);

		/* Take a trial step */

		gsl_vector_memcpy(state->x_trial, x);
		gsl_vector_add(state->x_trial, dx);

		pnorm = ScaledEnorm(diag, dx);

		if (state->iter == 1) {
			if (pnorm < state->delta) {
				state->delta = pnorm;
			}
		}

		/* Evaluate function at x + p */

		if (int status = GSL_MULTIROOT_FN_EVAL(function, x_trial, f_trial); status == GSL_EDOM)
			gsl_vector_scale(x, gsl_vector_get(f, 0));
		else if (status != GSL_SUCCESS)
			return GSL_EBADFUNC;

		/* Set df = f_trial - f */
		gsl_vector_memcpy(df, f_trial);
		gsl_vector_sub(df, f);

		/* Compute the scaled actual reduction */

		fnorm1 = gsl_blas_dnrm2(f_trial);

		actred = fnorm1 < fnorm ? 1. - gsl_pow_2(fnorm1 / fnorm) : -1.;

		/* Compute rdx = R dx */

		gsl_blas_dgemv(CblasNoTrans, 1., r, dx, 0., rdx);

		/* Compute the scaled predicted reduction phi1p = |Q^T f + R dx| */

		gsl_vector_add(qtf, rdx);
		fnorm1p = gsl_blas_dnrm2(qtf);

		prered = fnorm1p < fnorm ? 1. - gsl_pow_2(fnorm1p / fnorm) : 0.;

		/* Compute the ratio of the actual to predicted reduction */

		if (prered > 0)
			ratio = actred / prered;
		else
			ratio = 0.;

		/* Update the step bound */

		if (ratio < p1) {
			state->ncsuc = 0;
			state->ncfail++;
			state->delta *= p5;
		} else {
			state->ncfail = 0;
			state->ncsuc++;

			if (ratio >= p5 || state->ncsuc > 1)
				state->delta = GSL_MAX(state->delta, pnorm / p5);
			if (fabs(ratio - 1) <= p1)
				state->delta = pnorm / p5;
		}

		/* Test for successful iteration */

		if (ratio >= p0001) {
			gsl_vector_memcpy(x, x_trial);
			gsl_vector_memcpy(f, f_trial);
			state->fnorm = fnorm1;
			state->iter++;
		}

		/* Determine the progress of the iteration */

		state->nslow1++;
		if (actred >= p001)
			state->nslow1 = 0;

		if (actred >= p1)
			state->nslow2 = 0;

		if (state->ncfail == 2) {
			Jacobian(function, x, f, GSL_SQRT_DBL_EPSILON, J);
			state->nslow2++;
			if (state->iter == 1) {
				if (scale)
					ComputeDiag(J, diag);
				double Dx = ScaledEnorm(diag, x);
				state->delta = Dx > 0. ? 100. * Dx : 100.;
			} else if (scale)
				UpdateDiag(J, diag);

			/* Factorize J into QR decomposition */

			gsl_linalg_QR_decomp(J, tau);
			gsl_linalg_QR_unpack(J, tau, q, r);

			return GSL_SUCCESS;
		}
		/* Compute qtdf = Q^T df, w = (Q^T df - R dx)/|dx|,  v = D^2 dx/|dx| */
		gsl_blas_dgemv(CblasTrans, 1., q, df, 0., qtdf);
		gsl_vector_memcpy(w, qtdf);
		gsl_vector_sub(w, rdx);
		gsl_vector_scale(w, 1. / pnorm);
		gsl_vector_memcpy(v, diag);
		gsl_vector_mul(v, diag);
		gsl_vector_mul(v, dx);
		gsl_vector_scale(v, 1. / pnorm);

		/* Rank-1 update of the jacobian Q'R' = Q(R + w v^T) */
		gsl_linalg_QR_update(q, r, w, v);
		/* No progress as measured by jacobian evaluations */
		if (state->nslow2 == 5)
			return GSL_ENOPROGJ;
		/* No progress as measured by function evaluations */
		if (state->nslow1 == 10)
			return GSL_ENOPROG;
		return GSL_SUCCESS;
	}
	void MultiFunctionSolver::HybridFree(void *vstate) {
		HybridState *state = static_cast<HybridState *>(vstate);

		gsl_vector_free(state->v);
		gsl_vector_free(state->w);
		gsl_vector_free(state->rdx);
		gsl_vector_free(state->qtdf);
		gsl_vector_free(state->df);
		gsl_vector_free(state->f_trial);
		gsl_vector_free(state->x_trial);
		gsl_vector_free(state->gradient);
		gsl_vector_free(state->newton);
		gsl_vector_free(state->qtf);
		gsl_vector_free(state->diag);
		gsl_vector_free(state->tau);
		gsl_matrix_free(state->r);
		gsl_matrix_free(state->q);
		gsl_matrix_free(state->J);
	}
	int MultiFunctionSolver::DnewtonAlloc(void *vstate, size_t n) {
		DnewtonState *state = static_cast<DnewtonState *>(vstate);
		if (state->dx = gsl_vector_calloc(n); state->dx == nullptr)
			GSL_ERROR("failed to allocate space for dx", GSL_ENOMEM);
		if (state->lu = gsl_matrix_calloc(n, n); state->lu == nullptr) {
			DnewtonFree(vstate);
			GSL_ERROR("failed to allocate space for lu", GSL_ENOMEM);
		}
		if (state->permutation = gsl_permutation_calloc(n); state->permutation == nullptr) {
			DnewtonFree(vstate);
			GSL_ERROR("failed to allocate space for permutation", GSL_ENOMEM);
		}
		if (state->J = gsl_matrix_calloc(n, n); state->J == nullptr) {
			DnewtonFree(vstate);
			GSL_ERROR("failed to allocate space for d", GSL_ENOMEM);
		}
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::DnewtonSet(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		DnewtonState *state = static_cast<DnewtonState *>(vstate);
		state->iteration_coefficient = 1.0;
		for (int status = GSL_MULTIROOT_FN_EVAL(function, x, f); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x, f)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_scale(x, gsl_vector_get(f, 0));
		}
		if (int status = Jacobian(function, x, f, GSL_SQRT_DBL_EPSILON, state->J); status != GSL_SUCCESS)
			return status;
		gsl_vector_set_zero(dx);
		return GSL_SUCCESS;
	}
	int MultiFunctionSolver::DnewtonIterate(void *vstate, gsl_multiroot_function *function, gsl_vector *x, gsl_vector *f, gsl_vector *dx) {
		DnewtonState *state = static_cast<DnewtonState *>(vstate);
		int signum;
		gsl_matrix_memcpy(state->lu, state->J);
		if (int status = gsl_linalg_LU_decomp(state->lu, state->permutation, &signum); status != GSL_SUCCESS)
			return status;
		gsl_vector_memcpy(state->dx, dx);
		if (int status = gsl_linalg_LU_solve(state->lu, state->permutation, f, dx); status != GSL_SUCCESS)
			return status;
		const double norm_dx = gsl_blas_dnrm2(dx) * state->iteration_coefficient;
		const double norm_x = gsl_blas_dnrm2(x);
		if (norm_dx > 0.5 * state->iteration_coefficient * norm_x + GSL_SQRT_DBL_EPSILON)
			gsl_blas_daxpy(-(0.5 * state->iteration_coefficient * norm_x + GSL_SQRT_DBL_EPSILON) / norm_dx, dx, x);
		else
			gsl_blas_daxpy(-state->iteration_coefficient, dx, x);
		for (int status = GSL_MULTIROOT_FN_EVAL(function, x, f); status != GSL_SUCCESS; status = GSL_MULTIROOT_FN_EVAL(function, x, f)) {
			if (status != GSL_EDOM)
				return GSL_EBADFUNC;
			gsl_vector_memcpy(dx, x);
			gsl_vector_scale(x, gsl_vector_get(f, 0));
			gsl_vector_sub(dx, x);
			state->iteration_coefficient *= min(0.9, gsl_vector_get(f, 0));
			if (state->iteration_coefficient < 0.1)
				state->iteration_coefficient = 0.1;
		}
		if (int status = Jacobian(function, x, f, GSL_SQRT_DBL_EPSILON, state->J); status != GSL_SUCCESS)
			return status;
		return GSL_SUCCESS;
	}
	void MultiFunctionSolver::DnewtonFree(void *vstate) {
		DnewtonState *state = static_cast<DnewtonState *>(vstate);
		gsl_vector_free(state->dx);
		gsl_matrix_free(state->J);
		gsl_matrix_free(state->lu);
		gsl_permutation_free(state->permutation);
	}

	static const gsl_multiroot_fsolver_type HybridType{"sbody_hybrid",
													   sizeof(HybridState),
													   &MultiFunctionSolver::HybridAlloc,
													   &MultiFunctionSolver::HybridSet,
													   &MultiFunctionSolver::HybridIterate,
													   &MultiFunctionSolver::HybridFree};

	static const gsl_multiroot_fsolver_type HybridsType{"sbody_hybrids",
														sizeof(HybridState),
														&MultiFunctionSolver::HybridAlloc,
														&MultiFunctionSolver::HybridScaleSet,
														&MultiFunctionSolver::HybridScaleIterate,
														&MultiFunctionSolver::HybridFree};

	static const gsl_multiroot_fsolver_type DnewtonType{"sbody_dnewton",
														sizeof(DnewtonState),
														&MultiFunctionSolver::DnewtonAlloc,
														&MultiFunctionSolver::DnewtonSet,
														&MultiFunctionSolver::DnewtonIterate,
														&MultiFunctionSolver::DnewtonFree};

	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_hybrid = &HybridType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_hybrids = &HybridsType;
	const gsl_multiroot_fsolver_type *gsl_multiroot_fsolver_sbody_dnewton = &DnewtonType;

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
	int MultiDerivativeSolver::Solve(double epsabs, double epsrel, int max_iteration) {
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
	int MultiFunctionMinimizer::Solve(double epsabs) {
		while (gsl_multimin_test_size(gsl_multimin_fminimizer_size(solver_), epsabs) != GSL_SUCCESS || gsl_multimin_test_size(gsl_multimin_fminimizer_minimum(solver_), epsabs) != GSL_SUCCESS)
			if (int status = gsl_multimin_fminimizer_iterate(solver_); status != GSL_SUCCESS)
				return status;
		return GSL_SUCCESS;
	}

	gsl_vector *MultiFunctionMinimizer::Root() {
		return gsl_multimin_fminimizer_x(solver_);
	}
	double MultiFunctionMinimizer::Value() {
		return gsl_multimin_fminimizer_minimum(solver_);
	}
	double MultiFunctionMinimizer::StepSize() {
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

	// Maths
	double SquareRoot(double x) {
		if (x <= 0.)
			return 0.;
		return sqrt(x);
	}
	double SignSquare(double x) {
		if (x < 0.)
			return -x * x;
		return x * x;
	}
	double SignSquareRoot(double x) {
		if (x < 0.)
			return -sqrt(-x);
		return sqrt(x);
	}
	double Dot(const double x[], const double y[], size_t dimension) {
		if (dimension == 3)
			return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
		return cblas_ddot(dimension, x, 1, y, 1);
	}
	double Dot(const double x[], size_t dimension) {
		if (dimension == 3)
			return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
		return cblas_ddot(dimension, x, 1, x, 1);
	}
	double Norm(const double x[], size_t dimension) {
		if (dimension == 3)
			return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
		return cblas_dnrm2(dimension, x, 1);
	}
	int Cross(const double x[], const double y[], double z[]) {
		z[0] = x[1] * y[2] - x[2] * y[1];
		z[1] = x[2] * y[0] - x[0] * y[2];
		z[2] = x[0] * y[1] - x[1] * y[0];
		return GSL_SUCCESS;
	}
	double DotCross(const double x[], const double y[], const double z[]) {
		return x[0] * (y[1] * z[2] - y[2] * z[1]) + x[1] * (y[2] * z[0] - y[0] * z[2]) + x[2] * (y[0] * z[1] - y[1] * z[0]);
	}
	double TriangleArea(double a, double b, double c) {
		const double s = 0.5 * (a + b + c);
		return sqrt(s * (s - a) * (s - b) * (s - c));
	}
	double TriangleArea(const double x[], const double y[], const double z[]) {
		const double a = sqrt(gsl_pow_2(y[0] - x[0]) + gsl_pow_2(y[1] - x[1]) + gsl_pow_2(y[2] - x[2])),
					 b = sqrt(gsl_pow_2(z[0] - y[0]) + gsl_pow_2(z[1] - y[1]) + gsl_pow_2(z[2] - y[2])),
					 c = sqrt(gsl_pow_2(x[0] - z[0]) + gsl_pow_2(x[1] - z[1]) + gsl_pow_2(x[2] - z[2])),
					 s = 0.5 * (a + b + c);
		return sqrt(s * (s - a) * (s - b) * (s - c));
	}
	int RotateAroundAxis(double x[], Axis axis, double angle) {
		const double sin_angle = sin(angle), cos_angle = cos(angle);
		switch (axis) {
		case X:
			angle = x[1] * cos_angle - x[2] * sin_angle;
			x[2] = x[1] * sin_angle + x[2] * cos_angle;
			x[1] = angle;
			return GSL_SUCCESS;
		case Y:
			angle = x[2] * cos_angle - x[0] * sin_angle;
			x[0] = x[2] * sin_angle + x[0] * cos_angle;
			x[2] = angle;
			return GSL_SUCCESS;
		case Z:
			angle = x[0] * cos_angle - x[1] * sin_angle;
			x[1] = x[0] * sin_angle + x[1] * cos_angle;
			x[0] = angle;
			return GSL_SUCCESS;
		default:
			return GSL_EINVAL;
		}
	}
	int CartesianToSpherical(const double cartesian[], double spherical[], size_t dimension) {
#ifndef GSL_RANGE_CHECK_OFF
		if (dimension != 4 && dimension != 8) {
			PrintlnError("CartesianToSpherical dimension = {} invalid", dimension);
			return GSL_EINVAL;
		}
#endif
		if (spherical[1] = Norm(cartesian + 1); spherical[1] == 0.)
			return GSL_EZERODIV;
		const double r_1 = 1. / spherical[1];
		spherical[0] = cartesian[0];
		spherical[2] = acos(cartesian[3] * r_1);
		if (dimension == 4) {
			spherical[3] = atan2(cartesian[2], cartesian[1]);
			return GSL_SUCCESS;
		}
		spherical[4] = cartesian[4];
		spherical[5] = SBody::Dot(cartesian + 1, cartesian + 5) * r_1;
		if (const double r_xy = Norm(cartesian + 1, 2); r_xy == 0.) {
			spherical[3] = atan2(cartesian[6], cartesian[5]);
			spherical[6] = GSL_SIGN(cartesian[3]) * Norm(cartesian + 5, 2) * r_1;
			spherical[7] = 0;
		} else {
			spherical[3] = atan2(cartesian[2], cartesian[1]);
			spherical[6] = (-cartesian[7] + cartesian[3] * r_1 * spherical[5]) / r_xy;
			spherical[7] = (cartesian[6] * cartesian[1] - cartesian[5] * cartesian[2]) / gsl_pow_2(r_xy);
		}
		return GSL_SUCCESS;
	}
	int CartesianToSpherical(double x[], size_t dimension) {
		double cartesian[dimension];
		copy(x, x + dimension, cartesian);
		return CartesianToSpherical(cartesian, x, dimension);
	}
	int SphericalToCartesian(const double spherical[], double cartesian[], size_t dimension) {
#ifndef GSL_RANGE_CHECK_OFF
		if (dimension != 4 && dimension != 8) {
			PrintlnError("SphericalToCartesian dimension = {}", dimension);
			return GSL_EINVAL;
		}
#endif
		const double sin_theta = abs(sin(spherical[2])), cos_theta = GSL_SIGN(spherical[2]) * cos(spherical[2]), sin_phi = sin(spherical[3]), cos_phi = cos(spherical[3]);
		cartesian[0] = spherical[0];
		cartesian[1] = spherical[1] * sin_theta * cos_phi;
		cartesian[2] = spherical[1] * sin_theta * sin_phi;
		cartesian[3] = spherical[1] * cos_theta;
		if (dimension == 8) {
			cartesian[4] = spherical[4];
			cartesian[5] = spherical[5] * sin_theta * cos_phi + spherical[1] * (cos_theta * cos_phi * spherical[6] - sin_theta * sin_phi * spherical[7]);
			cartesian[6] = spherical[5] * sin_theta * sin_phi + spherical[1] * (cos_theta * sin_phi * spherical[6] + sin_theta * cos_phi * spherical[7]);
			cartesian[7] = spherical[5] * cos_theta - spherical[1] * sin_theta * spherical[6];
		}
		return GSL_SUCCESS;
	}
	int SphericalToCartesian(double x[], size_t dimension) {
		double spherical[dimension];
		copy(x, x + dimension, spherical);
		return SphericalToCartesian(spherical, x, dimension);
	}
	double SphericalAngle(double cos_theta_x, double cos_theta_y, double delta_theta_xy, double delta_phi_xy) {
		const double a = gsl_pow_2(sin(0.5 * delta_theta_xy)) + cos_theta_x * cos_theta_y * gsl_pow_2(sin(0.5 * PhiDifference(delta_phi_xy)));
		return 2. * atan2(SquareRoot(a), SquareRoot(1. - a));
	}
	int OppositeSign(double x, double y) {
		return (x > 0. && y < 0.) || (x < 0. && y > 0.);
	}
	int MapTheta(const double theta_0, double y[]) {
		if (OppositeSign(theta_0, y[2])) {
			y[2] = -y[2];
			y[3] += M_PI;
			y[6] = -y[6];
		} else if (y[2] <= -M_PI_2)
			y[2] += M_PI;
		else if (y[2] > M_PI_2)
			y[2] -= M_PI;
		return GSL_SUCCESS;
	}
	double ModBy2Pi(double phi) {
		return phi - floor(phi / M_2PI) * M_2PI;
	}
	double PhiDifference(double phi) {
		return phi - floor(phi / M_2PI + 0.5) * M_2PI;
	}
	int AMinusPlusB(double a, double b, double &a_minus_b, double &a_plus_b) {
		a_minus_b = a - b;
		a_plus_b = a + b;
		return GSL_SUCCESS;
	}
	double LinearInterpolation(double x, double x0, double x1, double y0, double y1) {
		if (x0 == x1)
			return 0.5 * (y0 + y1);
		return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);
	}
	int LinearInterpolation(double x, double x0, double x1, const double y0[], const double y1[], double y[], size_t size) {
		memset(y, 0, sizeof(double) * size);
		if (x0 == x1) {
			cblas_daxpy(size, 0.5, y1, 1, y, 1);
			cblas_daxpy(size, 0.5, y0, 1, y, 1);
			return GSL_SUCCESS;
		}
		const double x01_1 = 1. / (x1 - x0);
		cblas_daxpy(size, (x - x0) * x01_1, y1, 1, y, 1);
		cblas_daxpy(size, (x1 - x) * x01_1, y0, 1, y, 1);
		return GSL_SUCCESS;
	}
	int InterpolateSphericalPositionToCartesian(double t, double t0, double t1, const double y0[], const double y1[], double y[]) {
		double c0[8], c1[8];
		SphericalToCartesian(y0, c0);
		SphericalToCartesian(y1, c1);
		return LinearInterpolation(t, t0, t1, c0, c1, y, 8);
	}
	double Flux(double luminosity, double magnification, double redshift) {
		return luminosity * magnification / redshift;
	}
	double FluxDensity(double spectral_density, double magnification) {
		return spectral_density * magnification;
	}
	int PolishQuadraticRoot(double a, double b, double roots[], int root_num) {
		for (int i = 0; i < root_num; ++i) {
			for (int j = 32; j > 0; --j) {
				double f = roots[i] + a;
				f = fma(f, roots[i], b);
				double df = fma(2., roots[i], a);
				const double diff = f * df / (df * df - f);
				if (abs(roots[i]) * GSL_DBL_EPSILON < abs(diff))
					roots[i] -= diff;
				else
					break;
			}
		}
		return root_num;
	}
	int PolishCubicRoot(double a, double b, double c, double roots[], int root_num) {
		for (int i = 0; i < root_num; ++i) {
			for (int j = 32; j > 0; --j) {
				double f = roots[i] + a;
				f = fma(f, roots[i], b);
				f = fma(f, roots[i], c);
				double df = fma(3., roots[i], 2. * a);
				df = fma(df, roots[i], b);
				double d2f = fma(6., roots[i], 2. * a);
				const double diff = f * df / (df * df - 0.5 * f * d2f);
				if (abs(roots[i]) * GSL_DBL_EPSILON < abs(diff))
					roots[i] -= diff;
				else
					break;
			}
		}
		return root_num;
	}
	int PolishQuarticRoot(double a, double b, double c, double d, double roots[], int root_num) {
		int polish_state = 0;
		for (int i = 0; i < root_num; ++i) {
			for (int j = 32; j > 0; --j) {
				double f = roots[i] + a;
				f = fma(f, roots[i], b);
				f = fma(f, roots[i], c);
				f = fma(f, roots[i], d);
				double df = fma(4., roots[i], 3. * a);
				df = fma(df, roots[i], 2. * b);
				df = fma(df, roots[i], c);
				double d2f = fma(12., roots[i], 6. * a);
				d2f = fma(d2f, roots[i], 2. * b);
				const double diff = f * df / (df * df - 0.5 * f * d2f);
				if (abs(roots[i]) * GSL_DBL_EPSILON < abs(diff))
					roots[i] -= diff;
				else {
					polish_state |= 1 << i;
					break;
				}
			}
		}
		return polish_state;
	}
	int PolySolveQuarticWithZero(double a, double b, double c, double offset, double roots[]) {
		if (int root_num = gsl_poly_solve_cubic(a, b, c, roots, roots + 1, roots + 2); root_num == 1) {
			if (roots[0] < 0.) {
				roots[0] += offset;
				roots[1] = offset;
			} else {
				roots[1] = roots[0] + offset;
				roots[0] = offset;
			}
			return 2;
		}
		// root_num == 3
		roots[3] = offset;
		if (offset != 0.)
			for (int i = 0; i < 3; ++i)
				roots[i] += offset;
		for (int i = 3; i > 0 && roots[i] < roots[i - 1]; --i)
			swap(roots[i], roots[i - 1]);
		return 4;
	}
	int PolySolveQuartic(double a, double b, double c, double d, double roots[]) {
		if (d == 0.) {
			int root_num = PolySolveQuarticWithZero(a, b, c, 0., roots);
			PolishQuarticRoot(a, b, c, d, roots, root_num);
			return root_num;
		}
		// Now solve x^4 + ax^3 + bx^2 + cx + d = 0.
		const double a_4 = 0.25 * a;
		const double a2_16 = gsl_pow_2(a_4);
		// Let x = y - a/4:
		// Mathematica: Expand[(y - a/4)^4 + a*(y - a/4)^3 + b*(y - a/4)^2 + c*(y - a/4) + d]
		// We now solve the depressed quartic y^4 + py^2 + qy + r = 0.
		const double p = b - 6. * a2_16;
		const double q = c - 2. * b * a_4 + 2. * a * a2_16;
		const double r = d - c * a_4 + b * a2_16 - 3. * a2_16 * a2_16;
		if (r == 0.) {
			int root_num = PolySolveQuarticWithZero(0., p, q, -a_4, roots);
			PolishQuarticRoot(a, b, c, d, roots, root_num);
			return root_num;
		}
		// Biquadratic case:
		if (q == 0.) {
			if (int root_num = gsl_poly_solve_quadratic(1, p, r, roots, roots + 1); root_num == 0)
				return 0;
			if (roots[0] >= 0.) { // roots[1] >= roots[0]
				double root_of_root_0 = sqrt(roots[0]);
				double root_of_root_1 = sqrt(roots[1]);
				roots[0] = -root_of_root_1 - a_4;
				roots[1] = -root_of_root_0 - a_4;
				roots[2] = root_of_root_0 - a_4;
				roots[3] = root_of_root_1 - a_4;
				PolishQuarticRoot(a, b, c, d, roots, 4);
				return 4;
			}
			if (roots[1] >= 0.) {
				double root_of_root_1 = sqrt(roots[1]);
				roots[0] = -root_of_root_1 - a_4;
				roots[1] = root_of_root_1 - a_4;
				PolishQuarticRoot(a, b, c, d, roots, 2);
				return 2;
			}
			return 0;
		}

		// Now split the depressed quartic into two quadratics:
		// y^4 + py^2 + qy + r = (y^2 + sy + u)(y^2 - sy + v) = y^4 + (v+u-s^2)y^2 + s(v - u)y + uv
		// So p = v+u-s^2, q = s(v - u), r = uv.
		// Then (v+u)^2 - (v-u)^2 = 4uv = 4r = (p+s^2)^2 - q^2/s^2.
		// Multiply through by s^2 to get s^2(p+s^2)^2 - q^2 - 4rs^2 = 0, which is a cubic in s^2.
		// Then we let z = s^2, to get
		// z^3 + 2pz^2 + (p^2 - 4r)z - q^2 = 0.
		int z_root_num = gsl_poly_solve_cubic(2. * p, p * p - 4. * r, -q * q, roots, roots + 1, roots + 2);
		// z = s^2, so s = sqrt(z).
		// Hence we require a root > 0, and for the sake of sanity we should take the largest one:
		double largest_root = z_root_num == 1 ? roots[0] : roots[2];
		PolishCubicRoot(2. * p, p * p - 4. * r, -q * q, &largest_root, 1);
		if (largest_root <= 0.) // No real roots:
			return 0;
		const double s = sqrt(largest_root);
		// s is nonzero, because we took care of the biquadratic case.
		const double v = 0.5 * (p + largest_root + q / s);
		const double u = v - q / s;
		// We have q != 0., so u != v
		if (u > v) {
			// Now solve y^2 - sy + v = 0:
			if (gsl_poly_solve_quadratic(1., -s, v, roots, roots + 1) == 0)
				return 0;
			swap(roots[0], roots[1]); // roots[1] > 0 are more convinced, put it to the front
		} else {
			// Now solve y^2 + sy + u = 0:
			if (gsl_poly_solve_quadratic(1., s, u, roots, roots + 1) == 0) // roots[0] < 0 are more convinced
				return 0;
		}
		roots[0] -= a_4;
		roots[1] -= a_4;
		if (int status = PolishQuarticRoot(a, b, c, d, roots, 2); status == 1) { // roots[1] is not convinced, solve the cubic to get another root:
			int rest_root_num = gsl_poly_solve_cubic(a + roots[0], b + (a + roots[0]) * roots[0], -d / roots[0], roots + 1, roots + 2, roots + 3);
			if (rest_root_num == 3 && u < v)
				swap(roots[1], roots[3]); // we want the largest root in this case
			PolishQuarticRoot(a, b, c, d, roots + 1, 1);
		}
		const double prod_r2_r3 = d / (roots[0] * roots[1]);
		const double sum_r2_r3 = prod_r2_r3 * (roots[0] + roots[1]) / (roots[0] * roots[1]);
		if (gsl_poly_solve_quadratic(1., sum_r2_r3, prod_r2_r3, roots + 2, roots + 3) == 0) {
			if (roots[0] > roots[1])
				swap(roots[0], roots[1]);
			return 2;
		}
		PolishQuarticRoot(a, b, c, d, roots + 2, 2);
		sort(roots, roots + 4);
		// const double prod_r1_r2 = d / (roots[0] * roots[3]);
		// const double sum_r1_r2 = prod_r1_r2 * (roots[0] + roots[3]) / (roots[0] * roots[3]);
		// if (gsl_pow_2(sum_r1_r2) < 4. * prod_r1_r2) { // Re-verify the roots:
		// 	swap(roots[1], roots[3]);
		// 	return 2;
		// }
		return 4;
	}
	double EllipticIntegral(int p5, double y, double x, double a5, double b5, double a1, double b1, double a2, double b2, double a3, double b3, double a4, double b4) {
		if (x == y)
			return 0.;
		if (p5 != 0 && p5 != -2 && p5 != -4)
			return GSL_NAN;
		const double d12 = a1 * b2 - a2 * b1, d13 = a1 * b3 - a3 * b1, d14 = a1 * b4 - a4 * b1;
		const double d23 = a2 * b3 - a3 * b2, d24 = a2 * b4 - a4 * b2, d34 = a3 * b4 - a4 * b3;
		const double X1 = sqrt(a1 + b1 * x), X2 = sqrt(a2 + b2 * x), X3 = sqrt(a3 + b3 * x), X4 = sqrt(a4 + b4 * x), X52 = a5 + b5 * x;
		const double Y1 = sqrt(a1 + b1 * y), Y2 = sqrt(a2 + b2 * y), Y3 = sqrt(a3 + b3 * y), Y4 = sqrt(a4 + b4 * y), Y52 = a5 + b5 * y;
		const double U2_12 = gsl_pow_2((X1 * X2 * Y3 * Y4 + Y1 * Y2 * X3 * X4) / (x - y));
		const double U2_13 = U2_12 - d14 * d23;
		const double U2_14 = U2_12 - d13 * d24;
		const double I1 = 2. * gsl_sf_ellint_RF(U2_12, U2_13, U2_14, GSL_PREC_DOUBLE);
		if (p5 == 0)
			return I1;
		const double b52 = gsl_pow_2(b5);
		const double d15 = a1 * b5 - a5 * b1, d15_1 = 1. / d15, d25 = a2 * b5 - a5 * b2, d35 = a3 * b5 - a5 * b3, d45 = a4 * b5 - a5 * b4;
		const double W2 = U2_12 - d13 * d14 * d25 * d15_1;
		const double Q2 = X52 * Y52 * W2 / gsl_pow_2(X1 * Y1);
		const double P2 = Q2 + d25 * d35 * d45 * d15_1;
		const double RC_P2_Q2 = X1 * Y1 == 0. ? 0. : Carlson_RC(P2, Q2);
		const double I3 = 2. * (d12 * d13 * d14 * d15_1 / 3. * Carlson_RJ(U2_12, U2_13, U2_14, W2) + RC_P2_Q2);
		if (p5 == -2)
			return (b5 * I3 - b1 * I1) * d15_1;
		const double I2 = 2. * (d12 * d13 * gsl_sf_ellint_RD(U2_12, U2_13, U2_14, GSL_PREC_DOUBLE) / 3. + X1 * Y1 / (X4 * Y4 * sqrt(U2_14)));
		return -0.5 * d15_1 * (b1 / d15 + b2 / d25 + b3 / d35 + b4 / d45) * b5 * I3 + b52 * d24 * d34 / (2. * d15 * d25 * d35 * d45) * I2 + gsl_pow_2(b1 * d15_1) * (1. - d12 * d13 * b52 / (2. * b1 * b1 * d25 * d35)) * I1 - b52 / (d15 * d25 * d35) * (X1 * X2 * X3 / (X4 * X52) - Y1 * Y2 * Y3 / (Y4 * Y52));
	}

	double EllipticIntegral2Complex(int p5, double y, double x, double a5, double b5, double f, double g, double h, double a1, double b1, double a4, double b4) {
		if (x == y)
			return 0.;
		if (p5 != 0 && p5 != -2 && p5 != -4)
			return GSL_NAN;
		const double b12 = gsl_pow_2(b1), b52 = gsl_pow_2(b5);
		const double X1 = sqrt(a1 + b1 * x), X4 = sqrt(a4 + b4 * x);
		const double Y1 = sqrt(a1 + b1 * y), Y4 = sqrt(a4 + b4 * y);
		const double xi = sqrt(f + (g + h * x) * x), eta = sqrt(f + (g + h * y) * y);
		const double M2 = gsl_pow_2(X1 * Y4 + X4 * Y1) * (gsl_pow_2((xi + eta) / (x - y)) - h);
		const double c2_11 = 2 * (f * b12 - g * a1 * b1 + h * a1 * a1), c2_44 = 2 * (f * b4 * b4 - g * a4 * b4 + h * a4 * a4);
		const double c2_14 = 2. * f * b1 * b4 - g * (a1 * b4 + a4 * b1) + 2. * h * a1 * a4, c2_15 = 2. * f * b1 * b5 - g * (a1 * b5 + a5 * b1) + 2. * h * a1 * a5;
		const double c11 = sqrt(c2_11), c44 = sqrt(c2_44), c11_c44 = c11 * c44;
		const double L2m = max(0., M2 + c2_14 - c11_c44), L2p = max(0., M2 + c2_14 + c11_c44);
		const double I1 = 4. * gsl_sf_ellint_RF(M2, L2m, L2p, GSL_PREC_DOUBLE);
		if (p5 == 0)
			return I1;
		const double X52 = a5 + b5 * x, Y52 = a5 + b5 * y;
		const double c2_55 = 2 * (f * b52 - g * a5 * b5 + h * a5 * a5), c55 = sqrt(c2_55);
		const double d14 = a1 * b4 - a4 * b1, d15 = a1 * b5 - a5 * b1, d15_1 = 1. / d15, d45 = a4 * b5 - a5 * b4;
		const double W2p = M2 + d14 * (c2_15 + c11 * c55) * d15_1;
		const double U = (X1 * X4 * eta + Y1 * Y4 * xi) / (x - y), U2 = gsl_pow_2(U);
		const double W2 = U2 - 0.5 * c2_11 * d45 * d15_1;
		const double Q2 = X52 * Y52 * W2 / gsl_pow_2(X1 * Y1);
		const double P2 = Q2 + 0.5 * c2_55 * d45 * d15_1;
		const double RC_P2_Q2 = X1 * Y1 == 0. ? 0. : Carlson_RC(P2, Q2);
		const double I3 = 2. * (c11 / (3. * c55) * (4. * (W2p - M2) * Carlson_RJ(M2, L2m, L2p, W2p) - 1.5 * I1 + 3. * Carlson_RC(U2, W2)) + RC_P2_Q2);
		if (p5 == -2)
			return (b5 * I3 - b1 * I1) * d15_1;
		const double I2 = 2. * (c11 / (3. * c44) * (4. * (c2_14 + c11_c44) * gsl_sf_ellint_RD(M2, L2m, L2p, GSL_PREC_DOUBLE) - 1.5 * I1 + 3. / U) + X1 * Y1 / (X4 * Y4 * U));
		return -0.5 * d15_1 * ((b1 * b5) / d15 + 2. * b5 * (g * b5 - 2. * h * a5) / c2_55 + (b4 * b5) / d45) * I3 + b52 * 0.5 * c2_44 / (d15 * d45 * c2_55) * I2 + gsl_pow_2(b1 * d15_1) * (1. - 0.5 * c2_11 * b52 / (b12 * c2_55)) * I1 - b52 / (0.5 * d15 * c2_55) * (X1 * xi / (X4 * X52) - Y1 * eta / (Y4 * Y52));
	}
	double EllipticIntegral4Complex(int p5, double y, double x, double a5, double b5, double f1, double g1, double h1, double f2, double g2, double h2) {
		if (x == y)
			return 0.;
		if (p5 != 0 && p5 != -2 && p5 != -4)
			return GSL_NAN;
		const double x_y_1 = 1. / (x - y);
		const double xi12 = f1 + (g1 + h1 * x) * x, eta12 = f1 + (g1 + h1 * y) * y;
		const double xi22 = f2 + (g2 + h2 * x) * x, eta22 = f2 + (g2 + h2 * y) * y;
		const double xi1 = sqrt(xi12), eta1 = sqrt(eta12);
		const double xi2 = sqrt(xi22), eta2 = sqrt(eta22);
		const double xi1p = (g1 + 2. * h1 * x) / (2. * xi1), eta_1p = (g1 + 2. * h1 * y) / (2. * eta1);
		const double B = xi1p * xi2 - eta_1p * eta2;
		const double theta1 = xi12 + eta12 - h1 * gsl_pow_2(x - y);
		const double theta2 = xi22 + eta22 - h2 * gsl_pow_2(x - y);
		const double zeta1 = sqrt(2. * xi1 * eta1 + theta1), zeta2 = sqrt(2. * xi2 * eta2 + theta2);
		const double U = (xi1 * eta2 + xi2 * eta1) * x_y_1, U2 = gsl_pow_2(U);
		const double M = zeta1 * zeta2 * x_y_1;
		const double M2 = gsl_pow_2(M);
		const double delta11_2 = 4. * f1 * h1 - gsl_pow_2(g1), delta12_2 = 2. * (f1 * h2 + f2 * h1) - g1 * g2, delta22_2 = 4. * f2 * h2 - gsl_pow_2(g2);
		const double Delta = sqrt(gsl_pow_2(delta12_2) - delta11_2 * delta22_2);
		const double Delta_m = delta12_2 - Delta, Delta_p = delta12_2 + Delta;
		const double L2m = M2 + Delta_m, L2p = M2 + Delta_p;
		const double RF = gsl_sf_ellint_RF(M2, L2m, L2p, GSL_PREC_DOUBLE);
		if (p5 == 0)
			return 4. * RF;
		const double G = 2. * Delta * Delta_p * gsl_sf_ellint_RD(M2, L2m, L2p, GSL_PREC_DOUBLE) / 3. + Delta / (2. * U) + (delta12_2 * theta1 - delta11_2 * theta2) / (4. * xi1 * eta1 * U);
		const double Sigma = G - Delta_p * RF + B;
		const double alpha15 = 2. * f1 * b5 - g1 * a5, beta15 = g1 * b5 - 2. * h1 * a5;
		const double alpha25 = 2. * f2 * b5 - g2 * a5, beta25 = g2 * b5 - 2. * h2 * a5;
		const double gamma1 = 0.5 * (alpha15 * b5 - beta15 * a5), gamma2 = 0.5 * (alpha25 * b5 - beta25 * a5);
		const double gamma1_1 = 1. / gamma1, gamma2_1 = 1. / gamma2;
		const double Lambda = delta11_2 * gamma2 * gamma1_1;
		const double Omega2 = M2 + Lambda;
		const double psi = 0.5 * (alpha15 * beta25 - alpha25 * beta15);
		const double xi5 = a5 + b5 * x, eta5 = a5 + b5 * y;
		const double X = 0.5 * x_y_1 * (xi5 * (alpha15 + beta15 * y) * eta2 / eta1 + eta5 * (alpha15 + beta15 * x) * xi2 / xi1);
		const double S = 0.5 * (M2 + delta12_2) - U2, S2 = gsl_pow_2(S);
		const double mu = gamma1 * xi5 * eta5 / (xi1 * eta1);
		const double T = mu * S + 2. * gamma1 * gamma2, T2 = gsl_pow_2(T);
		const double V2 = gsl_pow_2(mu) * (S2 + Lambda * U2);
		const double a = S * Omega2 / U + 2. * Lambda * U, a2 = gsl_pow_2(a);
		const double b2 = (S2 / U2 + Lambda) * gsl_pow_2(Omega2);
		const double H = delta11_2 * psi * gsl_pow_2(gamma1_1) * (Carlson_RJ(M2, L2m, L2p, Omega2) / 3. + 0.5 * Carlson_RC(a2, b2)) - X * Carlson_RC(T2, V2);
		if (p5 == -2)
			return -2. * (b5 * H + beta15 * RF * gamma1_1);
		return b5 * (beta15 * gamma1_1 + beta25 * gamma2_1) * H + gsl_pow_2(beta15 * gamma1_1) * RF + gsl_pow_2(b5) * (Sigma - b5 * (xi1 * xi2 / xi5 - eta1 * eta2 / eta5)) * gamma1_1 * gamma2_1;
	}
	double Carlson_RC(double x, double y, gsl_mode_t mode) {
		if (y < 0.)
			return sqrt(x / (x - y)) * gsl_sf_ellint_RC(x - y, -y, mode);
		return gsl_sf_ellint_RC(x, y, mode);
	}
	double Carlson_RJ(double x, double y, double z, double p, gsl_mode_t mode) {
		if (p <= 0.) {
			if ((x >= y && z >= x) || (x <= y && z <= x))
				swap(x, y);
			else if ((z >= y && x >= z) || (z <= y && x <= z))
				swap(y, z);
			const double y_1 = 1. / y, y_p_1 = 1. / (y - p), gamma = y + (z - y) * (y - x) * y_p_1; // make sure gamma > 0.
			return ((gamma - y) * gsl_sf_ellint_RJ(x, y, z, gamma, mode) - 3. * (gsl_sf_ellint_RF(x, y, z, mode) - Carlson_RC(x * z * y_1, p * gamma * y_1))) * y_p_1;
		}
		return gsl_sf_ellint_RJ(x, y, z, p, mode);
	}
} // namespace SBody
