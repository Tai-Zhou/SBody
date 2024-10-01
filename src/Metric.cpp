/**
 * @file Metric.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Metric.h"

#include <cmath>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_vector.h>

#include "IO.h"
#include "Utility.h"

using namespace std;

namespace SBody {
	KerrFastTraceParameters::KerrFastTraceParameters(double a, double u_obs, double u_tgt, double mu_obs, double mu_tgt, double theta_tgt, double phi_tgt) : a(a), a2(gsl_pow_2(a)), u_plus(1. / (1. + sqrt(1. - a2))), u_minus(1. / (1. - sqrt(1. - a2))), u_r(-0.5 / sqrt(1. - a2)), u_obs(u_obs), u_tgt(u_tgt), mu_obs(mu_obs), mu_tgt(mu_tgt), theta_tgt(theta_tgt), phi_tgt(phi_tgt) {}
	int Metric::LocalInertialFrame(const double position[], TimeSystem time, gsl_matrix *coordinate) {
		GslBlock collector;
		gsl_matrix *metric = collector.MatrixAlloc(4, 4), *product = collector.MatrixAlloc(4, 4), *product_LU = collector.MatrixAlloc(4, 4);
		gsl_matrix_set_identity(product);
		MetricTensor(position, metric);
		gsl_vector *coordinate_row, *product_row;
		gsl_permutation *permutation = collector.PermutationAlloc(4);
		int signum;
		for (int i = 0; i < 4; ++i) {
			coordinate_row = collector.VectorAllocRowFromMatrix(coordinate, i);
			product_row = collector.VectorAllocRowFromMatrix(product, i);
			gsl_vector_set_basis(coordinate_row, i);
			if (i == 0) {
				if (time == T)
					copy(position + 5, position + 8, coordinate_row->data + 1);
				else
					copy(position + 4, position + 8, coordinate_row->data);
			} else {
				gsl_matrix_memcpy(product_LU, product);
				gsl_linalg_LU_decomp(product_LU, permutation, &signum);
				gsl_linalg_LU_svx(product_LU, permutation, coordinate_row);
			}
			gsl_vector_scale(coordinate_row, 1. / sqrt(abs(DotProduct(position, coordinate_row->data, coordinate_row->data, 4))));
			gsl_blas_dsymv(CblasUpper, 1., metric, coordinate_row, 0., product_row);
		}
		return isnan(gsl_matrix_get(coordinate, 3, 0)) ? GSL_EDOM : GSL_SUCCESS;
	}
	double Metric::Redshift(const double y[], const double photon[], TimeSystem object_time, TimeSystem photon_time) {
		const double u[4] = {1., y[5], y[6], y[7]}, v[4] = {1., photon[5], photon[6], photon[7]};
		if (object_time == T) {
			if (photon_time == T)
				return -DotProduct(y, u, v, 4) / (y[4] * photon[4]);
			return -DotProduct(y, u, photon + 4, 4) / y[4];
		}
		if (photon_time == T)
			return -DotProduct(y, y + 4, v, 4) / photon[4];
		return -DotProduct(y, y + 4, photon + 4, 4);
	}

	Newton::Newton(int PN) : PN_(PN) {};
	string Newton::Name() {
		return "Newton";
	}
	int Newton::MetricTensor(const double position[], gsl_matrix *metric) { // FIXME: Check
		gsl_matrix_set_zero(metric);
		gsl_matrix_set(metric, 0, 0, -1.);
		gsl_matrix_set(metric, 1, 1, 1.);
		gsl_matrix_set(metric, 2, 2, gsl_pow_2(position[1]));
		gsl_matrix_set(metric, 3, 3, gsl_pow_2(position[1] * sin(position[2])));
		return GSL_SUCCESS;
	}
	double Newton::DotProduct(const double position[], const double x[], const double y[], const size_t dimension) {
		if (dimension == 3)
			return x[1] * y[1] + gsl_pow_2(position[1]) * x[2] * y[2] + gsl_pow_2(position[1] * sin(position[2])) * x[3] * y[3];
		return -x[0] * y[0] + x[1] * y[1] + gsl_pow_2(position[1]) * x[2] * y[2] + gsl_pow_2(position[1] * sin(position[2])) * x[3] * y[3];
	}
	double Newton::DistanceSquare(const double x[], const double y[], const size_t dimension) {
		if (dimension == 3)
			return gsl_pow_2(x[1] - y[1]) + gsl_pow_2(x[1] * (x[2] - y[2])) + gsl_pow_2(x[1] * sin(x[2]) * (x[3] - y[3]));
		return -gsl_pow_2(x[0] - y[0]) + gsl_pow_2(x[1] - y[1]) + gsl_pow_2(x[1] * (x[2] - y[2])) + gsl_pow_2(x[1] * sin(x[2]) * (x[3] - y[3]));
	}
	int Newton::LagrangianToHamiltonian(double y[]) {
		return GSL_FAILURE;
	}
	int Newton::HamiltonianToLagrangian(double y[]) {
		return GSL_FAILURE;
	}
	int Newton::FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_target, const double theta_target, const double phi_target, double &alpha, double &beta, double photon[]) {
		const double sin_theta_target = abs(sin(theta_target)), cos_theta_target = GSL_SIGN(theta_target) * cos(theta_target), sin_phi_target = sin(phi_target), cos_phi_target = cos(phi_target);
		photon[1] = r_target * sin_theta_target * cos_phi_target;
		photon[2] = r_target * sin_theta_target * sin_phi_target;
		photon[3] = r_target * cos_theta_target;
		const double dx = r_observer * sin_theta_observer - photon[1], dz = r_observer * cos_theta_observer - photon[3];
		photon[0] = photon[8] = sqrt(gsl_pow_2(dx) + gsl_pow_2(photon[2]) + gsl_pow_2(dz));
		if (photon[0] == 0.)
			return GSL_EZERODIV;
		alpha = photon[2];
		beta = photon[3] * sin_theta_observer - photon[1] * cos_theta_observer;
		photon[4] = 1.;
		const double distance_1 = 1. / photon[0];
		photon[5] = dx * distance_1;
		photon[6] = photon[2] * distance_1;
		photon[7] = dz * distance_1;
		return CartesianToSpherical(photon);
	}
	double Newton::Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		const double r_1 = 1. / y[1], r_2 = gsl_pow_2(r_1), r_3 = gsl_pow_3(r_1), r_4 = gsl_pow_2(r_2);
		const double v2 = gsl_pow_2(y[5]) + gsl_pow_2(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(sin(y[2]) * y[7])), v4 = gsl_pow_2(v2), v6 = gsl_pow_3(v2), v8 = gsl_pow_2(v4);
		const double rdot = y[5], rdot2 = gsl_pow_2(rdot);
		double E = 0.5 * v2 - r_1;
		if (PN_ & 1)
			E += 0.5 * r_2 + 0.375 * v4 + 1.5 * r_1 * v2;
		if (PN_ & 2)
			E += -0.5 * r_3 + 0.3125 * v6 + 1.75 * r_2 * v2 + 0.5 * r_2 * rdot2 + 2.625 * r_1 * v4;
		if (PN_ & 8)
			E += 0.375 * r_4 + 1.25 * r_3 * v2 + 1.5 * r_3 * rdot2 + 0.2734375 * v8 + 8.4375 * r_2 * v4 + 0.75 * r_2 * v2 * rdot2 + 3.4375 * r_1 * v6;
		return E;
	}
	double Newton::AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		const double r_1 = 1. / y[1], r_2 = gsl_pow_2(r_1), r_3 = gsl_pow_3(r_1);
		const double v_tan2 = gsl_pow_2(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(sin(y[2]) * y[7]));
		const double v2 = gsl_pow_2(y[5]) + v_tan2, v4 = gsl_pow_2(v2), v6 = gsl_pow_3(v2);
		const double rdot = y[5], rdot2 = gsl_pow_2(rdot);
		double eff = 1.;
		if (PN_ & 1)
			eff += 3. * r_1 + 0.5 * v2;
		if (PN_ & 2)
			eff += 3.5 * r_2 + 0.375 * v4 + 3.5 * r_1 * v2;
		if (PN_ & 8)
			eff += 2.5 * r_3 + 0.3125 * v6 + 11.25 * r_2 * v2 + 0.5 * r_2 * rdot2 + 4.125 * r_1 * v4;
		return y[1] * sqrt(v_tan2) * eff;
	}
	double Newton::CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) {
		return gsl_pow_2(AngularMomentum(y, time, dynamics));
	}
	double Newton::Redshift(const double y[], const double photon[], TimeSystem object_time, TimeSystem photon_time) {
		const double delta_epsilon = 1. - 2. / y[1];
		return (1. - DotProduct(y, y + 4, photon + 4, 3) / sqrt(delta_epsilon)) / sqrt(delta_epsilon - DotProduct(y, y + 4, y + 4, 3));
	}
	int Newton::NormalizeTimelikeGeodesic(double y[]) {
		y[4] = 1;
		if (gsl_pow_2(y[5]) + gsl_pow_2(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(sin(y[2]) * y[7])) >= 1)
			return GSL_FAILURE;
		return GSL_SUCCESS;
	}
	int Newton::NormalizeNullGeodesic(double y[], double frequency) {
		y[4] = 1;
		const double v_1 = 1. / sqrt(gsl_pow_2(y[5]) + gsl_pow_2(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(sin(y[2]) * y[7])));
		for (int i = 5; i < 8; ++i)
			y[i] *= v_1;
		return GSL_SUCCESS;
	}
	unique_ptr<Integrator> Newton::GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion) {
		if (time != T || dynamics != LAGRANGIAN || motion != GEODESIC) {
			PrintlnError("{}::GetIntegrator() {}, {}, {} invaild", Name(), time, dynamics, motion);
			return nullptr;
		}
		return make_unique<Integrator>(NewtonTLagrangianGeodesic, Jacobian, const_cast<int *>(&this->PN_));
	}
	PN1::PN1(double fSP) : Newton(1), PN1_(fSP) {};
	unique_ptr<Integrator> PN1::GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion) {
		if (time != T || dynamics != LAGRANGIAN || motion != GEODESIC) {
			PrintlnError("{}::GetIntegrator() {}, {}, {} invaild", Name(), time, dynamics, motion);
			return nullptr;
		}
		return make_unique<Integrator>(PN1TLagrangianGeodesic, Jacobian, const_cast<double *>(&this->PN1_));
	}

	Schwarzschild::Schwarzschild() {}
	string Schwarzschild::Name() {
		return "Schwarzschild";
	}
	int Schwarzschild::MetricTensor(const double position[], gsl_matrix *metric) {
		gsl_matrix_set_zero(metric);
		gsl_matrix_set(metric, 0, 0, -(1. - 2. / position[1]));
		gsl_matrix_set(metric, 1, 1, position[1] / (position[1] - 2.));
		gsl_matrix_set(metric, 2, 2, gsl_pow_2(position[1]));
		gsl_matrix_set(metric, 3, 3, gsl_pow_2(position[1] * sin(position[2])));
		return position[1] == 2. ? GSL_EZERODIV : GSL_SUCCESS;
	}
	double Schwarzschild::DotProduct(const double position[], const double x[], const double y[], const size_t dimension) {
		if (dimension == 3)
			return position[1] * x[1] * y[1] / (position[1] - 2.) + gsl_pow_2(position[1]) * x[2] * y[2] + gsl_pow_2(position[1] * sin(position[2])) * x[3] * y[3];
		return -(1. - 2. / position[1]) * x[0] * y[0] + position[1] * x[1] * y[1] / (position[1] - 2.) + gsl_pow_2(position[1]) * x[2] * y[2] + gsl_pow_2(position[1] * sin(position[2])) * x[3] * y[3];
	}
	double Schwarzschild::DistanceSquare(const double x[], const double y[], const size_t dimension) {
		if (dimension == 3)
			return x[1] * gsl_pow_2(x[1] - y[1]) / (x[1] - 2.) + gsl_pow_2(x[1] * (x[2] - y[2])) + gsl_pow_2(x[1] * sin(x[2]) * PhiDifference(x[3] - y[3]));
		return -(1. - 2. / x[1]) * gsl_pow_2(x[0] - y[0]) + x[1] * gsl_pow_2(x[1] - y[1]) / (x[1] - 2.) + gsl_pow_2(x[1] * (x[2] - y[2])) + gsl_pow_2(x[1] * sin(x[2]) * PhiDifference(x[3] - y[3]));
	}
	int Schwarzschild::LagrangianToHamiltonian(double y[]) {
		const double dt_dtau = 1. / y[4], r2 = gsl_pow_2(y[1]);
		y[4] = (y[4] - 1. + 2. / y[1]) * dt_dtau; // 1 + p_t
		y[5] *= y[1] / (y[1] - 2.) * dt_dtau;
		y[6] *= r2 * dt_dtau;
		y[7] *= r2 * gsl_pow_2(sin(y[2])) * dt_dtau;
		return GSL_SUCCESS;
	}
	int Schwarzschild::HamiltonianToLagrangian(double y[]) {
		const double r_1 = 1. / y[1], g11_1 = 1. - 2. * r_1;
		y[4] = -g11_1 / (y[4] - 1.);
		y[5] *= g11_1 * y[4];
		y[6] *= gsl_pow_2(r_1) * y[4];
		y[7] *= gsl_pow_2(r_1 / sin(y[2])) * y[4];
		return GSL_SUCCESS;
	}
	int Schwarzschild::FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_target, const double theta_target, const double phi_target, double &alpha, double &beta, double photon[]) {
		const double sin_theta_target = abs(sin(theta_target)), cos_theta_target = GSL_SIGN(theta_target) * cos(theta_target), sin_phi_target = sin(phi_target), cos_phi_target = cos(phi_target);
		const double cos_observer_target = sin_theta_observer * sin_theta_target * cos_phi_target + cos_theta_observer * cos_theta_target;
		const double delta_phi = acos(cos_observer_target), sin_observer_target = sqrt(1. - gsl_pow_2(cos_observer_target));
		const double u0 = 1. / r_observer, u1 = 1. / r_target;
		const double g11_1 = 1. - 2. * u1;
		if (delta_phi == 0.) {
			alpha = 0.;
			beta = 0.;
			photon[0] = (r_target - r_observer + 4. * log(r_target * u0) + 4. * (u1 - u0)) * r_observer / (r_observer - 2.);
			photon[1] = r_target;
			photon[2] = theta_target;
			photon[3] = phi_target;
			photon[4] = r_target * (r_observer - 2.) / (r_observer * (r_target - 2.));
			photon[5] = (r_target - 2.) * u1;
			photon[6] = 0.;
			photon[7] = 0.;
			photon[8] = r_target - r_observer + 2. * log(r_target * u0);
			return GSL_SUCCESS;
		}
		double impact_parameter_upper_limit = r_target / sqrt(g11_1), turning_phi;
		if (double x0, x1, x2; u1 * 3. < 1. && gsl_poly_solve_cubic(-0.5, 0., 0.5 / gsl_pow_2(impact_parameter_upper_limit), &x0, &x1, &x2) == 3)
			turning_phi = M_SQRT1_2 * EllipticIntegral(0, u0, u1, 0., 0., -x0, 1., max(x1, u1), -1., x2, -1.); // prevent x1 < u1 < x1 + EPSILON
		else {
			impact_parameter_upper_limit = M_SQRT27 - EPSILON;
			turning_phi = M_PI; // make delta_phi <= turning_phi
		}
		const double integrate_parameters[4] = {u0, u1, delta_phi, turning_phi};
		gsl_function impact_function{
			[](double x, void *params) -> double {
				const double *p = static_cast<double *>(params);
				if (x == 0.)
					return -p[2];
				if (double x0, x1, x2; gsl_poly_solve_cubic(-0.5, 0., 0.5 / gsl_pow_2(x), &x0, &x1, &x2) == 3) {
					if (x1 < p[1]) // x1 < u1 < x1 + EPSILON
						return M_SQRT1_2 * EllipticIntegral(0, p[0], p[1], 0., 0., -x0, 1., p[1], -1., x2, -1.) - p[2];
					else if (p[2] > p[3]) // u1 -> turning point -> u1 -> u0
						return M_SQRT2 * EllipticIntegral(0., p[1], x1, 0., 0., -x0, 1., x1, -1., x2, -1.) + M_SQRT1_2 * EllipticIntegral(0, p[0], p[1], 0., 0., -x0, 1., x1, -1., x2, -1.) - p[2];
					else // u1 -> u0
						return M_SQRT1_2 * EllipticIntegral(0, p[0], p[1], 0., 0., -x0, 1., x1, -1., x2, -1.) - p[2];
				} else // impact_parameter < sqrt(27)
					return M_SQRT1_2 * EllipticIntegral2Complex(0, p[0], p[1], 0., 0., -0.5 / (gsl_pow_2(x) * x0), x0 - 0.5, 1., -x0, 1.) - p[2];
			},
			const_cast<double *>(integrate_parameters)};
		FunctionSolver impact_solver;
		if (int status = impact_solver.Set(&impact_function, delta_phi > turning_phi ? M_SQRT27 + EPSILON : 0., impact_parameter_upper_limit); status != GSL_SUCCESS)
			return status;
		if (int status = impact_solver.Solve(); status != GSL_SUCCESS)
			return status;
		alpha = impact_solver.Root() / sin_observer_target * sin_theta_target * sin_phi_target;
		beta = impact_solver.Root() / sin_observer_target * (cos_theta_target * sin_theta_observer - sin_theta_target * cos_phi_target * cos_theta_observer);
		if (double x0, x1, x2; gsl_poly_solve_cubic(-0.5, 0., 0.5 / gsl_pow_2(impact_solver.Root()), &x0, &x1, &x2) == 3) {
			if (x1 < u1) { // x1 < u1 < x1 + EPSILON
				const double ellip_int_4 = EllipticIntegral(-4, u0, u1, 0., 1., u1, -1., -x0, 1., x2, -1.);
				photon[0] = -M_SQRT1_2 * ellip_int_4 / (impact_solver.Root() * (1. - 2. * u0));
				photon[8] = -M_SQRT1_2 * (ellip_int_4 + 2. * EllipticIntegral(-2, u0, u1, 0., 1., -x0, 1., u1, -1., x2, -1.) + 4. * EllipticIntegral(-2, u0, u1, 1., -2., -x0, 1., u1, -1., x2, -1.)) / impact_solver.Root();
			} else if (delta_phi > turning_phi) { // u1 -> turning point -> u1 -> u0
				const double ellip_int_4 = 2. * EllipticIntegral(-4, u1, x1, 0., 1., -x0, 1., x1, -1., x2, -1.) + EllipticIntegral(-4, u0, u1, 0., 1., -x0, 1., x1, -1., x2, -1.);
				photon[0] = -M_SQRT1_2 * ellip_int_4 / (impact_solver.Root() * (1. - 2. * u0));
				photon[8] = -M_SQRT1_2 * (ellip_int_4 + 4. * EllipticIntegral(-2, u1, x1, 0., 1., -x0, 1., x1, -1., x2, -1.) + 2. * EllipticIntegral(-2, u0, u1, 0., 1., -x0, 1., x1, -1., x2, -1.) + 8. * EllipticIntegral(-2, u1, x1, 1., -2., -x0, 1., x1, -1., x2, -1.) + 4. * EllipticIntegral(-2, u0, u1, 1., -2., -x0, 1., x1, -1., x2, -1.)) / impact_solver.Root();
			} else { // u1 -> u0
				const double ellip_int_4 = EllipticIntegral(-4, u0, u1, 0., 1., -x0, 1., x1, -1., x2, -1.);
				photon[0] = -M_SQRT1_2 * ellip_int_4 / (impact_solver.Root() * (1. - 2. * u0));
				photon[8] = -M_SQRT1_2 * (ellip_int_4 + 2. * EllipticIntegral(-2, u0, u1, 0., 1., -x0, 1., x1, -1., x2, -1.) + 4. * EllipticIntegral(-2, u0, u1, 1., -2., -x0, 1., x1, -1., x2, -1.)) / impact_solver.Root();
			}
		} else { // impact_parameter < sqrt(27)
			const double ellip_int_4 = EllipticIntegral2Complex(-4, u0, u1, 0., 1., -0.5 / (gsl_pow_2(impact_solver.Root()) * x0), x0 - 0.5, 1., -x0, 1.);
			photon[0] = -M_SQRT1_2 * ellip_int_4 / (impact_solver.Root() * (1. - 2. * u0));
			photon[8] = -M_SQRT1_2 * (ellip_int_4 + 2. * EllipticIntegral2Complex(-2, u0, u1, 0., 1., -0.5 / (gsl_pow_2(impact_solver.Root()) * x0), x0 - 0.5, 1., -x0, 1.) + 4. * EllipticIntegral2Complex(-2, u0, u1, 1., -2., -0.5 / (gsl_pow_2(impact_solver.Root()) * x0), x0 - 0.5, 1., -x0, 1.)) / impact_solver.Root();
		}
		photon[1] = r_target;
		photon[2] = theta_target;
		photon[3] = phi_target;
		photon[4] = g11_1 / (1. - 2. * u0); // dt/d\tau = 1. when r = r_observer
		if (delta_phi > turning_phi)
			photon[5] = -g11_1 * SquareRoot(1. - g11_1 * gsl_pow_2(u1 * impact_solver.Root()));
		else
			photon[5] = g11_1 * SquareRoot(1. - g11_1 * gsl_pow_2(u1 * impact_solver.Root()));
		// the direction of the angular momentum is [alpha * cos_theta_observer, beta, -alpha * sin_theta_observer],
		// the y component should have the same sign as the component perpendicular to [cos(phi), sin(phi), 0],
		// which is also the component along the [-sin(phi), cos(phi), 0] direction.
		if (alpha == 0.) {
			photon[6] = GSL_SIGN(beta * cos_phi_target - alpha * cos_theta_observer * sin_phi_target) * g11_1 * gsl_pow_2(u1) * impact_solver.Root();
			photon[7] = 0.;
		} else {
			photon[6] = GSL_SIGN(beta * cos_phi_target - alpha * cos_theta_observer * sin_phi_target) * g11_1 * gsl_pow_2(u1) * SquareRoot(gsl_pow_2(impact_solver.Root()) - gsl_pow_2(alpha * sin_theta_observer / sin_theta_target));
			photon[7] = -g11_1 * alpha * sin_theta_observer * gsl_pow_2(u1 / sin_theta_target);
		}
		return GSL_SUCCESS;
	}
	double Schwarzschild::Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			if (time == T)
				return (y[1] - 2.) / (y[1] * y[4]);
			// time == TAU
			return (1. - 2. / y[1]) * y[4];
		} // dynamics == HAMILTONIAN
		return 1. - y[4];
	}
	double Schwarzschild::AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			if (time == T)
				return gsl_pow_2(y[1] * sin(y[2])) * y[7] / y[4];
			// time == TAU
			return gsl_pow_2(y[1] * sin(y[2])) * y[7];
		} // dynamics == HAMILTONIAN
		return y[7];
	}
	double Schwarzschild::CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			if (time == T)
				return gsl_pow_4(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(y[7] * cos(y[2]) * sin(y[2]))) / gsl_pow_2(y[4]);
			// time == TAU
			return gsl_pow_4(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(y[7] * cos(y[2]) * sin(y[2])));
		} // dynamics == HAMILTONIAN
		return gsl_pow_2(y[6]) + gsl_pow_2(y[7] / tan(y[2]));
	}
	int Schwarzschild::NormalizeTimelikeGeodesic(double y[]) {
		const double g11_1 = 1. - 2. / y[1];
		if (g11_1 <= 0)
			return 1;
		y[4] = sqrt(g11_1 - (gsl_pow_2(y[5]) / g11_1 + gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
		return isnan(y[4]) ? GSL_EDOM : GSL_SUCCESS;
	}
	int Schwarzschild::NormalizeNullGeodesic(double y[], double frequency) {
		const double g11_1 = 1. - 2. / y[1];
		if (g11_1 <= 0)
			return 1;
		const double coefficient = GSL_SIGN(frequency) * g11_1 / sqrt(gsl_pow_2(y[5]) + g11_1 * (gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
		if (isnan(coefficient))
			return GSL_EDOM;
		y[4] = frequency;
		y[5] *= coefficient;
		y[6] *= coefficient;
		y[7] *= coefficient;
		return GSL_SUCCESS;
	}
	unique_ptr<Integrator> Schwarzschild::GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion) {
		if (time == T) {
			if (dynamics == LAGRANGIAN) {
				if (motion == GEODESIC)
					return make_unique<Integrator>(SchwarzschildTLagrangianGeodesic, Jacobian);
				else if (motion == CIRCULAR)
					return make_unique<Integrator>(SchwarzschildTLagrangianCircular, Jacobian);
				else if (motion == HELICAL)
					return make_unique<Integrator>(SchwarzschildTLagrangianHelical, Jacobian);
			} else if (dynamics == HAMILTONIAN && motion == GEODESIC)
				return make_unique<Integrator>(SchwarzschildTHamiltonianGeodesic, Jacobian);
		} else if (time == TAU && dynamics == LAGRANGIAN && motion == GEODESIC)
			return make_unique<Integrator>(SchwarzschildTauLagrangianGeodesic, Jacobian);
		PrintlnError("{}::GetIntegrator() {}, {}, {} invaild", Name(), time, dynamics, motion);
		return nullptr;
	}
	ReissnerNordstrom::ReissnerNordstrom(double charge) : r_Q_(0.5 * sqrt(M_1_PI) * charge), r_Q2_(r_Q_ * r_Q_), r_Q4_(r_Q2_ * r_Q2_) {}
	string ReissnerNordstrom::Name() {
		return "Reissner-Nordstrom";
	}
	int ReissnerNordstrom::MetricTensor(const double position[], gsl_matrix *metric) {
		const double r_1 = 1. / position[1];
		gsl_matrix_set_zero(metric);
		gsl_matrix_set(metric, 0, 0, -1. + (2. - r_Q2_ * r_1) * r_1);
		gsl_matrix_set(metric, 1, 1, -1. / gsl_matrix_get(metric, 0, 0));
		gsl_matrix_set(metric, 2, 2, gsl_pow_2(position[1]));
		gsl_matrix_set(metric, 3, 3, gsl_pow_2(position[1] * sin(position[2])));
		return position[1] == 2. ? GSL_EZERODIV : GSL_SUCCESS;
	}
	double ReissnerNordstrom::DotProduct(const double position[], const double x[], const double y[], const size_t dimension) {
		const double r_1 = 1. / position[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
		if (dimension == 3)
			return x[1] * y[1] / (g11_1) + gsl_pow_2(position[1]) * x[2] * y[2] + gsl_pow_2(position[1] * sin(position[2])) * x[3] * y[3];
		return -g11_1 * x[0] * y[0] + x[1] * y[1] / g11_1 + gsl_pow_2(position[1]) * x[2] * y[2] + gsl_pow_2(position[1] * sin(position[2])) * x[3] * y[3];
	}
	double ReissnerNordstrom::DistanceSquare(const double x[], const double y[], const size_t dimension) {
		const double r_1 = 1. / x[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
		if (dimension == 3)
			return gsl_pow_2(x[1] - y[1]) / g11_1 + gsl_pow_2(x[1]) * (gsl_pow_2(x[2] - y[2]) + gsl_pow_2(sin(x[2]) * PhiDifference(x[3] - y[3])));
		return -g11_1 * gsl_pow_2(x[0] - y[0]) + gsl_pow_2(x[1] - y[1]) / g11_1 + gsl_pow_2(x[1]) * (gsl_pow_2(x[2] - y[2]) + gsl_pow_2(sin(x[2]) * PhiDifference(x[3] - y[3])));
	}
	int ReissnerNordstrom::LagrangianToHamiltonian(double y[]) {
		const double dt_dtau = 1. / y[4], r_1 = 1. / y[1], r2 = gsl_pow_2(y[1]);
		const double rs_rQ2 = (2. - r_Q2_ * r_1) * r_1;
		y[4] = (y[4] - 1. + rs_rQ2) * dt_dtau; // 1 + p_t
		y[5] *= dt_dtau / (1. - rs_rQ2);
		y[6] *= r2 * dt_dtau;
		y[7] *= r2 * gsl_pow_2(sin(y[2])) * dt_dtau;
		return GSL_SUCCESS;
	}
	int ReissnerNordstrom::HamiltonianToLagrangian(double y[]) {
		const double r_1 = 1. / y[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
		y[4] = -g11_1 / (y[4] - 1.);
		y[5] *= g11_1 * y[4];
		y[6] *= gsl_pow_2(r_1) * y[4];
		y[7] *= gsl_pow_2(r_1 / sin(y[2])) * y[4];
		return GSL_SUCCESS;
	}
	int ReissnerNordstrom::FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_target, const double theta_target, const double phi_target, double &alpha, double &beta, double photon[]) {
		return GSL_FAILURE;
	}
	double ReissnerNordstrom::Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			const double r_1 = 1. / y[1];
			if (time == T)
				return (1. - (2. - r_Q2_ * r_1) * r_1) / y[4];
			// time == TAU
			return (1. - (2. - r_Q2_ * r_1) * r_1) * y[4];
		} // dynamics == HAMILTONIAN
		return 1. - y[4];
	}
	double ReissnerNordstrom::AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			if (time == T)
				return gsl_pow_2(y[1] * sin(y[2])) * y[7] / y[4];
			// time == TAU
			return gsl_pow_2(y[1] * sin(y[2])) * y[7];
		} // dynamics == HAMILTONIAN
		return y[7];
	}
	double ReissnerNordstrom::CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			if (time == T)
				return gsl_pow_4(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(y[7] * cos(y[2]) * sin(y[2]))) / gsl_pow_2(y[4]);
			// time == TAU
			return gsl_pow_4(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(y[7] * cos(y[2]) * sin(y[2])));
		} // dynamics == HAMILTONIAN
		return gsl_pow_2(y[6]) + gsl_pow_2(y[7] / tan(y[2]));
	}
	int ReissnerNordstrom::NormalizeTimelikeGeodesic(double y[]) {
		const double r_1 = 1. / y[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
		if (g11_1 <= 0)
			return 1;
		y[4] = sqrt(g11_1 - (gsl_pow_2(y[5]) / g11_1 + gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
		return isnan(y[4]) ? GSL_EDOM : GSL_SUCCESS;
	}
	int ReissnerNordstrom::NormalizeNullGeodesic(double y[], double frequency) {
		const double r_1 = 1. / y[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
		if (g11_1 <= 0)
			return 1;
		const double coefficient = GSL_SIGN(frequency) * g11_1 / sqrt(gsl_pow_2(y[5]) + g11_1 * (gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
		if (isnan(coefficient))
			return GSL_EDOM;
		y[4] = frequency;
		y[5] *= coefficient;
		y[6] *= coefficient;
		y[7] *= coefficient;
		return GSL_SUCCESS;
	}
	unique_ptr<Integrator> ReissnerNordstrom::GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion) {
		if (time == T)
			if (dynamics == LAGRANGIAN)
				if (motion == GEODESIC)
					return make_unique<Integrator>(LagrangianGeodesic, Jacobian);
		PrintlnError("{}::GetIntegrator() {}, {}, {} invaild", Name(), time, dynamics, motion);
		return nullptr;
	}
	int ReissnerNordstrom::LagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		ReissnerNordstrom *metric = static_cast<ReissnerNordstrom *>(params);
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = y[5]; // dr/dt
		dydt[2] = y[6]; // d\theta/dt
		dydt[3] = y[7]; // d\phi/dt
		const double r = y[1], r_1 = 1. / r, r_2 = gsl_pow_2(r_1), r_Q2 = metric->r_Q2_;
		const double sin_theta = abs(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
		const double g11_1 = 1. - (2. - r_Q2 * r_1) * r_1, g11 = 1. / g11_1;
		const double m_r_Q_r2 = (1. - r_Q2 * r_1) * r_2, m_r_Q_r2_g11 = m_r_Q_r2 * g11;
		dydt[4] = 2. * m_r_Q_r2_g11 * y[5] * y[4];
		dydt[5] = -g11_1 * m_r_Q_r2 + 3. * m_r_Q_r2_g11 * gsl_pow_2(y[5]) + r * g11_1 * (gsl_pow_2(y[6]) + gsl_pow_2(sin_theta * y[7]));
		dydt[6] = -2. * (g11_1 * r_1 - m_r_Q_r2) * g11 * y[5] * y[6] + sin_theta * cos_theta * gsl_pow_2(y[7]);
		if (y[7] == 0.)
			dydt[7] = 0.;
		else if (sin_theta == 0.)
			return GSL_EZERODIV;
		else
			dydt[7] = -2. * ((r_1 - m_r_Q_r2_g11) * y[5] + cos_theta / sin_theta * y[6]) * y[7];
		return GSL_SUCCESS;
	}
	Kerr::Kerr(double spin) : a_(spin), a2_(a_ * a_), a4_(a2_ * a2_) {}
	string Kerr::Name() {
		return "Kerr";
	}
	int Kerr::MetricTensor(const double position[], gsl_matrix *metric) {
		const double r = position[1], r2 = gsl_pow_2(r), Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(position[2])), sin2_theta = gsl_pow_2(sin(position[2]));
		if (rho2 == 0. || Delta == 0.)
			return GSL_EZERODIV;
		const double r_rho_2 = 2. * r / rho2;
		gsl_matrix_set_zero(metric);
		gsl_matrix_set(metric, 0, 0, -(1. - r_rho_2));
		gsl_matrix_set(metric, 0, 3, -r_rho_2 * a_ * sin2_theta);
		gsl_matrix_set(metric, 1, 1, rho2 / Delta);
		gsl_matrix_set(metric, 2, 2, rho2);
		gsl_matrix_set(metric, 3, 0, -r_rho_2 * a_ * sin2_theta);
		gsl_matrix_set(metric, 3, 3, (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * sin2_theta);
		return GSL_SUCCESS;
	}
	double Kerr::DotProduct(const double position[], const double x[], const double y[], const size_t dimension) {
		const double r = position[1], r2 = gsl_pow_2(r), Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(position[2])), r_rho_2 = 2. * r / rho2, sin2_theta = gsl_pow_2(sin(position[2]));
		if (dimension == 3)
			return (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * gsl_pow_2(sin2_theta)) * x[3] * y[3];
		return (r_rho_2 - 1.) * x[0] * y[0] - r_rho_2 * a_ * sin2_theta * (x[0] * y[3] + x[3] * y[0]) + (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * gsl_pow_2(sin2_theta)) * x[3] * y[3];
	}
	double Kerr::DistanceSquare(const double x[], const double y[], const size_t dimension) {
		const double r = x[1], r2 = gsl_pow_2(r), d0 = x[0] - y[0], d3 = PhiDifference(x[3] - y[3]);
		const double Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(x[2])), r_rho_2 = 2. * r / rho2, sin2_theta = gsl_pow_2(sin(x[2]));
		if (dimension == 3)
			return (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + ((r2 + a2_) * sin2_theta + 2. * r * a2_ * gsl_pow_2(sin2_theta) / rho2) * gsl_pow_2(d3);
		return (r_rho_2 - 1.) * gsl_pow_2(d0) - 2. * r_rho_2 * a_ * sin2_theta * d0 * d3 + (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * gsl_pow_2(sin2_theta)) * gsl_pow_2(d3);
	}
	int Kerr::LagrangianToHamiltonian(double y[]) {
		const double dt_dtau = 1. / y[4], r2 = gsl_pow_2(y[1]), sin2_theta = gsl_pow_2(sin(y[2]));
		const double Delta = r2 - 2. * y[1] + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		if (rho2 == 0. || Delta == 0.)
			return GSL_EZERODIV;
		const double r_rho_2 = 2. * y[1] / rho2;
		y[4] = (y[4] - 1. + r_rho_2 * (1. - a_ * sin2_theta * y[7])) * dt_dtau; // 1 + p_t
		y[5] *= rho2 / Delta * dt_dtau;
		y[6] *= rho2 * dt_dtau;
		y[7] = (-r_rho_2 * a_ + (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * y[7]) * sin2_theta * dt_dtau;
		return GSL_SUCCESS;
	}
	int Kerr::HamiltonianToLagrangian(double y[]) {
		const double pt = y[4] - 1., r2 = gsl_pow_2(y[1]), sin2_theta = gsl_pow_2(sin(y[2]));
		const double Delta = r2 - 2. * y[1] + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		if (rho2 == 0. || Delta == 0.)
			return GSL_EZERODIV;
		const double rho_2 = 1. / rho2, r_rho_2 = 2. * y[1] * rho_2;
		y[4] = -Delta / ((Delta + r_rho_2 * (a2_ + r2)) * pt + r_rho_2 * a_ * y[7]);
		y[5] *= Delta * y[4] * rho_2;
		y[6] *= y[4] * rho_2;
		if (y[7] == 0.)
			y[7] = -r_rho_2 * a_ * pt / Delta * y[4];
		else if (sin2_theta == 0.)
			return GSL_EZERODIV;
		else
			y[7] = (-r_rho_2 * a_ * pt + (1. - r_rho_2) / sin2_theta * y[7]) / Delta * y[4];
		return GSL_SUCCESS;
	}
	int Kerr::MuFSolver(int alpha_1, double M_plus, double M_minus, double mu_plus, double mu_minus, double q2, double I_u, double int_to_mu_plus, double int_mu_full_turn, double A, double x, double k, double a, double mu_obs, double &mu_f) {
		if (q2 == 0.) {
			const double mu0_mu_plus = mu_obs / mu_plus;
			mu_f = mu_plus / cosh(abs(a * mu_plus) * I_u - alpha_1 * GSL_SIGN(mu_obs) * log((1. + sqrt(1. - gsl_pow_2(mu0_mu_plus))) / mu0_mu_plus));
			return GSL_SUCCESS;
		}
		// Here the first turning point is determined by the `beta`.
		const int N = alpha_1 > 0 ? ceil((I_u - int_to_mu_plus) / int_mu_full_turn) : floor((I_u + int_to_mu_plus) / int_mu_full_turn); // s_mu = 1
		const int alpha_2 = N & 1 ? -alpha_1 : alpha_1;
		const int alpha_3 = 2 * floor((2 * N + 3 - alpha_1) / 4) - 1;
		const double I = (I_u - alpha_1 * int_to_mu_plus - alpha_3 * int_mu_full_turn) / alpha_2;
		if (M_minus > 0.) {
			// const double A = mu_plus * abs(param->a);
			// const double x = sqrt((M_plus - gsl_pow_2(param->mu_obs)) / (M_plus - M_minus));
			// const double k = sqrt(1. - gsl_pow_2(mu_minus / mu_plus));
			const double J = mu_plus * abs(a) * I;
			double sn, cn, dn;
			if (int status = gsl_sf_elljac_e(J, k, &sn, &cn, &dn); status != GSL_SUCCESS)
				return status;
			mu_f = mu_minus / dn;
		} else {
			// const double A = sqrt(M_plus - M_minus) * abs(param->a);
			// const double x = sqrt(1. - gsl_pow_2(param->mu_obs / mu_plus));
			// const double k = sqrt(1. + M_minus / (gsl_pow_2(mu_plus) - M_minus));
			const double J = sqrt(gsl_pow_2(mu_plus) - M_minus) * abs(a) * I;
			double sn, cn, dn;
			if (int status = gsl_sf_elljac_e(J, k, &sn, &cn, &dn); status != GSL_SUCCESS)
				return status;
			mu_f = mu_minus * cn;
		}
		return GSL_SUCCESS;
	}
	int Kerr::DifferenceOfUMuPhi(const gsl_vector *alpha_beta, void *params, gsl_vector *delta_u_mu_phi) {
		const auto *param = static_cast<KerrFastTraceParameters *>(params);
		// The photon in the observer's frame has the tetrad velocity: [1, r / R, beta / R, -alpha / R], where R = sqrt(r^2 + alpha^2 + beta^2).
		const double alpha = gsl_vector_get(alpha_beta, 0);
		const double beta = gsl_vector_get(alpha_beta, 1);
		const double l = gsl_vector_get(alpha_beta, 0), l2 = gsl_pow_2(l); // FIXME:
		const double q2 = gsl_vector_get(alpha_beta, 1);				   // FIXME:
		const double c = param->a2 - l2 - q2, d = 2. * (gsl_pow_2(param->a - l) + q2), e = -param->a * q2;
		// U=1−2u+[a^2−q^2−l^2]u^2+2[(a−l)^2+q^2]u^3−a^2q^2u^4
		double u0, u1, u2, u3, I_u_0, I_u_1 = GSL_NAN, mu_f_0, mu_f_1, phi_f_0, phi_f_1;
		double int_u_plus_0, int_u_plus_1;
		double int_u_minus_0, int_u_minus_1;
		double int_u_0, int_u_1;
		double int_u2_0, int_u2_1;
		if (q2 == 0.) {			 // U=1−2u+[a^2−l^2]u^2+2(a−l)^2u^3
			if (param->a == l) { // U=1−2u
				I_u_0 = sqrt(1. - 2. * param->u_obs) - sqrt(1. - 2. * param->u_tgt);
				I_u_1 = 0.;
				const double x_0 = sqrt(1. - 2. * param->u_obs), x_1 = sqrt(1. - 2. * param->u_tgt);
				const double sqrt_2_u_plus = sqrt(1. - 2. * param->u_plus);
				const double sqrt_2_u_minus = sqrt(1. - 2. * param->u_minus);
				int_u_plus_0 = param->u_plus * 2. / sqrt_2_u_plus * (log(sqrt_2_u_plus + x_1) - log(sqrt_2_u_plus + x_0) + log(sqrt_2_u_plus - x_1) - log(sqrt_2_u_plus - x_0));		// TODO: Check
				int_u_minus_0 = param->u_minus * 2. / sqrt_2_u_minus * (log(sqrt_2_u_minus + x_1) - log(sqrt_2_u_minus + x_0) + log(sqrt_2_u_minus - x_1) - log(sqrt_2_u_minus - x_0)); // TODO: Check
			} else {
				const double B_1 = 0.5 / gsl_pow_2(param->a - l);
				if (int root_num_U = gsl_poly_solve_cubic((param->a2 - l2) * B_1, -2. * B_1, B_1, &u0, &u1, &u2); root_num_U == 1) {
					const double f = 2. * u0 / gsl_pow_2(param->a - l);
					I_u_0 = EllipticIntegral2Complex(0, param->u_obs, param->u_tgt, 1., 0., f, f / u0, 1., -u0, 1.);
					I_u_1 = 0.;
					int_u_plus_0 = EllipticIntegral2Complex(-2, param->u_obs, param->u_tgt, -1., 1. / param->u_plus, f, f / u0, 1., -u0, 1.);
					int_u_minus_0 = EllipticIntegral2Complex(-2, param->u_obs, param->u_tgt, -1., 1. / param->u_minus, f, f / u0, 1., -u0, 1.);
				} else {
					I_u_0 = EllipticIntegral(0, param->u_obs, u1, 1., 0., -u0, 1., u1, -1., u2, -1.);
					int_u_plus_0 = EllipticIntegral(-2, param->u_obs, u1, -1., 1. / param->u_plus, -u0, 1., u1, -1., u2, -1.);
					int_u_minus_0 = EllipticIntegral(-2, param->u_obs, u1, -1., 1. / param->u_minus, -u0, 1., u1, -1., u2, -1.);
					if (u1 >= param->u_tgt) {
						I_u_1 = I_u_0 * 2.;
						I_u_0 -= EllipticIntegral(0, param->u_tgt, u1, 1., 0., -u0, 1., u1, -1., u2, -1.);
						I_u_1 -= I_u_0;
						int_u_plus_1 = int_u_plus_0 * 2.;
						int_u_plus_0 -= EllipticIntegral(-2, param->u_tgt, u1, -1., 1. / param->u_plus, -u0, 1., u1, -1., u2, -1.);
						int_u_plus_1 -= int_u_plus_0;
						int_u_minus_1 = int_u_minus_0 * 2.;
						int_u_minus_0 -= EllipticIntegral(-2, param->u_tgt, u1, -1., 1. / param->u_minus, -u0, 1., u1, -1., u2, -1.);
						int_u_minus_1 -= int_u_minus_0;
					}
				}
			}
		} else {
			const double A_1 = -1. / e, q2_sign = GSL_SIGN(q2);
			if (int root_num_U = PolySolveQuartic(d * A_1, c * A_1, 2. * A_1, A_1, &u0, &u1, &u2, &u3); root_num_U == 0) {
				const double sqrt_e = sqrt(e), sqrt_e_1 = 1. / sqrt_e;
				const double c_sqrt_e = c * sqrt_e_1, d_e = d / e;
				double h1_params[2] = {c_sqrt_e, 2. * c_sqrt_e - sqrt_e * gsl_pow_2(d_e)};
				gsl_function_fdf h1_function{
					[](double x, void *params) -> double {
						const double *param = static_cast<double *>(params);
						const double h2 = gsl_pow_2(x), h4 = gsl_pow_2(h2);
						return h2 * h4 - x * (param[0] * (h4 + 1.) - param[1] * h2) - h4 - h2 + 1;
					},
					[](double x, void *params) -> double {
						const double *param = static_cast<double *>(params);
						const double h2 = gsl_pow_2(x), h4 = gsl_pow_2(h2);
						return x * (6. * h4 - 4. * h2 - 2.) - param[0] * (5. * h4 + 1.) + param[1] * 3. * h2;
					},
					[](double x, void *params, double *f, double *df) {
						const double *param = static_cast<double *>(params);
						const double h2 = gsl_pow_2(x), h4 = gsl_pow_2(h2);
						*f = h2 * h4 - x * (param[0] * (h4 + 1.) - param[1] * h2) - h4 - h2 + 1;
						*df = x * (6. * h4 - 4. * h2 - 2.) - param[0] * (5. * h4 + 1.) + param[1] * 3. * h2;
					},
					h1_params};
				DerivativeSolver h1_solver(&h1_function, 0.5); // f(0) = 1, f(1) = -sqrt(e) * (d/e)^2 < 0
				if (int status = h1_solver.Solve(); status != GSL_SUCCESS)
					return status;
				const double h1 = h1_solver.Root(), h2 = 1. / h1, g1 = d_e / (h2 - h1);
				I_u_0 = EllipticIntegral4Complex(0, param->u_obs, param->u_tgt, 1., 0., sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				I_u_1 = 0.;
				int_u_plus_0 = EllipticIntegral4Complex(-2, param->u_obs, param->u_tgt, -1., 1. / param->u_plus, sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				int_u_minus_0 = EllipticIntegral4Complex(-2, param->u_obs, param->u_tgt, -1., 1. / param->u_minus, sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
			} else if (root_num_U == 2) {
				const double f = A_1 / (u0 * u1); // 1. / (−a^2*q^2*u_1*u_4)
				I_u_0 = EllipticIntegral2Complex(0, param->u_obs, u1, 1., 0., f, (1. / u0 + 1. / u1) * f, 1., -u0, 1., q2_sign * u1, -q2_sign);
				int_u_plus_0 = EllipticIntegral2Complex(-2, param->u_obs, u1, -1., 1. / param->u_plus, f, (1. / u0 + 1. / u1) * f, 1., -u0, 1., q2_sign * u1, -q2_sign);
				int_u_minus_0 = EllipticIntegral2Complex(-2, param->u_obs, u1, -1., 1. / param->u_minus, f, (1. / u0 + 1. / u1) * f, 1., -u0, 1., q2_sign * u1, -q2_sign);
				if (u1 >= param->u_tgt) {
					I_u_1 = I_u_0 * 2;
					I_u_0 -= EllipticIntegral2Complex(0, param->u_tgt, u1, 1., 0., f, (1. / u0 + 1. / u1) * f, 1., -u0, 1., q2_sign * u1, -q2_sign);
					I_u_1 -= I_u_0;
					int_u_plus_1 = int_u_plus_0 * 2.;
					int_u_plus_0 -= EllipticIntegral2Complex(-2, param->u_tgt, u1, -1., 1. / param->u_plus, f, (1. / u0 + 1. / u1) * f, 1., -u0, 1., q2_sign * u1, -q2_sign);
					int_u_plus_1 -= int_u_plus_0;
					int_u_minus_1 = int_u_minus_0 * 2.;
					int_u_minus_0 -= EllipticIntegral2Complex(-2, param->u_tgt, u1, -1., 1. / param->u_minus, f, (1. / u0 + 1. / u1) * f, 1., -u0, 1., q2_sign * u1, -q2_sign);
					int_u_minus_1 -= int_u_minus_0;
				}
			} else { // root_num_U == 4
				I_u_0 = EllipticIntegral(0, param->u_obs, u1, 1., 0., -u0, 1., u1, -1., u2, -1., u3, -1.);
				int_u_plus_0 = EllipticIntegral(-2, param->u_obs, u1, -1., 1. / param->u_plus, -u0, 1., u1, -1., u2, -1., u3, -1.);
				int_u_minus_0 = EllipticIntegral(-2, param->u_obs, u1, -1., 1. / param->u_minus, -u0, 1., u1, -1., u2, -1., u3, -1.);
				if (u1 >= param->u_tgt) {
					I_u_1 = I_u_0 * 2;
					I_u_0 -= EllipticIntegral(0, param->u_tgt, u1, 1., 0., -u0, 1., u1, -1., u2, -1., u3, -1.);
					I_u_1 -= I_u_0;
					int_u_plus_1 = int_u_plus_0 * 2.;
					int_u_plus_0 -= EllipticIntegral(-2, param->u_tgt, u1, -1., 1. / param->u_plus, -u0, 1., u1, -1., u2, -1., u3, -1.);
					int_u_plus_1 -= int_u_plus_0;
					int_u_minus_1 = int_u_minus_0 * 2.;
					int_u_minus_0 -= EllipticIntegral(0, param->u_tgt, u1, -1., 1. / param->u_minus, -u0, 1., u1, -1., u2, -1., u3, -1.);
					int_u_minus_1 -= int_u_minus_0;
				}
			}
		}
		// M = q^2 + (a^2-q^2-l^2)*mu^2-a^2*mu^4
		double M_plus, M_minus;
		if (int root_num_M = gsl_poly_solve_quadratic(-param->a2, c, q2, &M_minus, &M_plus); root_num_M == 0)
			return GSL_FAILURE;
		double mu_plus, mu_minus;
		double int_to_mu_plus, int_mu_full_turn = gsl_sf_ellint_Kcomp(sqrt(1. - gsl_pow_2(mu_minus / mu_plus)), GSL_PREC_DOUBLE); // FIXME: check equal
		if (M_minus > 0.) {
			if (param->mu_obs >= 0.) {
				mu_plus = sqrt(M_plus);
				mu_minus = sqrt(M_minus);
				int_to_mu_plus = EllipticIntegral(0, param->mu_obs, mu_plus, 1., 0., mu_plus, 1., mu_minus, 1., -mu_minus, 1., mu_plus, -1.) / abs(param->a);
				int_mu_full_turn = EllipticIntegral(0, mu_minus, mu_plus, 1., 0., mu_plus, 1., mu_minus, 1., -mu_minus, 1., mu_plus, -1.) / abs(param->a);
			} else {
				mu_plus = -sqrt(M_minus);
				mu_minus = -sqrt(M_plus);
				int_to_mu_plus = EllipticIntegral(0, param->mu_obs, mu_plus, 1., 0., -mu_minus, 1., mu_plus, -1., -mu_plus, -1., -mu_minus, -1.) / abs(param->a);
				int_mu_full_turn = EllipticIntegral(0, mu_minus, mu_plus, 1., 0., -mu_minus, 1., mu_plus, -1., -mu_plus, -1., -mu_minus, -1.) / abs(param->a);
			}
		} else {
			mu_plus = sqrt(M_plus);
			mu_minus = -mu_plus;
			int_to_mu_plus = EllipticIntegral2Complex(0, param->mu_obs, mu_plus, 1., 0., 1., 0., -M_minus, -mu_minus, 1., mu_plus, -1.) / abs(param->a);
			int_mu_full_turn = EllipticIntegral2Complex(0, mu_minus, mu_plus, 1., 0., 1., 0., -M_minus, -mu_minus, 1., mu_plus, -1.) / abs(param->a);
		}
		const int alpha_1 = beta >= 0. ? 1. : -1.;
		const double n = M_minus > 0. ? (M_plus + M_minus) / (1 - M_plus) : 1. / (1. + gsl_pow_2(mu_plus)) - 1.;
		double A;
		double x;
		double k;
		if (M_minus > 0.) {
			A = mu_plus * abs(param->a);
			x = sqrt((M_plus - gsl_pow_2(param->mu_obs)) / (M_plus - M_minus));
			k = sqrt(1. - gsl_pow_2(mu_minus / mu_plus));
		} else {
			A = sqrt(M_plus - M_minus) * abs(param->a);
			x = sqrt(1. - gsl_pow_2(param->mu_obs / mu_plus));
			k = sqrt(1. + M_minus / (gsl_pow_2(mu_plus) - M_minus));
		}
		if (int status = MuFSolver(alpha_1, M_plus, M_minus, mu_plus, mu_minus, q2, I_u_0, int_to_mu_plus, int_mu_full_turn, A, x, k, param->a, param->mu_obs, mu_f_0); status != GSL_SUCCESS)
			return status;
		const double phi_mu_0 = l * (-I_u_0 + gsl_sf_ellint_P(x, k, n, GSL_PREC_DOUBLE) / (A * (1. - M_plus)));
		const double phi_u_0 = param->u_r * ((l / param->u_plus + 2. * (param->a - l)) * int_u_plus_0 - (l / param->u_minus + 2. * (param->a - l)) * int_u_minus_0);
		gsl_vector_set(delta_u_mu_phi, 1, mu_f_0 - param->mu_tgt);
		gsl_vector_set(delta_u_mu_phi, 2, phi_u_0 + phi_mu_0 - param->phi_tgt);
		if (I_u_1 == 0.) {
			gsl_vector_set(delta_u_mu_phi, 0, 0.);
			return GSL_SUCCESS;
		}
		gsl_vector_set(delta_u_mu_phi, 0, param->u_tgt - u1);
		if (isnan(I_u_1))
			return GSL_SUCCESS;
		if (int status = MuFSolver(alpha_1, M_plus, M_minus, mu_plus, mu_minus, q2, I_u_1, int_to_mu_plus, int_mu_full_turn, A, x, k, param->a, param->mu_obs, mu_f_1); status != GSL_SUCCESS)
			return status;
		const double phi_mu_1 = l * (-I_u_1 + gsl_sf_ellint_P(x, k, n, GSL_PREC_DOUBLE) / (A * (1. - M_plus)));
		const double phi_u_1 = param->u_r * ((l / param->u_plus + 2. * (param->a - l)) * int_u_plus_1 - (l / param->u_minus + 2. * (param->a - l)) * int_u_minus_1);
		const double square_diff_0 = gsl_pow_2(acos(mu_f_0) - param->theta_tgt) + gsl_pow_2(phi_u_0 + phi_mu_0 - param->phi_tgt);
		const double square_diff_1 = gsl_pow_2(acos(mu_f_1) - param->theta_tgt) + gsl_pow_2(phi_u_1 + phi_mu_1 - param->phi_tgt);
		if (square_diff_0 > square_diff_1) {
			gsl_vector_set(delta_u_mu_phi, 1, mu_f_1 - param->mu_tgt);
			gsl_vector_set(delta_u_mu_phi, 2, phi_u_1 + phi_mu_1 - param->phi_tgt);
		}
		return GSL_SUCCESS;
	}
	int Kerr::FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_target, const double theta_target, const double phi_target, double &alpha, double &beta, double photon[]) { // TODO:
		const double sin_theta_target = abs(sin(theta_target)), cos_theta_target = GSL_SIGN(theta_target) * cos(theta_target), sin_phi_target = sin(phi_target), cos_phi_target = cos(phi_target);
		const double sin2_theta_target = gsl_pow_2(sin_theta_target), cos2_theta_target = gsl_pow_2(cos_theta_target);
		const double sin2_theta_observer = gsl_pow_2(sin_theta_observer), cos2_theta_observer = gsl_pow_2(cos_theta_observer);
		const double cos_observer_target = sin_theta_observer * sin_theta_target * cos_phi_target + cos_theta_observer * cos_theta_target;
		const double delta_phi = acos(cos_observer_target), sin_observer_target = sqrt(1. - gsl_pow_2(cos_observer_target));
		const double r2_target = gsl_pow_2(r_target), r2_observer = gsl_pow_2(r_observer), u_observer = 1. / r_observer, u_target = 1. / r_target;
		const double Delta_target = r_target * (r_target - 2.) + a2_, rho2_target = r2_target + a2_ * cos2_theta_target;
		const double Delta_observer = r_observer * (r_observer - 2.) + a2_, rho2_observer = r2_observer + a2_ * cos2_theta_observer;
		const double rho_2_target = 1. / rho2_target, rho_2_observer = 1. / rho2_observer;
		const double g00_target = -(1. - 2. * r_target * rho_2_target), g03_target = -2. * r_target * rho_2_target * a_ * sin2_theta_target;
		const double g11_target = rho2_target / Delta_target, g22_target = rho2_target, g33_target = (r2_target + a2_) * sin2_theta_target - g03_target * a_ * sin2_theta_target;
		const double g00_observer = -(1. - 2. * r_observer * rho_2_observer), g03_observer = -2. * r_observer * rho_2_observer * a_ * sin2_theta_observer;
		const double g11_observer = rho2_observer / Delta_observer, g22_observer = rho2_observer, g33_observer = (r2_observer + a2_) * sin2_theta_observer - g03_observer * a_ * sin2_theta_observer;
		// The observer's frame
		const double e0_observer = sqrt(-1. / g00_observer), e1_observer = sqrt(1. / g11_observer), e2_observer = sqrt(1. / g22_observer), e3_observer = sqrt(1. / g33_observer);
		KerrFastTraceParameters fast_trace_parameters(a_, u_observer, u_target, cos_theta_observer, cos_theta_target, theta_target, phi_target);
		GslBlock collector;
		gsl_vector *alpha_beta_initial_value = collector.VectorAlloc(3);
		gsl_multiroot_function alpha_beta_function{DifferenceOfUMuPhi, 2, &fast_trace_parameters};
		MultiFunctionSolver alpha_beta_solver(&alpha_beta_function, alpha_beta_initial_value, 3);
		if (int status = alpha_beta_solver.Solve(); status != GSL_SUCCESS)
			return status;
		// const double R = sqrt(p[7] + gsl_pow_2(gsl_vector_get(x, 0)) + gsl_pow_2(gsl_vector_get(x, 1))), R_1 = 1. / R;
		// const double E = -p[8] * p[12] + p[9] * p[14] * gsl_vector_get(x, 0) * R_1, L = (p[1] == 0. || p[1] == M_PI) ? 0. : p[9] * p[12] - p[11] * p[14] * gsl_vector_get(x, 0) * R_1;
		alpha = gsl_vector_get(alpha_beta_solver.Root(), 0);
		beta = gsl_vector_get(alpha_beta_solver.Root(), 1);
		// [1., 1., beta / r2, -alpha / r2 / sin_theta] * r2 / (r2 + alpha2 + beta2)
		const double E = 1. - 2. * u_observer, L = -alpha * sin_theta_observer;
		photon[0] = 0.;
		photon[1] = r_target;
		photon[2] = theta_target;
		photon[3] = phi_target;
		photon[4] = (gsl_pow_2(g03_target * g03_target) - g00_target * g33_target) / (E * g33_target + L * g03_target);
		photon[7] = -(g00_target * L + g03_target * E) / (g03_target * L + g33_target * E);
		return GSL_FAILURE;
	}
	double Kerr::Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			if (time == T)
				return (1. - 2. * y[1] / (gsl_pow_2(y[1]) + a2_ * gsl_pow_2(cos(y[2]))) * (1. - a_ * gsl_pow_2(sin(y[2])) * y[7])) / y[4];
			// time == TAU
			return y[4] - 2. * y[1] / (gsl_pow_2(y[1]) + a2_ * gsl_pow_2(cos(y[2]))) * (y[4] - a_ * gsl_pow_2(sin(y[2])) * y[7]);
		} // dynamics == HAMILTONIAN
		return 1. - y[4];
	}
	double Kerr::AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			const double r2 = gsl_pow_2(y[1]), sin2_theta = gsl_pow_2(sin(y[2]));
			const double r_rho_2 = 2. * y[1] / (r2 + a2_ * gsl_pow_2(cos(y[2])));
			if (time == T)
				return (-r_rho_2 * a_ + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta / y[4];
			// time == TAU
			return (-r_rho_2 * a_ * y[4] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta;
		} // dynamics == HAMILTONIAN
		return y[7];
	}
	double Kerr::CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			const double r2 = gsl_pow_2(y[1]), sin2_theta = gsl_pow_2(sin(y[2])), cos2_theta = gsl_pow_2(cos(y[2]));
			const double rho2 = r2 + a2_ * cos2_theta;
			const double r_rho_2 = 2. * y[1] / rho2;
			if (time == T)
				return mu2 * cos2_theta * a2_ + (gsl_pow_2(rho2 * y[6]) + cos2_theta * (gsl_pow_2(-r_rho_2 * a_ + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta - a2_ * gsl_pow_2(r_rho_2 * (1. - a_ * sin2_theta * y[7]) - 1.))) / gsl_pow_2(y[4]);
			// time == TAU
			return mu2 * cos2_theta * a2_ + gsl_pow_2(rho2 * y[6]) + cos2_theta * (gsl_pow_2(-r_rho_2 * a_ * y[4] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta - a2_ * gsl_pow_2(r_rho_2 * (y[4] - a_ * sin2_theta * y[7]) - y[4]));
		} // dynamics == HAMILTONIAN
		return gsl_pow_2(y[6]) + gsl_pow_2(cos(y[2])) * (a2_ * (mu2 - gsl_pow_2(1. - y[4])) + gsl_pow_2(y[7] / sin(y[2])));
	}
	int Kerr::NormalizeTimelikeGeodesic(double y[]) {
		const double r = y[1], r2 = gsl_pow_2(r), a2_r2 = a2_ + r2;
		const double sin2_theta = gsl_pow_2(sin(y[2])), sin4_theta = gsl_pow_2(sin2_theta);
		const double rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double r_rho_2 = 2. * r / rho2;
		// y[7] += 2. * a_ * r / (gsl_pow_2(a2_r2) - a2_ * Delta * sin2_theta);
		y[4] = sqrt(1. - r_rho_2 + 2. * r_rho_2 * a_ * sin2_theta * y[7] - (rho2 / (r2 - 2. * r + a2_) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + (a2_r2 * sin2_theta + r_rho_2 * a2_ * sin4_theta) * gsl_pow_2(y[7])));
		return isnan(y[4]) ? GSL_EDOM : GSL_SUCCESS;
	}
	int Kerr::NormalizeNullGeodesic(double y[], double frequency) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sin2_theta = gsl_pow_2(sin(y[2])), sin4_theta = gsl_pow_2(sin2_theta);
		const double rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double r_rho_2 = 2. * r / rho2;
		const double a = rho2 / (r2 - 2. * r + a2_) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2_ + r2) * sin2_theta + r_rho_2 * a2_ * sin4_theta) * gsl_pow_2(y[7]);
		const double b = -2. * r_rho_2 * a_ * sin2_theta * y[7];
		const double c = r_rho_2 - 1.;
		const double coefficient = GSL_SIGN(frequency) * 0.5 * (-b + sqrt(gsl_pow_2(b) - 4. * a * c)) / a;
		if (isnan(coefficient))
			return GSL_EDOM;
		y[4] = frequency;
		y[5] *= coefficient;
		y[6] *= coefficient;
		y[7] *= coefficient;
		return GSL_SUCCESS;
	}
	unique_ptr<Integrator> Kerr::GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion) {
		if (time == T) {
			if (dynamics == LAGRANGIAN) {
				if (motion == GEODESIC)
					return make_unique<Integrator>(KerrTLagrangianGeodesic, Jacobian, this);
				else if (motion == CIRCULAR)
					return make_unique<Integrator>(KerrTLagrangianCircular, Jacobian, this);
				else if (motion == HELICAL)
					return make_unique<Integrator>(KerrTLagrangianHelical, Jacobian, this);
			} else if (dynamics == HAMILTONIAN && motion == GEODESIC)
				return make_unique<Integrator>(KerrTHamiltonianGeodesic, Jacobian, this);
		} else if (time == TAU) {
			if (dynamics == LAGRANGIAN)
				return make_unique<Integrator>(KerrTauLagrangianGeodesic, Jacobian, this);
			else if (dynamics == HAMILTONIAN)
				return make_unique<Integrator>(KerrTauHamiltonianGeodesic, Jacobian, this);
		}
		PrintlnError("{}::GetIntegrator() {}, {}, {} invaild", Name(), time, dynamics, motion);
		return nullptr;
	}

	KerrNewman::KerrNewman(double spin, double charge) : Kerr(spin), r_Q_(0.5 * sqrt(M_1_PI) * charge), r_Q2_(r_Q_ * r_Q_), r_Q4_(r_Q2_ * r_Q2_) {}
	string KerrNewman::Name() {
		return "Kerr-Newman";
	}
	int KerrNewman::MetricTensor(const double position[], gsl_matrix *metric) {
		const double r = position[1], r2 = gsl_pow_2(r), Delta = r2 - 2. * r + a2_ + r_Q2_, rho2 = r2 + a2_ * gsl_pow_2(cos(position[2])), sin2_theta = gsl_pow_2(sin(position[2]));
		if (rho2 == 0. || Delta == 0.)
			return GSL_EZERODIV;
		const double r_rho_2 = (2. * r - r_Q2_) / rho2;
		gsl_matrix_set_zero(metric);
		gsl_matrix_set(metric, 0, 0, -(1. - r_rho_2));
		gsl_matrix_set(metric, 0, 3, -r_rho_2 * a_ * sin2_theta);
		gsl_matrix_set(metric, 1, 1, rho2 / Delta);
		gsl_matrix_set(metric, 2, 2, rho2);
		gsl_matrix_set(metric, 3, 0, -r_rho_2 * a_ * sin2_theta);
		gsl_matrix_set(metric, 3, 3, (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * sin2_theta);
		return GSL_SUCCESS;
	}
	double KerrNewman::DotProduct(const double position[], const double x[], const double y[], const size_t dimension) {
		const double r = position[1], r2 = gsl_pow_2(r), Delta = r2 - 2. * r + a2_ + r_Q2_, rho2 = r2 + a2_ * gsl_pow_2(cos(position[2])), r_rho_2 = (2. * r - r_Q2_) / rho2, sin2_theta = gsl_pow_2(sin(position[2]));
		if (dimension == 3)
			return (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * gsl_pow_2(sin2_theta)) * x[3] * y[3];
		return (r_rho_2 - 1.) * x[0] * y[0] - r_rho_2 * a_ * sin2_theta * (x[0] * y[3] + x[3] * y[0]) + (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * gsl_pow_2(sin2_theta)) * x[3] * y[3];
	}
	double KerrNewman::DistanceSquare(const double x[], const double y[], const size_t dimension) {
		const double r = x[1], r2 = gsl_pow_2(r), d0 = x[0] - y[0], d3 = PhiDifference(x[3] - y[3]);
		const double Delta = r2 - 2. * r + a2_ + r_Q2_, rho2 = r2 + a2_ * gsl_pow_2(cos(x[2])), r_rho_2 = (2. * r - r_Q2_) / rho2, sin2_theta = gsl_pow_2(sin(x[2]));
		if (dimension == 3)
			return (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + ((r2 + a2_) * sin2_theta + 2. * r * a2_ * gsl_pow_2(sin2_theta) / rho2) * gsl_pow_2(d3);
		return (r_rho_2 - 1.) * gsl_pow_2(d0) - 2. * r_rho_2 * a_ * sin2_theta * d0 * d3 + (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * gsl_pow_2(sin2_theta)) * gsl_pow_2(d3);
	}
	int KerrNewman::LagrangianToHamiltonian(double y[]) {
		const double dt_dtau = 1. / y[4], r2 = gsl_pow_2(y[1]), sin2_theta = gsl_pow_2(sin(y[2]));
		const double Delta = r2 - 2. * y[1] + a2_ + r_Q2_, rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		if (rho2 == 0. || Delta == 0.)
			return GSL_EZERODIV;
		const double r_rho_2 = (2. * y[1] - r_Q2_) / rho2;
		y[4] = (y[4] - 1. + r_rho_2 * (1. - a_ * sin2_theta * y[7])) * dt_dtau; // 1 + p_t
		y[5] *= rho2 / Delta * dt_dtau;
		y[6] *= rho2 * dt_dtau;
		y[7] = (-r_rho_2 * a_ + (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * y[7]) * sin2_theta * dt_dtau;
		return GSL_SUCCESS;
	}
	int KerrNewman::HamiltonianToLagrangian(double y[]) {
		const double pt = y[4] - 1., r2 = gsl_pow_2(y[1]), sin2_theta = gsl_pow_2(sin(y[2]));
		const double Delta = r2 - 2. * y[1] + a2_ + r_Q2_, rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		if (rho2 == 0. || Delta == 0.)
			return GSL_EZERODIV;
		const double rho_2 = 1. / rho2, r_rho_2 = (2. * y[1] - r_Q2_) / rho2;
		y[4] = -Delta / ((Delta + r_rho_2 * (a2_ + r2)) * pt + r_rho_2 * a_ * y[7]);
		y[5] *= Delta * y[4] * rho_2;
		y[6] *= y[4] * rho_2;
		if (y[7] == 0.)
			y[7] = (-r_rho_2 * a_ * pt + (1. - r_rho_2) / sin2_theta * y[7]) / Delta * y[4];
		else if (sin2_theta == 0.)
			return GSL_EZERODIV;
		else
			y[7] = (-r_rho_2 * a_ * pt + (1. - r_rho_2) / sin2_theta * y[7]) / Delta * y[4];
		return GSL_SUCCESS;
	}
	int KerrNewman::FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_target, const double theta_target, const double phi_target, double &alpha, double &beta, double photon[]) {
		return GSL_FAILURE;
	}
	double KerrNewman::Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			if (time == T)
				return (1. - (2. * y[1] - r_Q2_) / (gsl_pow_2(y[1]) + a2_ * gsl_pow_2(cos(y[2]))) * (1. - a_ * gsl_pow_2(sin(y[2])) * y[7])) / y[4];
			// time == TAU
			return y[4] - (2. * y[1] - r_Q2_) / (gsl_pow_2(y[1]) + a2_ * gsl_pow_2(cos(y[2]))) * (y[4] - a_ * gsl_pow_2(sin(y[2])) * y[7]);
		} // dynamics == HAMILTONIAN
		return 1. - y[4];
	}
	double KerrNewman::AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			const double r2 = gsl_pow_2(y[1]), sin2_theta = gsl_pow_2(sin(y[2]));
			const double r_rho_2 = (2. * y[1] - r_Q2_) / (r2 + a2_ * gsl_pow_2(cos(y[2])));
			if (time == T)
				return (-r_rho_2 * a_ + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta / y[4];
			// time == TAU
			return (-r_rho_2 * a_ * y[4] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta;
		} // dynamics == HAMILTONIAN
		return y[7];
	}
	double KerrNewman::CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) {
		return GSL_NAN;
	}
	int KerrNewman::NormalizeTimelikeGeodesic(double y[]) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sin2_theta = gsl_pow_2(sin(y[2]));
		const double rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double r_rho_2 = (2. * r - r_Q2_) / rho2;
		y[4] = sqrt(1. - r_rho_2 + 2. * r_rho_2 * a_ * sin2_theta * y[7] - (rho2 / (r2 - 2. * r + a2_ + r_Q2_) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * sin2_theta) * gsl_pow_2(y[7])));
		return isnan(y[4]) ? GSL_EDOM : GSL_SUCCESS;
	}
	int KerrNewman::NormalizeNullGeodesic(double y[], double frequency) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sin2_theta = gsl_pow_2(sin(y[2]));
		const double rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double r_rho_2 = (2. * r - r_Q2_) / rho2;
		const double a = rho2 / (r2 - 2. * r + a2_ + r_Q2_) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * sin2_theta) * gsl_pow_2(y[7]);
		const double b = -2. * r_rho_2 * a_ * sin2_theta * y[7];
		const double c = r_rho_2 - 1.;
		const double coefficient = GSL_SIGN(frequency) * 0.5 * (-b + sqrt(gsl_pow_2(b) - 4. * a * c)) / a;
		if (isnan(coefficient))
			return GSL_EDOM;
		y[4] = frequency;
		y[5] *= coefficient;
		y[6] *= coefficient;
		y[7] *= coefficient;
		return GSL_SUCCESS;
	}
	unique_ptr<Integrator> KerrNewman::GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion) {
		return nullptr;
	}

	KerrTaubNUT::KerrTaubNUT(double spin, double NUT) : Kerr(spin), l_(NUT), l2_(l_ * l_), l4_(l2_ * l2_) {}
	string KerrTaubNUT::Name() {
		return "Kerr-Taub-NUT";
	}
	int KerrTaubNUT::MetricTensor(const double position[], gsl_matrix *metric) {
		return GSL_FAILURE;
	}
	double KerrTaubNUT::DotProduct(const double position[], const double x[], const double y[], const size_t dimension) {
		const double r = position[1], sin_theta = sin(position[2]), cos_theta = cos(position[2]);
		const double r2 = gsl_pow_2(r), sin2_theta = gsl_pow_2(sin_theta);
		const double Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + gsl_pow_2(l_ + a_ * cos_theta), chi = a_ * sin2_theta - 2. * l_ * cos_theta;
		const double rho_2 = 1. / rho2;
		if (dimension == 3)
			return (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + (gsl_pow_2(rho2 + a_ * chi) * sin2_theta - chi * chi * Delta) * rho_2 * x[3] * y[3];
		return (a2_ * sin2_theta - Delta) * rho_2 * x[0] * y[0] - 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * rho_2 * (x[0] * y[3] + x[3] * y[0]) + (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + (gsl_pow_2(rho2 + a_ * chi) * sin2_theta - chi * chi * Delta) * rho_2 * x[3] * y[3];
	}
	double KerrTaubNUT::DistanceSquare(const double x[], const double y[], const size_t dimension) {
		const double r = x[1], sin_theta = sin(x[2]), cos_theta = cos(x[2]), d0 = x[0] - y[0], d3 = PhiDifference(x[3] - y[3]);
		const double r2 = gsl_pow_2(r), sin2_theta = gsl_pow_2(sin_theta);
		const double Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + gsl_pow_2(l_ + a_ * cos_theta), chi = a_ * sin2_theta - 2. * l_ * cos_theta;
		const double rho_2 = 1. / rho2;
		if (dimension == 3)
			return (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + (gsl_pow_2(rho2 + a_ * chi) * sin2_theta - gsl_pow_2(chi) * Delta) * rho_2 * gsl_pow_2(d3);
		return (a2_ * sin2_theta - Delta) * rho_2 * gsl_pow_2(d0) - 4. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * rho_2 * d0 * d3 + (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + (gsl_pow_2(rho2 + a_ * chi) * sin2_theta - gsl_pow_2(chi) * Delta) * rho_2 * gsl_pow_2(d3);
	}
	int KerrTaubNUT::LagrangianToHamiltonian(double y[]) {
		const double dt_dtau = 1. / y[4], r = y[1], r2 = gsl_pow_2(r);
		const double sin2_theta = gsl_pow_2(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + gsl_pow_2(l_ + a_ * cos_theta);
		if (rho2 == 0. || Delta == 0.)
			return GSL_EZERODIV;
		const double rho_2 = 1. / rho2;
		y[4] = (y[4] * rho2 - Delta + a2_ * sin2_theta - 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]) * rho_2 * dt_dtau; // 1 + p_t
		y[5] *= rho2 / Delta * dt_dtau;
		y[6] *= rho2 * dt_dtau;
		y[7] = (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) + (gsl_pow_2(r2 + l2_ + a2_) * sin2_theta - gsl_pow_2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * y[7]) * rho_2 * dt_dtau;
		return GSL_SUCCESS;
	}
	int KerrTaubNUT::HamiltonianToLagrangian(double y[]) {
		const double pt = y[4] - 1., r = y[1], r2 = gsl_pow_2(r);
		const double sin2_theta = gsl_pow_2(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + gsl_pow_2(l_ + a_ * cos_theta);
		if (rho2 == 0. || Delta == 0.)
			return GSL_EZERODIV;
		const double rho_2 = 1. / rho2;
		y[4] = Delta * rho2 * sin2_theta / ((Delta * gsl_pow_2(a_ * sin2_theta - 2. * l_ * cos_theta) - gsl_pow_2(r2 + l2_ + a2_) * sin2_theta) * pt - 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]);
		y[5] *= Delta * rho_2 * y[4];
		y[6] *= rho_2 * y[4];
		y[7] = (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * pt + (Delta - a2_ * sin2_theta) * y[7]) / (Delta * rho2 * sin2_theta) * y[4];
		return GSL_SUCCESS;
	}
	int KerrTaubNUT::FastTrace(const double r_observer, const double theta_observer, const double sin_theta_observer, const double cos_theta_observer, const double r_target, const double theta_target, const double phi_target, double &alpha, double &beta, double photon[]) {
		return GSL_FAILURE;
	}
	double KerrTaubNUT::Energy(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			const double r = y[1], r2 = gsl_pow_2(r);
			const double sin2_theta = gsl_pow_2(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
			const double Delta = r2 - 2. * r - l2_ + a2_;
			if (time == T)
				return (Delta - a2_ * sin2_theta + 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]) / ((r2 + gsl_pow_2(l_ + a_ * cos_theta)) * y[4]);
			// time == TAU
			return ((Delta - a2_ * sin2_theta) * y[4] + 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]) / (r2 + gsl_pow_2(l_ + a_ * cos_theta));
		} // dynamics == HAMILTONIAN
		return 1. - y[4];
	}
	double KerrTaubNUT::AngularMomentum(const double y[], TimeSystem time, DynamicalSystem dynamics) {
		if (dynamics == LAGRANGIAN) {
			const double r = y[1], r2 = gsl_pow_2(r);
			const double sin2_theta = gsl_pow_2(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
			const double Delta = r2 - 2. * r - l2_ + a2_;
			if (time == T)
				return (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) + (gsl_pow_2(r2 + l2_ + a2_) * sin2_theta - gsl_pow_2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * y[7]) / ((r2 + gsl_pow_2(l_ + a_ * cos_theta)) * y[4]);
			// time == TAU
			return (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[4] + (gsl_pow_2(r2 + l2_ + a2_) * sin2_theta - gsl_pow_2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * y[7]) / (r2 + gsl_pow_2(l_ + a_ * cos_theta));
		} // dynamics == HAMILTONIAN
		return y[7];
	}
	double KerrTaubNUT::CarterConstant(const double y[], const double mu2, TimeSystem time, DynamicalSystem dynamics) {
		return GSL_NAN;
	}
	int KerrTaubNUT::NormalizeTimelikeGeodesic(double y[]) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sin2_theta = gsl_pow_2(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		const double rho2 = r2 + gsl_pow_2(l_ + a_ * cos_theta);
		// y[7] += 2. * a_ * r / (gsl_pow_2(a2_ + r2) - a2_ * Delta * sin2_theta);
		y[4] = sqrt(((Delta - a2_ * sin2_theta) + 4. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7] - (gsl_pow_2(r2 + l2_ + a2_) * sin2_theta - gsl_pow_2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * gsl_pow_2(y[7])) / rho2 - rho2 * (gsl_pow_2(y[5]) / Delta + gsl_pow_2(y[6])));
		return isnan(y[4]) ? GSL_EDOM : GSL_SUCCESS;
	}
	int KerrTaubNUT::NormalizeNullGeodesic(double y[], double frequency) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sin2_theta = gsl_pow_2(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		const double l_a_cos_theta = l_ + a_ * cos_theta;
		const double rho2 = r2 + gsl_pow_2(l_a_cos_theta), rho_2 = 1. / rho2;
		const double chi = a_ * sin2_theta - 2. * l_ * cos_theta, rho2_a_chi = r2 + l2_ + a2_;
		const double a = rho2 / Delta * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + rho_2 * (gsl_pow_2(rho2_a_chi) * sin2_theta - gsl_pow_2(chi) * Delta) * gsl_pow_2(y[7]);
		const double b = -4. * rho_2 * ((r + l2_) * chi + l_ * cos_theta * rho2_a_chi) * y[7];
		const double c = -rho_2 * (Delta - a2_ * sin2_theta);
		const double coefficient = GSL_SIGN(frequency) * 0.5 * (-b + sqrt(gsl_pow_2(b) - 4. * a * c)) / a;
		if (isnan(coefficient))
			return GSL_EDOM;
		y[4] = frequency;
		y[5] *= coefficient;
		y[6] *= coefficient;
		y[7] *= coefficient;
		return GSL_SUCCESS;
	}
	unique_ptr<Integrator> KerrTaubNUT::GetIntegrator(TimeSystem time, DynamicalSystem dynamics, MotionMode motion) {
		if (motion != GEODESIC)
			PrintlnError("Motion not supported!");
		if (time == T) {
			if (dynamics == LAGRANGIAN)
				return make_unique<Integrator>(KerrTaubNutTLagrangianGeodesic, Jacobian, this);
			else if (dynamics == HAMILTONIAN)
				return make_unique<Integrator>(KerrTaubNutTHamiltonianGeodesic, Jacobian, this);
		} else if (time == TAU) {
			return make_unique<Integrator>(KerrTaubNutTauLagrangianGeodesic, Jacobian, this);
		}
		PrintlnError("{}::GetIntegrator() {}, {}, {} invaild", Name(), time, dynamics, motion);
		return nullptr;
	}
	int NewtonTLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		const int PN = *static_cast<int *>(params);
		const double r_1 = 1. / y[1], r_2 = gsl_pow_2(r_1), sin_theta = abs(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
		const double v_tan2 = gsl_pow_2(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(sin_theta * y[7]));
		const double rdot = y[5], rdot2 = gsl_pow_2(y[5]);
		const double v2 = rdot2 + v_tan2;
		dydt[0] = y[4];
		dydt[1] = y[5];
		dydt[2] = y[6];
		dydt[3] = y[7];
		dydt[4] = 0.;
		double A = 1., B = 0.;
		if (PN & 1) {
			A += v2 - 4. * r_1;
			B += -4. * rdot;
		}
		if (PN & 2) {
			A += r_1 * (-2. * rdot2 + 9. * r_1);
			B += 2. * r_1 * rdot;
		}
		if (PN & 8) {
			A += -16. * gsl_pow_3(r_1) + rdot2 * r_2;
			B += -4. * rdot * r_2;
		}
		const double B_r2 = B * r_2;
		dydt[5] = r_1 * v_tan2 - r_2 * A - B_r2 * y[5];
		dydt[6] = sin_theta * cos_theta * gsl_pow_2(y[7]) - (2. * y[5] * r_1 + B_r2) * y[6];
		dydt[7] = -(2. * cos_theta / sin_theta * y[6] + 2. * y[5] * r_1 + B_r2) * y[7];
		return GSL_SUCCESS;
	}
	int PN1TLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		const double PN1 = *static_cast<double *>(params);
		const double r_1 = 1. / y[1], r_2 = gsl_pow_2(r_1), sin_theta = abs(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
		const double v_tan2 = gsl_pow_2(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(sin_theta * y[7]));
		const double v2 = gsl_pow_2(y[5]) + v_tan2;
		dydt[0] = y[4];
		dydt[1] = y[5];
		dydt[2] = y[6];
		dydt[3] = y[7];
		dydt[4] = 0.;
		const double A = 1. + (v2 - 4. * r_1) * PN1, B_r2 = -4. * y[5] * PN1 * r_2;
		dydt[5] = r_1 * v_tan2 - r_2 * A - B_r2 * y[5];
		dydt[6] = sin_theta * cos_theta * gsl_pow_2(y[7]) - (2. * y[5] * r_1 + B_r2) * y[6];
		dydt[7] = -(2. * cos_theta / sin_theta * y[6] + 2. * y[5] * r_1 + B_r2) * y[7];
		return GSL_SUCCESS;
	}
	int SchwarzschildTLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = y[5]; // dr/dt
		dydt[2] = y[6]; // d\theta/dt
		dydt[3] = y[7]; // d\phi/dt
		const double r = y[1], sin_theta = abs(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
		const double rm2 = r - 2., rm3 = r - 3.;
		const double rm2r_1 = 1. / (rm2 * r);
		// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = 2. * y[5] * rm2r_1 * y[4];
		// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = -rm2 / gsl_pow_3(r) + 3. * rm2r_1 * gsl_pow_2(y[5]) + rm2 * (gsl_pow_2(y[6]) + gsl_pow_2(sin_theta * y[7]));
		// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = -2. * rm3 * rm2r_1 * y[5] * y[6] + sin_theta * cos_theta * gsl_pow_2(y[7]);
		// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		if (y[7] == 0.)
			dydt[7] = 0.;
		else if (sin_theta == 0.)
			return GSL_EZERODIV;
		else
			dydt[7] = -2. * (rm3 * rm2r_1 * y[5] + cos_theta / sin_theta * y[6]) * y[7];
		return GSL_SUCCESS;
	}
	int SchwarzschildTLagrangianCircular(double t, const double y[], double dydt[], void *params) {
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = 0.;	// dr/dt
		dydt[2] = 0.;	// d\theta/dt
		dydt[3] = y[7]; // d\phi/dt
		dydt[4] = 0.;
		dydt[5] = 0.;
		dydt[6] = 0.;
		dydt[7] = 0.;
		return GSL_SUCCESS;
	}
	int SchwarzschildTLagrangianHelical(double t, const double y[], double dydt[], void *params) {
		const double r = y[1];
		if (r <= 2.)
			return GSL_FAILURE;
		const double sin2_theta = gsl_pow_2(sin(y[2]));
		const double g00 = -1. + 2. / r, g11 = -1. / g00, g33 = gsl_pow_2(r) * sin2_theta;
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = y[5]; // dr/dt
		dydt[2] = 0.;	// d\theta/dt = 0.
		dydt[3] = y[7]; // d\phi/dt
		// d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt
		dydt[4] = sin2_theta * y[5] * (r * (-g00 - g11 * gsl_pow_2(y[5]) - gsl_pow_2(y[4])) + 1. + gsl_pow_2(g11) * gsl_pow_2(y[5])) / (y[4] * (g33 + gsl_pow_2(g33 * y[7] / y[4])));
		// d^2r/dt^2 = 0.
		dydt[5] = 0.;
		// d^2\theta/dt^2 = 0.
		dydt[6] = 0.;
		// d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
		dydt[7] = y[7] * (-2. / r * y[5] + dydt[4] / y[4]);
		return GSL_SUCCESS;
	}
	int SchwarzschildTHamiltonianGeodesic(double t, const double y[], double dydt[], void *params) {
		const double r_1 = 1. / y[1], r_2 = gsl_pow_2(r_1), g11_1 = 1. - 2. * r_1, E = 1. - y[4], L2 = gsl_pow_2(y[7]);
		const double sin_1_theta = 1. / abs(sin(y[2])), sin_2_theta = gsl_pow_2(sin_1_theta);
		//[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
		dydt[0] = g11_1 / E;						  // d\tau/dt
		dydt[1] = g11_1 * y[5] * dydt[0];			  // dr/dt
		dydt[2] = y[6] * r_2 * dydt[0];				  // d\theta/dt
		dydt[3] = y[7] * sin_2_theta * r_2 * dydt[0]; // d\phi/dt
		dydt[4] = 0.;
		dydt[5] = (-(gsl_pow_2(y[5]) + gsl_pow_2(E) / gsl_pow_2(g11_1)) + (gsl_pow_2(y[6]) + L2 * sin_2_theta) * r_1) * r_2 * dydt[0];
		dydt[6] = sin_2_theta * L2 * GSL_SIGN(y[2]) * cos(y[2]) * sin_1_theta * r_2 * dydt[0];
		dydt[7] = 0.;
		return GSL_SUCCESS;
	}
	int SchwarzschildTauLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		dydt[0] = y[4]; // dt/d\tau
		dydt[1] = y[5]; // dr/d\tau
		dydt[2] = y[6]; // d\theta/d\tau
		dydt[3] = y[7]; // d\phi/d\tau
		const double r = y[1], sin_theta = abs(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
		const double rm2 = r - 2., r_1 = 1. / r;
		const double rm2r_1 = 1. / (rm2 * r);
		// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = -2. * y[5] * rm2r_1 * y[4];
		// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = -rm2 * gsl_pow_3(r_1) * gsl_pow_2(y[4]) + rm2r_1 * gsl_pow_2(y[5]) + rm2 * (gsl_pow_2(y[6]) + gsl_pow_2(sin_theta * y[7]));
		// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = -2. * r_1 * y[5] * y[6] + sin_theta * cos_theta * gsl_pow_2(y[7]);
		// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = -2. * (r_1 * y[5] + cos_theta / sin_theta * y[6]) * y[7];
		return GSL_SUCCESS;
	}
	int KerrTLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = static_cast<Kerr *>(params);
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = y[5]; // dr/dt
		dydt[2] = y[6]; // d\theta/dt
		dydt[3] = y[7]; // d\phi/dt
		const double r = y[1], r2 = gsl_pow_2(r), r4 = gsl_pow_4(r), a = kerr->a_, a2 = kerr->a2_, a4 = kerr->a4_, a2_r2 = a2 + r2;
		const double sin_theta = abs(sin(y[2])), sin2_theta = gsl_pow_2(sin_theta), sin4_theta = gsl_pow_2(sin2_theta), cos_theta = GSL_SIGN(y[2]) * cos(y[2]), cos2_theta = gsl_pow_2(cos_theta), sin_theta_cos_theta = sin_theta * cos_theta, cot_theta = cos_theta / sin_theta;
		const double Delta = r2 - 2. * r + a2, Delta_1 = 1. / Delta;
		const double rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2), r2_a2_cos2_theta = r2 - a2 * cos2_theta;
		const double dydt4 = 2. * rho_4 * (Delta_1 * a2_r2 * r2_a2_cos2_theta * y[5] - 2. * a2 * r * sin_theta_cos_theta * y[6] * (1. - a * sin2_theta * y[7]) - Delta_1 * a * (2. * r4 + r2 * rho2 + a2 * r2_a2_cos2_theta) * sin2_theta * y[5] * y[7]);
		// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = dydt4 * y[4];
		// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = (-Delta * r2_a2_cos2_theta * gsl_pow_2(1. - a * sin2_theta * y[7]) - (r * (a2 * sin2_theta - r) + a2 * cos2_theta) * Delta_1 * rho4 * gsl_pow_2(y[5]) + 2. * a2 * sin_theta_cos_theta * rho4 * y[5] * y[6] + r * Delta * rho4 * gsl_pow_2(y[6]) + Delta * sin2_theta * r * rho4 * gsl_pow_2(y[7])) * rho_6 + dydt4 * y[5];
		// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = (2. * a2 * r * sin_theta_cos_theta - a2 * sin_theta_cos_theta * rho4 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho4 * y[5] * y[6] - 4. * a * r * sin_theta_cos_theta * a2_r2 * y[7] + sin_theta_cos_theta * (2. * a4 * r * sin4_theta + 4. * a2 * r * sin2_theta * rho2 + a2_r2 * rho4) * gsl_pow_2(y[7])) * rho_6 + dydt4 * y[6];
		// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = (-2. * a * r2_a2_cos2_theta * Delta_1 * y[5] + 4. * a * r * cot_theta * y[6] - 2. * Delta_1 * (r * rho4 - 2. * r2 * rho2 - r2_a2_cos2_theta * a2 * sin2_theta) * y[5] * y[7] - 2. * cot_theta * (rho4 + 2. * a2 * r * sin2_theta) * y[6] * y[7]) * rho_4 + dydt4 * y[7];
		return GSL_SUCCESS;
	}
	int KerrTLagrangianCircular(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = static_cast<Kerr *>(params);
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = 0.;	// dr/dt
		dydt[2] = 0.;	// d\theta/dt
		dydt[3] = y[7]; // d\phi/dt
		dydt[4] = 0.;
		dydt[5] = 0.;
		dydt[6] = 0.;
		dydt[7] = 0.;
		return GSL_SUCCESS;
	}
	int KerrTLagrangianHelical(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = static_cast<Kerr *>(params);
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = y[5]; // dr/dt
		dydt[2] = 0.;	// d\theta/dt = 0.
		dydt[3] = y[7]; // d\phi/dt
		const double r = y[1], r2 = gsl_pow_2(r), a = kerr->a_, a2 = kerr->a2_;
		const double sin2_theta = gsl_pow_2(sin(y[2]));
		const double Delta_1 = 1. / (r2 - 2. * r + a2);
		const double rho2 = r2 + a2 * gsl_pow_2(cos(y[2])), rho_2 = 1. / rho2;
		const double r_rho_2 = 2. * y[1] * rho_2;
		const double g00 = -(1. - r_rho_2), g11 = rho2 * Delta_1, g03 = -r_rho_2 * a * sin2_theta, g33 = (r2 + a2 - g03 * a) * sin2_theta;
		const double L_1 = y[4] / (g03 + g33 * y[7]), L_2 = gsl_pow_2(L_1);
		const double A = g33 * (g33 * L_2 + 1.), A_1 = 1. / A, B = 2. * g03 * (g33 * L_2 + 1.), C = g00 + g11 * gsl_pow_2(y[5]) + gsl_pow_2(g03) * L_2;
		const double dg00_dr = (2. - 4. * r2 * rho_2) * rho_2, dg11_dr = 2. * (r - g11 * (r - 1.)) * Delta_1, dg03_dr = -dg00_dr * a * sin2_theta, dg33_dr = 2. * r * sin2_theta - dg03_dr * a * sin2_theta;
		const double dA_dr = dg33_dr * (2. * g33 * L_2 + 1.), dB_dr = 2. * (dg03_dr * (g33 * L_2 + 1.) + g03 * dg33_dr * L_2), dC_dr = dg00_dr + dg11_dr * gsl_pow_2(y[5]) + 2. * g03 * dg03_dr * L_2;
		// d^2r/dt^2 = 0.
		dydt[5] = 0.;
		// d^2\theta/dt^2 = 0.
		dydt[6] = 0.;
		// d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
		dydt[7] = y[5] * A_1 * (-y[7] * dA_dr + 0.5 * (-dB_dr + (B * dB_dr - 2. * (A * dC_dr + C * dA_dr)) / (2. * A * y[7] + B)));
		// d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt = d((g03 + g33 * d\phi/dt) / L)/dr * dr/dt
		dydt[4] = (y[5] * (dg03_dr + dg33_dr * y[7]) + g33 * dydt[7]) * L_1;
		return GSL_SUCCESS;
	}
	int KerrTHamiltonianGeodesic(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = static_cast<Kerr *>(params);
		const double r = y[1], r2 = gsl_pow_2(y[1]), a = kerr->a_, a2 = kerr->a2_, a2_r2 = a2 + r2, pr2 = gsl_pow_2(y[5]), ptheta2 = gsl_pow_2(y[6]);
		const double E = 1. - y[4], E2 = gsl_pow_2(E), delta_E2 = (2. - y[4]) * y[4], L2 = gsl_pow_2(y[7]);
		const double sin_theta = abs(sin(y[2])), sin2_theta = gsl_pow_2(sin_theta), sin_2_theta = 1. / sin2_theta, sin_4_theta = gsl_pow_2(sin_2_theta), cos_theta = GSL_SIGN(y[2]) * cos(y[2]), cos2_theta = gsl_pow_2(cos_theta);
		const double Delta = a2_r2 - 2. * r, Delta_1 = 1. / Delta, Delta_2 = gsl_pow_2(Delta_1);
		const double rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho_4 = gsl_pow_2(rho_2);
		const double Q = ptheta2 + cos2_theta * (a2 * delta_E2 + L2 * sin_2_theta);
		const double R = -a2_r2 * r2 * delta_E2 + E2 * 2. * r * a2 + 2. * r * r2 - Delta * Q - (r2 - 2. * r) * L2 - 4. * r * a * E * y[7]; // R = gsl_pow_2(E * a2_r2 - a * y[7]) - Delta * (r2 + gsl_pow_2(y[7] - a * E) + Q);
		//[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
		dydt[0] = rho2 * Delta / (a2_r2 * (E * a2_r2 - a * y[7]) + a * Delta * (y[7] - a * E * sin2_theta));			// d\tau/dt
		dydt[1] = Delta * rho_2 * y[5] * dydt[0];																		// dr/dt
		dydt[2] = rho_2 * y[6] * dydt[0];																				// d\theta/dt
		dydt[3] = (2. * E * a * r - y[7] * (a2 - Delta * (1. + cos2_theta * sin_2_theta))) * Delta_1 * rho_2 * dydt[0]; // d\phi/dt
		dydt[4] = 0.;
		dydt[5] = ((a2 * cos2_theta + r * a2 * sin2_theta - r2) * pr2 + ((-2. * delta_E2 * r2 - a2 * delta_E2 + 3. * r - L2 - Q) * r + a2 * E2 + L2 - 2. * a * E * y[7] + Q) * rho2 * Delta_1 - ((r - 1.) * rho2 + Delta * r) * R * Delta_2) * rho_4 * dydt[0];
		dydt[6] = (-(Delta * pr2 - R * Delta_1) * a2 * rho_2 + a2 * delta_E2 + L2 * sin_2_theta + cos2_theta * sin_4_theta * L2) * sin_theta * cos_theta * rho_2 * dydt[0];
		dydt[7] = 0.;
		return GSL_SUCCESS;
	}
	int KerrTauLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = static_cast<Kerr *>(params);
		dydt[0] = y[4]; // dt/d\tau
		dydt[1] = y[5]; // dr/d\tau
		dydt[2] = y[6]; // d\theta/d\tau
		dydt[3] = y[7]; // d\phi/d\tau
		const double r = y[1], r2 = gsl_pow_2(r), a = kerr->a_, a2 = kerr->a2_, a4 = kerr->a4_, a2_r2 = a2 + r2;
		const double sin_theta = abs(sin(y[2])), sin2_theta = gsl_pow_2(sin_theta), sin4_theta = gsl_pow_4(sin_theta), cos_theta = GSL_SIGN(y[2]) * cos(y[2]), cos2_theta = gsl_pow_2(cos_theta), sin_theta_cos_theta = sin_theta * cos_theta, cot_theta = cos_theta / sin_theta;
		const double Delta = r2 - 2. * r + a2, Delta_1 = 1. / Delta;
		const double rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2), r2_a2_cos2_theta = r2 - a2 * cos2_theta;
		dydt[4] = -2. * rho_4 * (Delta_1 * a2_r2 * r2_a2_cos2_theta * y[5] * (y[4] - a * sin2_theta * y[7]) - 2. * a2 * r * sin_theta_cos_theta * y[6] * (y[4] - a * sin2_theta * y[7]) - 2. * Delta_1 * r2 * rho2 * a * sin2_theta * y[5] * y[7]);
		// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = (-Delta * r2_a2_cos2_theta * gsl_pow_2(y[4] - a * sin2_theta * y[7]) - (r * (a2 * sin2_theta - r) + a2 * cos2_theta) * Delta_1 * rho4 * gsl_pow_2(y[5]) + 2. * a2 * sin_theta_cos_theta * rho4 * y[5] * y[6] + r * Delta * rho4 * gsl_pow_2(y[6]) + Delta * sin2_theta * r * rho4 * gsl_pow_2(y[7])) * rho_6;
		// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = (2. * a2 * r * sin_theta_cos_theta * gsl_pow_2(y[4]) - a2 * sin_theta_cos_theta * rho4 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho4 * y[5] * y[6] - 4. * a * r * sin_theta_cos_theta * a2_r2 * y[4] * y[7] + sin_theta_cos_theta * (2. * a4 * r * sin4_theta + 4. * a2 * r * sin2_theta * rho2 + a2_r2 * rho4) * gsl_pow_2(y[7])) * rho_6;
		// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = (-2. * a * r2_a2_cos2_theta * Delta_1 * y[4] * y[5] + 4. * a * r * cot_theta * y[4] * y[6] - 2. * Delta_1 * (r * rho4 - 2. * r2 * rho2 - r2_a2_cos2_theta * a2 * sin2_theta) * y[5] * y[7] - 2. * cot_theta * (rho4 + 2. * a2 * r * sin2_theta) * y[6] * y[7]) * rho_4;
		return GSL_SUCCESS;
	}
	int KerrTauHamiltonianGeodesic(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = static_cast<Kerr *>(params);
		const double r = y[1], r2 = gsl_pow_2(y[1]), a = kerr->a_, a2 = kerr->a2_, a2_r2 = a2 + r2, pr2 = gsl_pow_2(y[5]), ptheta2 = gsl_pow_2(y[6]);
		const double E = 1. - y[4], E2 = gsl_pow_2(E), delta_E2 = (2. - y[4]) * y[4], L2 = gsl_pow_2(y[7]);
		const double sin_theta = abs(sin(y[2])), sin2_theta = gsl_pow_2(sin_theta), sin_2_theta = 1. / sin2_theta, sin_4_theta = gsl_pow_2(sin_2_theta), cos_theta = GSL_SIGN(y[2]) * cos(y[2]), cos2_theta = gsl_pow_2(cos_theta);
		const double Delta = a2_r2 - 2. * r, Delta_1 = 1. / Delta, Delta_2 = gsl_pow_2(Delta_1);
		const double rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho_4 = gsl_pow_2(rho_2);
		const double Q = ptheta2 + cos2_theta * (a2 * delta_E2 + L2 * sin_2_theta);
		const double R = -a2_r2 * r2 * delta_E2 + E2 * 2. * r * a2 + 2. * r * r2 - Delta * Q - (r2 - 2. * r) * L2 - 4. * r * a * E * y[7]; // R = gsl_pow_2(E * a2_r2 - a * y[7]) - Delta * (r2 + gsl_pow_2(y[7] - a * E) + Q);
		//[t,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
		dydt[0] = (a2_r2 * (E * a2_r2 - a * y[7]) * Delta_1 + a * (y[7] - a * E * sin2_theta)) * rho_2;		  // dt/d\tau
		dydt[1] = Delta * rho_2 * y[5];																		  // dr/d\tau
		dydt[2] = rho_2 * y[6];																				  // d\theta/d\tau
		dydt[3] = (2. * E * a * r - y[7] * (a2 - Delta * (1. + cos2_theta * sin_2_theta))) * Delta_1 * rho_2; // d\phi/d\tau
		dydt[4] = 0.;
		dydt[5] = ((a2 * cos2_theta + r * a2 * sin2_theta - r2) * pr2 + ((-2. * delta_E2 * r2 - a2 * delta_E2 + 3. * r - L2 - Q) * r + a2 * E2 + L2 - 2. * a * E * y[7] + Q) * rho2 * Delta_1 - ((r - 1.) * rho2 + Delta * r) * R * Delta_2) * rho_4;
		dydt[6] = (-(Delta * pr2 - R * Delta_1) * a2 * rho_2 + a2 * delta_E2 + L2 * sin_2_theta + cos2_theta * sin_4_theta * L2) * sin_theta * cos_theta * rho_2;
		dydt[7] = 0.;
		return GSL_SUCCESS;
	}
	int KerrTaubNutTLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		KerrTaubNUT *kerr_taub_nut = static_cast<KerrTaubNUT *>(params);
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = y[5]; // dr/dt
		dydt[2] = y[6]; // d\theta/dt
		dydt[3] = y[7]; // d\phi/dt
		const double r = y[1], r2 = gsl_pow_2(r), a = kerr_taub_nut->a_, a2 = kerr_taub_nut->a2_, l = kerr_taub_nut->l_, l2 = kerr_taub_nut->l2_;
		const double sin_theta = abs(sin(y[2])), sin_1_theta = 1. / sin_theta, sin2_theta = gsl_pow_2(sin_theta), cos_theta = GSL_SIGN(y[2]) * cos(y[2]), cos2_theta = gsl_pow_2(cos_theta);
		const double Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
		const double l_a_cos_theta = l + a * cos_theta, l_a_cos_theta2 = gsl_pow_2(l_a_cos_theta);
		const double rho2 = r2 + l_a_cos_theta2, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2);
		const double chi = a * sin2_theta - 2. * l * cos_theta, rho2_a_chi = r2 + l2 + a2;
		const double rho2_r_Delta = r * (r + 2. * l * l_a_cos_theta) - l_a_cos_theta2, rho2_a_cos_theta = (2. * r + l * l_a_cos_theta) * l_a_cos_theta - r2 * l;
		const double dydt4 = 2. * rho_4 * (Delta_1 * rho2_a_chi * rho2_r_Delta * y[5] * (1. - chi * y[7]) - sin_1_theta * chi * rho2_a_cos_theta * y[6] * (1. - chi * y[7]) - 2. * Delta_1 * rho2 * r * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * y[5] * y[7] - sin_1_theta * rho4 * l * (1. + cos2_theta) * y[6] * y[7]);
		// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = dydt4 * y[4];
		// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = -Delta * rho2_r_Delta * rho_6 * gsl_pow_2(1. - chi * y[7]) + (rho2_r_Delta - r * a2 * sin2_theta) * Delta_1 * rho_2 * gsl_pow_2(y[5]) + 2. * a * sin_theta * l_a_cos_theta * rho_2 * y[5] * y[6] + r * Delta * rho_2 * gsl_pow_2(y[6]) + Delta * r * sin2_theta * rho_2 * gsl_pow_2(y[7]) + dydt4 * y[5];
		// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = rho2_a_cos_theta * sin_theta * rho_6 * (a - 2. * rho2_a_chi * y[7]) - a * sin_theta * l_a_cos_theta * rho_2 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho_2 * y[5] * y[6] - (rho2 * (Delta * chi * l_a_cos_theta - gsl_pow_2(rho2_a_chi) * cos_theta) - 2. * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * rho2_a_chi * l_a_cos_theta) * sin_theta * rho_6 * gsl_pow_2(y[7]) + dydt4 * y[6];
		// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = -2. * a * rho2_r_Delta * Delta_1 * rho_4 * y[5] * (1. - chi * y[7]) + 2. * rho2_a_cos_theta * rho_4 * sin_1_theta * y[6] * (1. - chi * y[7]) - 2. * (1. - a2 * sin2_theta * Delta_1) * r * rho_2 * y[5] * y[7] - 2. * cos_theta * sin_1_theta * y[6] * y[7] + dydt4 * y[7];
		return GSL_SUCCESS;
	}
	int KerrTaubNutTLagrangianHelical(double t, const double y[], double dydt[], void *params) {
		KerrTaubNUT *kerr_taub_nut = static_cast<KerrTaubNUT *>(params);
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = y[5]; // dr/dt
		dydt[2] = 0.;	// d\theta/dt = 0.
		dydt[3] = y[7]; // d\phi/dt
		const double r = y[1], r2 = gsl_pow_2(r), a = kerr_taub_nut->a_, a2 = kerr_taub_nut->a2_, l = kerr_taub_nut->l_, l2 = kerr_taub_nut->l2_;
		const double sin2_theta = gsl_pow_2(sin(y[2])), cos_theta = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
		const double l_a_cos_theta2 = gsl_pow_2(l + a * cos_theta);
		const double rho2 = r2 + l_a_cos_theta2, rho_2 = 1. / rho2;
		const double chi = a * sin2_theta - 2. * l * cos_theta, rho2_a_chi = r2 + l2 + a2;
		const double g00 = -rho_2 * (Delta - a2 * sin2_theta), g11 = rho2 * Delta_1, g03 = -g00 * chi - a * sin2_theta, g33 = -g03 * chi + rho2_a_chi * sin2_theta;
		const double L_1 = y[4] / (g03 + g33 * y[7]), L_2 = gsl_pow_2(L_1);
		const double A = g33 * (g33 * L_2 + 1.), A_1 = 1. / A, B = 2. * g03 * (g33 * L_2 + 1.), C = g00 + g11 * gsl_pow_2(y[5]) + gsl_pow_2(g03) * L_2;
		const double dg00_dr = -2. * rho_2 * ((r - 1.) + g00 * r), dg11_dr = 2. * Delta_1 * (r - g11 * (r - 1.)), dg03_dr = -dg00_dr * chi, dg33_dr = 2. * r * sin2_theta - dg03_dr * chi;
		const double dA_dr = dg33_dr * (2. * g33 * L_2 + 1.), dB_dr = 2. * (dg03_dr * (g33 * L_2 + 1.) + g03 * dg33_dr * L_2), dC_dr = dg00_dr + dg11_dr * gsl_pow_2(y[5]) + 2. * g03 * dg03_dr * L_2;
		// d^2r/dt^2 = 0.
		dydt[5] = 0.;
		// d^2\theta/dt^2 = 0.
		dydt[6] = 0.;
		// d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
		dydt[7] = y[5] * A_1 * (-y[7] * dA_dr + 0.5 * (-dB_dr + (B * dB_dr - 2. * (A * dC_dr + C * dA_dr)) / (2. * A * y[7] + B)));
		// d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt = d((g03 + g33 * d\phi/dt) / L)/dr * dr/dt
		dydt[4] = (y[5] * (dg03_dr + dg33_dr * y[7]) + g33 * dydt[7]) * L_1;
		return GSL_SUCCESS;
	}
	int KerrTaubNutTHamiltonianGeodesic(double t, const double y[], double dydt[], void *params) { // TODO:
		return GSL_FAILURE;
	}
	int KerrTaubNutTauLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		KerrTaubNUT *kerr_taub_nut = static_cast<KerrTaubNUT *>(params);
		dydt[0] = y[4]; // dt/d\tau
		dydt[1] = y[5]; // dr/d\tau
		dydt[2] = y[6]; // d\theta/d\tau
		dydt[3] = y[7]; // d\phi/d\tau
		const double r = y[1], r2 = gsl_pow_2(r), a = kerr_taub_nut->a_, a2 = kerr_taub_nut->a2_, l = kerr_taub_nut->l_, l2 = kerr_taub_nut->l2_;
		const double sin_theta = abs(sin(y[2])), sin_1_theta = 1. / sin_theta, sin2_theta = gsl_pow_2(sin_theta), cos_theta = GSL_SIGN(y[2]) * cos(y[2]), cos2_theta = gsl_pow_2(cos_theta);
		const double Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
		const double l_a_cos_theta = l + a * cos_theta, l_a_cos_theta2 = gsl_pow_2(l_a_cos_theta);
		const double rho2 = r2 + l_a_cos_theta2, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2);
		const double chi = a * sin2_theta - 2. * l * cos_theta, rho2_a_chi = r2 + l2 + a2;
		const double rho2_r_Delta = r * (r + 2. * l * l_a_cos_theta) - l_a_cos_theta2, rho2_a_cos_theta = (2. * r + l * l_a_cos_theta) * l_a_cos_theta - r2 * l;
		// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = -2. * rho_4 * (Delta_1 * rho2_a_chi * rho2_r_Delta * y[5] * (y[4] - chi * y[7]) - sin_1_theta * chi * rho2_a_cos_theta * y[6] * (y[4] - chi * y[7]) - 2. * Delta_1 * rho2 * r * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * y[5] * y[7] - sin_1_theta * rho4 * l * (1. + cos2_theta) * y[6] * y[7]);
		// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = -Delta * rho2_r_Delta * rho_6 * gsl_pow_2(y[4] - chi * y[7]) + (rho2_r_Delta - r * a2 * sin2_theta) * Delta_1 * rho_2 * gsl_pow_2(y[5]) + 2. * a * sin_theta * l_a_cos_theta * rho_2 * y[5] * y[6] + r * Delta * rho_2 * gsl_pow_2(y[6]) + Delta * r * sin2_theta * rho_2 * gsl_pow_2(y[7]);
		// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = rho2_a_cos_theta * sin_theta * rho_6 * y[4] * (a * y[4] - 2. * rho2_a_chi * y[7]) - a * sin_theta * l_a_cos_theta * rho_2 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho_2 * y[5] * y[6] - (rho2 * (Delta * chi * l_a_cos_theta - gsl_pow_2(rho2_a_chi) * cos_theta) - 2. * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * rho2_a_chi * l_a_cos_theta) * sin_theta * rho_6 * gsl_pow_2(y[7]);
		// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = -2. * a * rho2_r_Delta * Delta_1 * rho_4 * y[5] * (y[4] - chi * y[7]) + 2. * rho2_a_cos_theta * rho_4 * sin_1_theta * y[6] * (y[4] - chi * y[7]) - 2. * (1. - a2 * sin2_theta * Delta_1) * r * rho_2 * y[5] * y[7] - 2. * cos_theta * sin_1_theta * y[6] * y[7];
		return GSL_SUCCESS;
	}
	int Jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
		return GSL_SUCCESS;
	}
} // namespace SBody
