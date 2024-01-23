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
#include <gsl/gsl_vector.h>

#include "IO.h"
#include "Utility.h"

namespace SBody {
	int Metric::LocalInertialFrame(const double position[], gsl_matrix *coordinate, const double timelike[]) {
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
				if (timelike != nullptr) {
					gsl_vector_set(coordinate_row, 1, timelike[1]);
					gsl_vector_set(coordinate_row, 2, timelike[2]);
					gsl_vector_set(coordinate_row, 3, timelike[3]);
				}
			} else {
				gsl_matrix_memcpy(product_LU, product);
				gsl_linalg_LU_decomp(product_LU, permutation, &signum);
				gsl_linalg_LU_svx(product_LU, permutation, coordinate_row);
			}
			gsl_vector_scale(coordinate_row, 1. / sqrt(abs(DotProduct(position, coordinate_row->data, coordinate_row->data, 4))));
			gsl_blas_dsymv(CblasUpper, 1., metric, coordinate_row, 0., product_row);
		}
		return isnan(gsl_matrix_get(coordinate, 3, 0)) ? 1 : 0;
	}
	double Metric::Redshift(const double y[], const double photon[], time_system time) {
		if (time == T) {
			const double u[4] = {1., y[5], y[6], y[7]}, v[4] = {1., photon[5], photon[6], photon[7]};
			return -DotProduct(y, u, v, 4) / (y[4] * photon[4]);
		}
		return -DotProduct(y, y + 4, photon + 4, 4);
	}

	Newton::Newton(int PN) : PN_(PN){};
	std::string Newton::Name() {
		return "Newton";
	}
	int Newton::MetricTensor(const double position[], gsl_matrix *metric) { // FIXME: Check
		gsl_matrix_set_zero(metric);
		gsl_matrix_set(metric, 0, 0, -1.);
		gsl_matrix_set(metric, 1, 1, 1.);
		gsl_matrix_set(metric, 2, 2, gsl_pow_2(position[1]));
		gsl_matrix_set(metric, 3, 3, gsl_pow_2(position[1] * sin(position[2])));
		return 0;
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
	int Newton::BaseToHamiltonian(double y[]) {
		return GSL_FAILURE;
	}
	int Newton::HamiltonianToBase(double y[]) {
		return GSL_FAILURE;
	}
	double Newton::Energy(const double y[], coordinate_system coordinate) {
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
	double Newton::AngularMomentum(const double y[], coordinate_system coordinate) {
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
	double Newton::CarterConstant(const double y[], const double mu2, coordinate_system coordinate) { // FIXME: not verified!
		return gsl_pow_2(AngularMomentum(y, coordinate));
	}
	double Newton::Redshift(const double y[], const double photon[], time_system time) {
		const double delta_epsilon = 1. - 2. / y[1];
		return (1. - DotProduct(y, y + 4, photon + 4, 3) / sqrt(delta_epsilon)) / sqrt(delta_epsilon - gsl_pow_2(y[5]) - gsl_pow_2(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(sin(y[2]) * y[7])));
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
	std::unique_ptr<Integrator> Newton::GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion) {
		if (time != T || coordinate != LAGRANGIAN || motion != GEODESIC) {
			PrintlnError("Newton::GetIntegrator() {}, {}, {} invaild", time, coordinate, motion);
			SBodyExit(GSL_EINVAL);
		}
		return std::make_unique<Integrator>(
			[](double t, const double y[], double dydt[], void *params) -> int {
				int PN = *static_cast<int *>(params);
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
			},
			[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; },
			const_cast<int *>(&this->PN_));
	}
	PN1::PN1(double fSP) : Newton(1), PN1_(fSP){};
	std::unique_ptr<Integrator> PN1::GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion) {
		if (time != T || coordinate != LAGRANGIAN || motion != GEODESIC) {
			PrintlnError("PN1::GetIntegrator() {}, {}, {} invaild", time, coordinate, motion);
			SBodyExit(GSL_EINVAL);
		}
		return std::make_unique<Integrator>(
			[](double t, const double y[], double dydt[], void *params) -> int {
				const double PN1 = *static_cast<const double *>(params);
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
			},
			[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; },
			const_cast<double *>(&this->PN1_));
	}

	Schwarzschild::Schwarzschild() {}
	std::string Schwarzschild::Name() {
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
	int Schwarzschild::BaseToHamiltonian(double y[]) {
		const double dtdtau = 1. / y[4], r2 = gsl_pow_2(y[1]);
		y[4] = (y[4] - 1. + 2. / y[1]) * dtdtau; // 1 + p_t
		y[5] *= y[1] / (y[1] - 2.) * dtdtau;
		y[6] *= r2 * dtdtau;
		y[7] *= r2 * gsl_pow_2(sin(y[2])) * dtdtau;
		return 0;
	}
	int Schwarzschild::HamiltonianToBase(double y[]) {
		const double r_1 = 1. / y[1], g00 = 1. - 2. * r_1;
		y[4] = -g00 / (y[4] - 1.);
		y[5] *= g00 * y[4];
		y[6] *= gsl_pow_2(r_1) * y[4];
		y[7] *= gsl_pow_2(r_1 / sin(y[2])) * y[4];
		return 0;
	}
	double Schwarzschild::Energy(const double y[], coordinate_system coordinate) {
		if (coordinate == BASE)
			return (y[1] - 2.) / (y[1] * y[4]);
		else if (coordinate == LAGRANGIAN)
			return 1.; // TODO:
		else
			return 1. - y[4];
	}
	double Schwarzschild::AngularMomentum(const double y[], coordinate_system coordinate) {
		if (coordinate == BASE)
			return gsl_pow_2(y[1]) * y[7] / y[4];
		else if (coordinate == LAGRANGIAN)
			return 1.; // TODO:
		else
			return y[7];
	}
	double Schwarzschild::CarterConstant(const double y[], const double mu2, coordinate_system coordinate) {
		if (coordinate == BASE)
			return gsl_pow_4(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(y[7] * cos(y[2]) * sin(y[2]))) / gsl_pow_2(y[4]);
		else if (coordinate == LAGRANGIAN)
			return 1.; // TODO:
		else
			return gsl_pow_2(y[6]) + gsl_pow_2(y[7] / tan(y[2]));
	}
	int Schwarzschild::NormalizeTimelikeGeodesic(double y[]) {
		const double g00 = 1. - 2. / y[1];
		if (g00 <= 0)
			return 1;
		y[4] = sqrt(g00 - (gsl_pow_2(y[5]) / g00 + gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
		return std::isnan(y[4]);
	}
	int Schwarzschild::NormalizeNullGeodesic(double y[], double frequency) {
		const double g00 = 1. - 2. / y[1];
		if (g00 <= 0)
			return 1;
		const double coefficient = GSL_SIGN(frequency) * g00 / sqrt(gsl_pow_2(y[5]) + g00 * (gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
		y[4] = frequency;
		y[5] *= coefficient;
		y[6] *= coefficient;
		y[7] *= coefficient;
		return 0;
	}
	std::unique_ptr<Integrator> Schwarzschild::GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion) {
		if (time == T) {
			if (coordinate == LAGRANGIAN) {
				if (motion == GEODESIC)
					return std::make_unique<Integrator>(
						[](double t, const double y[], double dydt[], void *params) -> int {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = y[6]; // d\theta/dt
							dydt[3] = y[7]; // d\phi/dt
							const double r = y[1], sint = abs(sin(y[2])), cost = GSL_SIGN(y[2]) * cos(y[2]);
							const double r2m = r - 2., r3m = r - 3.;
							const double r2mr_1 = 1. / (r2m * r);
							// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
							dydt[4] = 2. * y[5] * r2mr_1 * y[4];
							// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[5] = -r2m / gsl_pow_3(r) + 3. * r2mr_1 * gsl_pow_2(y[5]) + r2m * (gsl_pow_2(y[6]) + gsl_pow_2(sint * y[7]));
							// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[6] = -2. * r3m * r2mr_1 * y[5] * y[6] + sint * cost * gsl_pow_2(y[7]);
							// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[7] = -2. * (r3m * r2mr_1 * y[5] + cost / sint * y[6]) * y[7];
							return GSL_SUCCESS;
						},
						[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; });
				if (motion == CIRCULAR)
					return std::make_unique<Integrator>(
						[](double t, const double y[], double dydt[], void *params) -> int {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = y[6]; // d\theta/dt
							dydt[3] = y[7]; // d\phi/dt
							dydt[4] = 0.;
							dydt[5] = 0.;
							dydt[6] = 0.;
							dydt[7] = 0.;
							return GSL_SUCCESS;
						},
						[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; });
				if (motion == RIAF)
					return std::make_unique<Integrator>(
						[](double t, const double y[], double dydt[], void *params) -> int { return GSL_FAILURE; },
						[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; });
				if (motion == HELICAL)
					return std::make_unique<Integrator>(
						[](double t, const double y[], double dydt[], void *params) -> int {
							const double g00 = 1. - 2. / y[1], sin_theta = abs(sin(y[2]));
							if (g00 <= 0)
								return GSL_FAILURE;
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = 0.;	// d\theta/dt = 0.
							dydt[3] = y[7]; // d\phi/dt
							// d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt
							dydt[4] = ((1. + gsl_pow_2(y[5] / g00)) / gsl_pow_2(y[1]) + y[1] * gsl_pow_2(sin_theta * y[7])) * y[5] / sqrt(g00 - (gsl_pow_2(y[5]) / g00 + gsl_pow_2(y[1] * sin_theta * y[7])));
							// d^2r/dt^2 = 0.
							dydt[5] = 0.;
							// d^2\theta/dt^2 = 0.
							dydt[6] = 0.;
							// d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
							dydt[7] = -2. * y[7] / y[1] * y[5];
							return GSL_SUCCESS;
						},
						[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; });
			}
			if (coordinate == HAMILTONIAN && motion == GEODESIC)
				return std::make_unique<Integrator>(
					[](double t, const double y[], double dydt[], void *params) -> int {
						const double r_1 = 1. / y[1], r_2 = gsl_pow_2(r_1), g00 = 1. - 2. * r_1, E = 1. - y[4], L2 = gsl_pow_2(y[7]);
						const double sint_1 = 1. / abs(sin(y[2])), sint_2 = gsl_pow_2(sint_1);
						//[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
						dydt[0] = g00 / E;						 // d\tau/dt
						dydt[1] = g00 * y[5] * dydt[0];			 // dr/dt
						dydt[2] = y[6] * r_2 * dydt[0];			 // d\theta/dt
						dydt[3] = y[7] * sint_2 * r_2 * dydt[0]; // d\phi/dt
						dydt[4] = 0.;
						dydt[5] = (-(gsl_pow_2(y[5]) + gsl_pow_2(E) / gsl_pow_2(g00)) + (gsl_pow_2(y[6]) + L2 * sint_2) * r_1) * r_2 * dydt[0];
						dydt[6] = sint_2 * L2 * GSL_SIGN(y[2]) * cos(y[2]) * sint_1 * r_2 * dydt[0];
						dydt[7] = 0.;
						return GSL_SUCCESS;
					},
					[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; });
		}
		if (time == TAU && coordinate == LAGRANGIAN && motion == GEODESIC)
			return std::make_unique<Integrator>(
				[](double t, const double y[], double dydt[], void *params) -> int {
					dydt[0] = y[4]; // dt/d\tau
					dydt[1] = y[5]; // dr/d\tau
					dydt[2] = y[6]; // d\theta/d\tau
					dydt[3] = y[7]; // d\phi/d\tau
					const double r = y[1], sint = abs(sin(y[2])), cost = GSL_SIGN(y[2]) * cos(y[2]);
					const double r2m = r - 2., r_1 = 1. / r;
					const double r2mr_1 = 1. / (r2m * r);
					// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
					dydt[4] = -2. * y[5] * r2mr_1 * y[4];
					// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
					dydt[5] = -r2m * gsl_pow_3(r_1) * gsl_pow_2(y[4]) + r2mr_1 * gsl_pow_2(y[5]) + r2m * (gsl_pow_2(y[6]) + gsl_pow_2(sint * y[7]));
					// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
					dydt[6] = -2. * r_1 * y[5] * y[6] + sint * cost * gsl_pow_2(y[7]);
					// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
					dydt[7] = -2. * (r_1 * y[5] + cost / sint * y[6]) * y[7];
					return GSL_SUCCESS;
				},
				[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; });
		return std::make_unique<Integrator>(
			[](double t, const double y[], double dydt[], void *params) -> int { return GSL_FAILURE; },
			[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; });
	}

	Kerr::Kerr(double spin) : a_(spin), a2_(a_ * a_), a4_(a2_ * a2_) {}
	std::string Kerr::Name() {
		return "Kerr";
	}
	int Kerr::MetricTensor(const double position[], gsl_matrix *metric) {
		const double r2 = gsl_pow_2(position[1]), Delta = r2 - 2. * position[1] + a2_, rho2 = gsl_pow_2(position[1]) + a2_ * gsl_pow_2(cos(position[2])), mr_rho2 = 2. * position[1] / rho2, sin2_theta = gsl_pow_2(sin(position[2]));
		gsl_matrix_set_zero(metric);
		gsl_matrix_set(metric, 0, 0, -(1. - mr_rho2));
		gsl_matrix_set(metric, 0, 3, -mr_rho2 * a_ * sin2_theta);
		gsl_matrix_set(metric, 1, 1, rho2 / Delta);
		gsl_matrix_set(metric, 2, 2, rho2);
		gsl_matrix_set(metric, 3, 0, -mr_rho2 * a_ * sin2_theta);
		gsl_matrix_set(metric, 3, 3, (r2 + a2_ * (1. + mr_rho2 * sin2_theta)) * sin2_theta);
		return rho2 == 0. || Delta == 0. ? 1 : 0;
	}
	double Kerr::DotProduct(const double position[], const double x[], const double y[], const size_t dimension) {
		const double r = position[1], r2 = gsl_pow_2(r);
		const double sint = sin(position[2]), sint2 = gsl_pow_2(sint), cost = cos(position[2]), cost2 = gsl_pow_2(cost);
		const double Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * cost2;
		const double mr_rho2 = 2. * r / rho2;
		if (dimension == 3)
			return (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + ((r2 + a2_) * sint2 + mr_rho2 * a2_ * gsl_pow_2(sint2)) * x[3] * y[3];
		return (mr_rho2 - 1.) * x[0] * y[0] - mr_rho2 * a_ * sint2 * (x[0] * y[3] + x[3] * y[0]) + (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + ((r2 + a2_) * sint2 + mr_rho2 * a2_ * gsl_pow_2(sint2)) * x[3] * y[3];
	}
	double Kerr::DistanceSquare(const double x[], const double y[], const size_t dimension) {
		const double r = x[1], r2 = gsl_pow_2(r), d0 = x[0] - y[0], d3 = PhiDifference(x[3] - y[3]);
		const double sint = sin(x[2]), sint2 = gsl_pow_2(sint), cost = cos(x[2]), cost2 = gsl_pow_2(cost);
		const double Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * cost2;
		const double mr_rho2 = 2. * r / rho2;
		if (dimension == 3)
			return (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + ((r2 + a2_) * sint2 + 2. * r * a2_ * gsl_pow_2(sint2) / rho2) * gsl_pow_2(d3);
		return (mr_rho2 - 1.) * gsl_pow_2(d0) - 2. * mr_rho2 * a_ * sint2 * d0 * d3 + (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + ((r2 + a2_) * sint2 + mr_rho2 * a2_ * gsl_pow_2(sint2)) * gsl_pow_2(d3);
	}
	int Kerr::BaseToHamiltonian(double y[]) {
		const double dtdtau = 1. / y[4], r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2]));
		const double Delta = r2 - 2. * y[1] + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double mr_rho2 = 2. * y[1] / rho2;
		y[4] = (y[4] - 1. + mr_rho2 * (1. - a_ * sint2 * y[7])) * dtdtau; // 1 + p_t
		y[5] *= rho2 / Delta * dtdtau;
		y[6] *= rho2 * dtdtau;
		y[7] = (-mr_rho2 * a_ + (a2_ + r2 + mr_rho2 * a2_ * sint2) * y[7]) * sint2 * dtdtau;
		return 0;
	}
	int Kerr::HamiltonianToBase(double y[]) {
		const double pt = y[4] - 1., r2 = gsl_pow_2(y[1]);
		const double Delta = r2 - 2. * y[1] + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double rho_2 = 1. / rho2, mr_rho2 = 2. * y[1] * rho_2;
		y[4] = -Delta / ((Delta + mr_rho2 * (a2_ + r2)) * pt + mr_rho2 * a_ * y[7]);
		y[5] *= Delta * y[4] * rho_2;
		y[6] *= y[4] * rho_2;
		y[7] = (-mr_rho2 * a_ * pt + (1. - mr_rho2) / gsl_pow_2(sin(y[2])) * y[7]) / Delta * y[4];
		return 0;
	}
	double Kerr::Energy(const double y[], coordinate_system coordinate) {
		if (coordinate == BASE)
			return (1. - 2. * y[1] / (gsl_pow_2(y[1]) + a2_ * gsl_pow_2(cos(y[2]))) * (1. - a_ * gsl_pow_2(sin(y[2])) * y[7])) / y[4];
		else if (coordinate == LAGRANGIAN)
			return 1.; // TODO:
		else
			return 1. - y[4];
	}
	double Kerr::AngularMomentum(const double y[], coordinate_system coordinate) {
		const double r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2]));
		const double mr_rho2 = 2. * y[1] / (r2 + a2_ * gsl_pow_2(cos(y[2])));
		if (coordinate == BASE)
			return (-mr_rho2 * a_ + (a2_ + r2 + mr_rho2 * a2_ * sint2) * y[7]) * sint2 / y[4];
		else if (coordinate == LAGRANGIAN)
			return 1.; // TODO:
		else
			return y[7];
	}
	double Kerr::CarterConstant(const double y[], const double mu2, coordinate_system coordinate) {
		const double r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2])), cost2 = gsl_pow_2(cos(y[2]));
		const double rho2 = r2 + a2_ * cost2;
		const double mr_rho2 = 2. * y[1] / rho2;
		if (coordinate == BASE)
			return mu2 * cost2 * a2_ + (gsl_pow_2(rho2 * y[6]) + cost2 * (gsl_pow_2(-mr_rho2 * a_ + (a2_ + r2 + mr_rho2 * a2_ * sint2) * y[7]) * sint2 - a2_ * gsl_pow_2(mr_rho2 * (1. - a_ * sint2 * y[7]) - 1.))) / gsl_pow_2(y[4]);
		else if (coordinate == LAGRANGIAN)
			return 1.; // TODO:
		else
			return gsl_pow_2(y[6]) + gsl_pow_2(cos(y[2])) * (a2_ * (mu2 - gsl_pow_2(1. - y[4])) + gsl_pow_2(y[7] / sin(y[2])));
	}
	int Kerr::NormalizeTimelikeGeodesic(double y[]) {
		const double r = y[1], r2 = gsl_pow_2(r), a2r2 = a2_ + r2;
		const double sint2 = gsl_pow_2(sin(y[2])), sint4 = gsl_pow_2(sint2);
		const double Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double mr_rho2 = 2. * r / rho2;
		y[7] += 2. * a_ * r / (gsl_pow_2(a2r2) - a2_ * Delta * sint2);
		y[4] = sqrt(1. - mr_rho2 + 2. * mr_rho2 * a_ * sint2 * y[7] - (rho2 / (r2 - 2. * r + a2_) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + (a2r2 * sint2 + mr_rho2 * a2_ * sint4) * gsl_pow_2(y[7])));
		return std::isnan(y[4]);
	}
	int Kerr::NormalizeNullGeodesic(double y[], double frequency) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sint2 = gsl_pow_2(sin(y[2])), sint4 = gsl_pow_2(sint2);
		const double rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double mr_rho2 = 2. * r / rho2;
		const double effa = rho2 / (r2 - 2. * r + a2_) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2_ + r2) * sint2 + mr_rho2 * a2_ * sint4) * gsl_pow_2(y[7]);
		const double effb = -2. * mr_rho2 * a_ * sint2 * y[7];
		const double effc = mr_rho2 - 1.;
		const double eff = GSL_SIGN(frequency) * 0.5 * (-effb + sqrt(gsl_pow_2(effb) - 4. * effa * effc)) / effa;
		y[4] = frequency;
		y[5] *= eff;
		y[6] *= eff;
		y[7] *= eff;
		return 0;
	}
	std::unique_ptr<Integrator> Kerr::GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion) {
		if (motion != GEODESIC)
			PrintlnError("Motion not supported!");
		if (time == T) {
			if (coordinate == LAGRANGIAN)
				return std::make_unique<Integrator>(KerrTLagrangianGeodesic, Jacobian, this);
			else if (coordinate == HAMILTONIAN)
				return std::make_unique<Integrator>(KerrTHamiltonianGeodesic, Jacobian, this);
		} else if (time == TAU) {
			if (coordinate == LAGRANGIAN)
				return std::make_unique<Integrator>(KerrTauLagrangianGeodesic, Jacobian, this);
			else if (coordinate == HAMILTONIAN)
				return std::make_unique<Integrator>(KerrTauHamiltonianGeodesic, Jacobian, this);
		}
		PrintlnError("Time system not supported!");
		return std::make_unique<Integrator>(
			[](double t, const double y[], double dydt[], void *params) -> int { return GSL_FAILURE; },
			[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; });
	}

	KerrTaubNUT::KerrTaubNUT(double spin, double NUT) : Kerr(spin), l_(NUT), l2_(l_ * l_), l4_(l2_ * l2_) {}
	std::string KerrTaubNUT::Name() {
		return "Kerr-Taub-NUT";
	}
	int KerrTaubNUT::MetricTensor(const double position[], gsl_matrix *metric) {
		return 1;
	}
	double KerrTaubNUT::DotProduct(const double position[], const double x[], const double y[], const size_t dimension) {
		const double r = position[1], sint = sin(position[2]), cost = cos(position[2]);
		const double r2 = gsl_pow_2(r), sint2 = gsl_pow_2(sint);
		const double Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + gsl_pow_2(l_ + a_ * cost), chi = a_ * sint2 - 2. * l_ * cost;
		const double rho_2 = 1. / rho2;
		if (dimension == 3)
			return (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + (gsl_pow_2(rho2 + a_ * chi) * sint2 - chi * chi * Delta) * rho_2 * x[3] * y[3];
		return (a2_ * sint2 - Delta) * rho_2 * x[0] * y[0] - 2. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) * rho_2 * (x[0] * y[3] + x[3] * y[0]) + (rho2 / Delta) * x[1] * y[1] + rho2 * x[2] * y[2] + (gsl_pow_2(rho2 + a_ * chi) * sint2 - chi * chi * Delta) * rho_2 * x[3] * y[3];
	}
	double KerrTaubNUT::DistanceSquare(const double x[], const double y[], const size_t dimension) {
		const double r = x[1], sint = sin(x[2]), cost = cos(x[2]), d0 = x[0] - y[0], d3 = PhiDifference(x[3] - y[3]);
		const double r2 = gsl_pow_2(r), sint2 = gsl_pow_2(sint);
		const double Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + gsl_pow_2(l_ + a_ * cost), chi = a_ * sint2 - 2. * l_ * cost;
		const double rho_2 = 1. / rho2;
		if (dimension == 3)
			return (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + (gsl_pow_2(rho2 + a_ * chi) * sint2 - gsl_pow_2(chi) * Delta) * rho_2 * gsl_pow_2(d3);
		return (a2_ * sint2 - Delta) * rho_2 * gsl_pow_2(d0) - 4. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) * rho_2 * d0 * d3 + (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + (gsl_pow_2(rho2 + a_ * chi) * sint2 - gsl_pow_2(chi) * Delta) * rho_2 * gsl_pow_2(d3);
	}
	int KerrTaubNUT::BaseToHamiltonian(double y[]) {
		const double dtdtau = 1. / y[4], r = y[1], r2 = gsl_pow_2(r);
		const double sint2 = gsl_pow_2(sin(y[2])), cost = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		const double rho2 = r2 + gsl_pow_2(l_ + a_ * cost), rho_2 = 1. / rho2;
		y[4] = (y[4] * rho2 - Delta + a2_ * sint2 - 2. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) * y[7]) * rho_2 * dtdtau; // 1 + p_t
		y[5] *= rho2 / Delta * dtdtau;
		y[6] *= rho2 * dtdtau;
		y[7] = (-2. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) + (gsl_pow_2(r2 + l2_ + a2_) * sint2 - gsl_pow_2(a_ * sint2 - 2. * l_ * cost) * Delta) * y[7]) * rho_2 * dtdtau;
		return 0;
	}
	int KerrTaubNUT::HamiltonianToBase(double y[]) {
		const double pt = y[4] - 1., r = y[1], r2 = gsl_pow_2(r);
		const double sint2 = gsl_pow_2(sin(y[2])), cost = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		const double rho2 = r2 + gsl_pow_2(l_ + a_ * cost), rho_2 = 1. / rho2;
		y[4] = Delta * rho2 * sint2 / ((Delta * gsl_pow_2(a_ * sint2 - 2. * l_ * cost) - gsl_pow_2(r2 + l2_ + a2_) * sint2) * pt - 2. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) * y[7]);
		y[5] *= Delta * rho_2 * y[4];
		y[6] *= rho_2 * y[4];
		y[7] = (-2. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) * pt + (Delta - a2_ * sint2) * y[7]) / (Delta * rho2 * sint2) * y[4];
		return 0;
	}
	double KerrTaubNUT::Energy(const double y[], coordinate_system coordinate) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sint2 = gsl_pow_2(sin(y[2])), cost = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		if (coordinate == BASE)
			return (Delta - a2_ * sint2 + 2. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) * y[7]) / ((r2 + gsl_pow_2(l_ + a_ * cost)) * y[4]);
		else if (coordinate == LAGRANGIAN)
			return 1.; // TODO:
		else
			return 1. - y[4];
	}
	double KerrTaubNUT::AngularMomentum(const double y[], coordinate_system coordinate) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sint2 = gsl_pow_2(sin(y[2])), cost = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		if (coordinate == BASE)
			return (-2. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) + (gsl_pow_2(r2 + l2_ + a2_) * sint2 - gsl_pow_2(a_ * sint2 - 2. * l_ * cost) * Delta) * y[7]) / ((r2 + gsl_pow_2(l_ + a_ * cost)) * y[4]);
		else if (coordinate == LAGRANGIAN)
			return 1.; // TODO:
		else
			return y[7];
	}
	double KerrTaubNUT::CarterConstant(const double y[], const double mu2, coordinate_system coordinate) { // TODO:
		const double r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2])), cost2 = gsl_pow_2(cos(y[2]));
		const double rho2 = r2 + a2_ * cost2;
		const double mr_rho2 = 2. * y[1] / rho2;
		if (coordinate == BASE)
			return mu2 * cost2 * a2_ + (gsl_pow_2(rho2 * y[6]) + cost2 * (gsl_pow_2(-mr_rho2 * a_ + (a2_ + r2 + mr_rho2 * a2_ * sint2) * y[7]) * sint2 - a2_ * gsl_pow_2(mr_rho2 * (1. - a_ * sint2 * y[7]) - 1.))) / gsl_pow_2(y[4]);
		else if (coordinate == LAGRANGIAN)
			return 1.;
		else
			return gsl_pow_2(y[6]) + gsl_pow_2(cos(y[2])) * (a2_ * (mu2 - gsl_pow_2(1. - y[4])) + gsl_pow_2(y[7] / sin(y[2])));
	}
	int KerrTaubNUT::NormalizeTimelikeGeodesic(double y[]) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sint2 = gsl_pow_2(sin(y[2])), cost = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		const double rho2 = r2 + gsl_pow_2(l_ + a_ * cost);
		y[7] += 2. * a_ * r / (gsl_pow_2(a2_ + r2) - a2_ * Delta * sint2);
		y[4] = sqrt(((Delta - a2_ * sint2) + 4. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) * y[7] - (gsl_pow_2(r2 + l2_ + a2_) * sint2 - gsl_pow_2(a_ * sint2 - 2. * l_ * cost) * Delta) * gsl_pow_2(y[7])) / rho2 - rho2 * (gsl_pow_2(y[5]) / Delta + gsl_pow_2(y[6])));
		return std::isnan(y[4]);
	}
	int KerrTaubNUT::NormalizeNullGeodesic(double y[], double frequency) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sint2 = gsl_pow_2(sin(y[2])), cost = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		const double lacost = l_ + a_ * cost;
		const double rho2 = r2 + gsl_pow_2(lacost), rho_2 = 1. / rho2;
		const double chi = a_ * sint2 - 2. * l_ * cost, rho2achi = r2 + l2_ + a2_;
		const double effa = rho2 / Delta * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + rho_2 * (gsl_pow_2(rho2achi) * sint2 - gsl_pow_2(chi) * Delta) * gsl_pow_2(y[7]);
		const double effb = -4. * rho_2 * ((r + l2_) * chi + l_ * cost * rho2achi) * y[7];
		const double effc = -rho_2 * (Delta - a2_ * sint2);
		const double eff = GSL_SIGN(frequency) * 0.5 * (-effb + sqrt(gsl_pow_2(effb) - 4. * effa * effc)) / effa;
		y[4] = frequency;
		y[5] *= eff;
		y[6] *= eff;
		y[7] *= eff;
		return 0;
	}
	std::unique_ptr<Integrator> KerrTaubNUT::GetIntegrator(time_system time, coordinate_system coordinate, motion_mode motion) {
		if (motion != GEODESIC)
			PrintlnError("Motion not supported!");
		if (time == T) {
			if (coordinate == LAGRANGIAN)
				return std::make_unique<Integrator>(KerrTaubNutTLagrangianGeodesic, Jacobian, this);
			else if (coordinate == HAMILTONIAN)
				return std::make_unique<Integrator>(KerrTaubNutTHamiltonianGeodesic, Jacobian, this);
		} else if (time == TAU) {
			return std::make_unique<Integrator>(KerrTaubNutTauLagrangianGeodesic, Jacobian, this);
		}
		PrintlnError("Time system not supported!");
		return std::make_unique<Integrator>(
			[](double t, const double y[], double dydt[], void *params) -> int { return GSL_FAILURE; },
			[](double t, const double y[], double *dfdy, double dfdt[], void *params) -> int { return GSL_SUCCESS; });
	}

	int KerrTLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = reinterpret_cast<Kerr *>(params);
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = y[5]; // dr/dt
		dydt[2] = y[6]; // d\theta/dt
		dydt[3] = y[7]; // d\phi/dt
		const double r = y[1], r2 = gsl_pow_2(r), r4 = gsl_pow_4(r), a = kerr->a_, a2 = kerr->a2_, a4 = kerr->a4_, a2r2 = a2 + r2;
		const double sint = abs(sin(y[2])), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint), cost = GSL_SIGN(y[2]) * cos(y[2]), cost2 = gsl_pow_2(cost), sintcost = sint * cost, cott = cost / sint;
		const double Delta = r2 - 2. * r + a2, Delta_1 = 1. / Delta;
		const double rho2 = r2 + a2 * cost2, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2), r2a2cost2 = r2 - a2 * cost2;
		const double dydt4 = 2. * rho_4 * (Delta_1 * a2r2 * r2a2cost2 * y[5] - 2. * a2 * r * sintcost * y[6] * (1. - a * sint2 * y[7]) - Delta_1 * a * (2. * r4 + r2 * rho2 + a2 * r2a2cost2) * sint2 * y[5] * y[7]);
		// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = dydt4 * y[4];
		// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = (-Delta * r2a2cost2 * gsl_pow_2(1. - a * sint2 * y[7]) - (r * (a2 * sint2 - r) + a2 * cost2) * Delta_1 * rho4 * gsl_pow_2(y[5]) + 2. * a2 * sintcost * rho4 * y[5] * y[6] + r * Delta * rho4 * gsl_pow_2(y[6]) + Delta * sint2 * r * rho4 * gsl_pow_2(y[7])) * rho_6 + dydt4 * y[5];
		// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = (2. * a2 * r * sintcost - a2 * sintcost * rho4 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho4 * y[5] * y[6] - 4. * a * r * sintcost * a2r2 * y[7] + sintcost * (2. * a4 * r * sint4 + 4. * a2 * r * sint2 * rho2 + a2r2 * rho4) * gsl_pow_2(y[7])) * rho_6 + dydt4 * y[6];
		// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = (-2. * a * r2a2cost2 * Delta_1 * y[5] + 4. * a * r * cott * y[6] - 2. * Delta_1 * (r * rho4 - 2. * r2 * rho2 - r2a2cost2 * a2 * sint2) * y[5] * y[7] - 2. * cott * (rho4 + 2. * a2 * r * sint2) * y[6] * y[7]) * rho_4 + dydt4 * y[7];
		return GSL_SUCCESS;
	}
	int KerrTauLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = reinterpret_cast<Kerr *>(params);
		dydt[0] = y[4]; // dt/d\tau
		dydt[1] = y[5]; // dr/d\tau
		dydt[2] = y[6]; // d\theta/d\tau
		dydt[3] = y[7]; // d\phi/d\tau
		const double r = y[1], r2 = gsl_pow_2(r), a = kerr->a_, a2 = kerr->a2_, a4 = kerr->a4_, a2r2 = a2 + r2;
		const double sint = abs(sin(y[2])), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint), cost = GSL_SIGN(y[2]) * cos(y[2]), cost2 = gsl_pow_2(cost), sintcost = sint * cost, cott = cost / sint;
		const double Delta = r2 - 2. * r + a2, Delta_1 = 1. / Delta;
		const double rho2 = r2 + a2 * cost2, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2), r2a2cost2 = r2 - a2 * cost2;
		dydt[4] = -2. * rho_4 * (Delta_1 * a2r2 * r2a2cost2 * y[5] * (y[4] - a * sint2 * y[7]) - 2. * a2 * r * sintcost * y[6] * (y[4] - a * sint2 * y[7]) - 2. * Delta_1 * r2 * rho2 * a * sint2 * y[5] * y[7]);
		// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = (-Delta * r2a2cost2 * gsl_pow_2(y[4] - a * sint2 * y[7]) - (r * (a2 * sint2 - r) + a2 * cost2) * Delta_1 * rho4 * gsl_pow_2(y[5]) + 2. * a2 * sintcost * rho4 * y[5] * y[6] + r * Delta * rho4 * gsl_pow_2(y[6]) + Delta * sint2 * r * rho4 * gsl_pow_2(y[7])) * rho_6;
		// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = (2. * a2 * r * sintcost * gsl_pow_2(y[4]) - a2 * sintcost * rho4 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho4 * y[5] * y[6] - 4. * a * r * sintcost * a2r2 * y[4] * y[7] + sintcost * (2. * a4 * r * sint4 + 4. * a2 * r * sint2 * rho2 + a2r2 * rho4) * gsl_pow_2(y[7])) * rho_6;
		// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = (-2. * a * r2a2cost2 * Delta_1 * y[4] * y[5] + 4. * a * r * cott * y[4] * y[6] - 2. * Delta_1 * (r * rho4 - 2. * r2 * rho2 - r2a2cost2 * a2 * sint2) * y[5] * y[7] - 2. * cott * (rho4 + 2. * a2 * r * sint2) * y[6] * y[7]) * rho_4;
		return GSL_SUCCESS;
	}
	int KerrTHamiltonianGeodesic(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = reinterpret_cast<Kerr *>(params);
		const double r = y[1], r2 = gsl_pow_2(y[1]), a = kerr->a_, a2 = kerr->a2_, a2r2 = a2 + r2, pr2 = gsl_pow_2(y[5]), ptheta2 = gsl_pow_2(y[6]);
		const double E = 1. - y[4], E2 = gsl_pow_2(E), deltaE2 = (2. - y[4]) * y[4], L2 = gsl_pow_2(y[7]);
		const double sint = abs(sin(y[2])), sint2 = gsl_pow_2(sint), sint_2 = 1. / sint2, sint_4 = gsl_pow_2(sint_2), cost = GSL_SIGN(y[2]) * cos(y[2]), cost2 = gsl_pow_2(cost);
		const double Delta = a2r2 - 2. * r, Delta_1 = 1. / Delta, Delta_2 = gsl_pow_2(Delta_1);
		const double rho2 = r2 + a2 * cost2, rho_2 = 1. / rho2, rho_4 = gsl_pow_2(rho_2);
		const double Q = ptheta2 + cost2 * (a2 * deltaE2 + L2 * sint_2);
		const double R = -a2r2 * r2 * deltaE2 + E2 * 2. * r * a2 + 2. * r * r2 - Delta * Q - (r2 - 2. * r) * L2 - 4. * r * a * E * y[7]; // R = gsl_pow_2(E * a2r2 - a * y[7]) - Delta * (r2 + gsl_pow_2(y[7] - a * E) + Q);
		//[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
		dydt[0] = rho2 * Delta / (a2r2 * (E * a2r2 - a * y[7]) + a * Delta * (y[7] - a * E * sint2));		  // d\tau/dt
		dydt[1] = Delta * rho_2 * y[5] * dydt[0];															  // dr/dt
		dydt[2] = rho_2 * y[6] * dydt[0];																	  // d\theta/dt
		dydt[3] = (2. * E * a * r - y[7] * (a2 - Delta * (1. + cost2 * sint_2))) * Delta_1 * rho_2 * dydt[0]; // d\phi/dt
		dydt[4] = 0.;
		dydt[5] = ((a2 * cost2 + r * a2 * sint2 - r2) * pr2 + ((-2. * deltaE2 * r2 - a2 * deltaE2 + 3. * r - L2 - Q) * r + a2 * E2 + L2 - 2. * a * E * y[7] + Q) * rho2 * Delta_1 - ((r - 1.) * rho2 + Delta * r) * R * Delta_2) * rho_4 * dydt[0];
		dydt[6] = (-(Delta * pr2 - R * Delta_1) * a2 * rho_2 + a2 * deltaE2 + L2 * sint_2 + cost2 * sint_4 * L2) * sint * cost * rho_2 * dydt[0];
		dydt[7] = 0.;
		return GSL_SUCCESS;
	}
	int KerrTauHamiltonianGeodesic(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = reinterpret_cast<Kerr *>(params);
		const double r = y[1], r2 = gsl_pow_2(y[1]), a = kerr->a_, a2 = kerr->a2_, a2r2 = a2 + r2, pr2 = gsl_pow_2(y[5]), ptheta2 = gsl_pow_2(y[6]);
		const double E = 1. - y[4], E2 = gsl_pow_2(E), deltaE2 = (2. - y[4]) * y[4], L2 = gsl_pow_2(y[7]);
		const double sint = abs(sin(y[2])), sint2 = gsl_pow_2(sint), sint_2 = 1. / sint2, sint_4 = gsl_pow_2(sint_2), cost = GSL_SIGN(y[2]) * cos(y[2]), cost2 = gsl_pow_2(cost);
		const double Delta = a2r2 - 2. * r, Delta_1 = 1. / Delta, Delta_2 = gsl_pow_2(Delta_1);
		const double rho2 = r2 + a2 * cost2, rho_2 = 1. / rho2, rho_4 = gsl_pow_2(rho_2);
		const double Q = ptheta2 + cost2 * (a2 * deltaE2 + L2 * sint_2);
		const double R = -a2r2 * r2 * deltaE2 + E2 * 2. * r * a2 + 2. * r * r2 - Delta * Q - (r2 - 2. * r) * L2 - 4. * r * a * E * y[7]; // R = gsl_pow_2(E * a2r2 - a * y[7]) - Delta * (r2 + gsl_pow_2(y[7] - a * E) + Q);
		//[t,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
		dydt[0] = (a2r2 * (E * a2r2 - a * y[7]) * Delta_1 + a * (y[7] - a * E * sint2)) * rho_2;	// dt/d\tau
		dydt[1] = Delta * rho_2 * y[5];																// dr/d\tau
		dydt[2] = rho_2 * y[6];																		// d\theta/d\tau
		dydt[3] = (2. * E * a * r - y[7] * (a2 - Delta * (1. + cost2 * sint_2))) * Delta_1 * rho_2; // d\phi/d\tau
		dydt[4] = 0.;
		dydt[5] = ((a2 * cost2 + r * a2 * sint2 - r2) * pr2 + ((-2. * deltaE2 * r2 - a2 * deltaE2 + 3. * r - L2 - Q) * r + a2 * E2 + L2 - 2. * a * E * y[7] + Q) * rho2 * Delta_1 - ((r - 1.) * rho2 + Delta * r) * R * Delta_2) * rho_4;
		dydt[6] = (-(Delta * pr2 - R * Delta_1) * a2 * rho_2 + a2 * deltaE2 + L2 * sint_2 + cost2 * sint_4 * L2) * sint * cost * rho_2;
		dydt[7] = 0.;
		return GSL_SUCCESS;
	}
	int KerrTLagrangianHelical(double t, const double y[], double dydt[], void *params) {
		Kerr *kerr = reinterpret_cast<Kerr *>(params);
		const double g00 = 1. - 2. / y[1], sin_theta = abs(sin(y[2]));
		if (g00 <= 0)
			return GSL_FAILURE;
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = y[5]; // dr/dt
		dydt[2] = 0.;	// d\theta/dt = 0.
		dydt[3] = y[7]; // d\phi/dt
		// d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt
		dydt[4] = ((1. + gsl_pow_2(y[5] / g00)) / gsl_pow_2(y[1]) + y[1] * gsl_pow_2(sin_theta * y[7])) * y[5] / sqrt(g00 - (gsl_pow_2(y[5]) / g00 + gsl_pow_2(y[1] * sin_theta * y[7])));
		// d^2r/dt^2 = 0.
		dydt[5] = 0.;
		// d^2\theta/dt^2 = 0.
		dydt[6] = 0.;
		// d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
		dydt[7] = -2. * y[7] / y[1] * y[5];
		return GSL_SUCCESS;
	}
	int KerrTaubNutTLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		KerrTaubNUT *kerr_taub_nut = reinterpret_cast<KerrTaubNUT *>(params);
		dydt[0] = y[4]; // d\tau/dt
		dydt[1] = y[5]; // dr/dt
		dydt[2] = y[6]; // d\theta/dt
		dydt[3] = y[7]; // d\phi/dt
		const double r = y[1], r2 = gsl_pow_2(r), a = kerr_taub_nut->a_, a2 = kerr_taub_nut->a2_, l = kerr_taub_nut->l_, l2 = kerr_taub_nut->l2_;
		const double sint = abs(sin(y[2])), sint_1 = 1. / sint, sint2 = gsl_pow_2(sint), cost = GSL_SIGN(y[2]) * cos(y[2]), cost2 = gsl_pow_2(cost);
		const double Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
		const double lacost = l + a * cost, lacost2 = gsl_pow_2(lacost);
		const double rho2 = r2 + lacost2, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2);
		const double chi = a * sint2 - 2. * l * cost, rho2achi = r2 + l2 + a2;
		const double rho2rm = r * (r + 2. * l * lacost) - lacost2, rho2acost = (2. * r + l * lacost) * lacost - r2 * l;
		const double dydt4 = 2. * rho_4 * (Delta_1 * rho2achi * rho2rm * y[5] * (1. - chi * y[7]) - sint_1 * chi * rho2acost * y[6] * (1. - chi * y[7]) - 2. * Delta_1 * rho2 * r * ((r + l2) * a * sint2 + Delta * l * cost) * y[5] * y[7] - sint_1 * rho4 * l * (1. + cost2) * y[6] * y[7]);
		// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = dydt4 * y[4];
		// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = -Delta * rho2rm * rho_6 * gsl_pow_2(1. - chi * y[7]) + (rho2rm - r * a2 * sint2) * Delta_1 * rho_2 * gsl_pow_2(y[5]) + 2. * a * sint * lacost * rho_2 * y[5] * y[6] + r * Delta * rho_2 * gsl_pow_2(y[6]) + Delta * r * sint2 * rho_2 * gsl_pow_2(y[7]) + dydt4 * y[5];
		// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = rho2acost * sint * rho_6 * (a - 2. * rho2achi * y[7]) - a * sint * lacost * rho_2 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho_2 * y[5] * y[6] - (rho2 * (Delta * chi * lacost - gsl_pow_2(rho2achi) * cost) - 2. * ((r + l2) * a * sint2 + Delta * l * cost) * rho2achi * lacost) * sint * rho_6 * gsl_pow_2(y[7]) + dydt4 * y[6];
		// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = -2. * a * rho2rm * Delta_1 * rho_4 * y[5] * (1. - chi * y[7]) + 2. * rho2acost * rho_4 * sint_1 * y[6] * (1. - chi * y[7]) - 2. * (1. - a2 * sint2 * Delta_1) * r * rho_2 * y[5] * y[7] - 2. * cost * sint_1 * y[6] * y[7] + dydt4 * y[7];
		return GSL_SUCCESS;
	}
	int KerrTaubNutTauLagrangianGeodesic(double t, const double y[], double dydt[], void *params) {
		KerrTaubNUT *kerr_taub_nut = reinterpret_cast<KerrTaubNUT *>(params);
		dydt[0] = y[4]; // dt/d\tau
		dydt[1] = y[5]; // dr/d\tau
		dydt[2] = y[6]; // d\theta/d\tau
		dydt[3] = y[7]; // d\phi/d\tau
		const double r = y[1], r2 = gsl_pow_2(r), a = kerr_taub_nut->a_, a2 = kerr_taub_nut->a2_, l = kerr_taub_nut->l_, l2 = kerr_taub_nut->l2_;
		const double sint = abs(sin(y[2])), sint_1 = 1. / sint, sint2 = gsl_pow_2(sint), cost = GSL_SIGN(y[2]) * cos(y[2]), cost2 = gsl_pow_2(cost);
		const double Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
		const double lacost = l + a * cost, lacost2 = gsl_pow_2(lacost);
		const double rho2 = r2 + lacost2, rho_2 = 1. / rho2, rho4 = gsl_pow_2(rho2), rho_4 = gsl_pow_2(rho_2), rho_6 = gsl_pow_3(rho_2);
		const double chi = a * sint2 - 2. * l * cost, rho2achi = r2 + l2 + a2;
		const double rho2rm = r * (r + 2. * l * lacost) - lacost2, rho2acost = (2. * r + l * lacost) * lacost - r2 * l;
		// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
		dydt[4] = -2. * rho_4 * (Delta_1 * rho2achi * rho2rm * y[5] * (y[4] - chi * y[7]) - sint_1 * chi * rho2acost * y[6] * (y[4] - chi * y[7]) - 2. * Delta_1 * rho2 * r * ((r + l2) * a * sint2 + Delta * l * cost) * y[5] * y[7] - sint_1 * rho4 * l * (1. + cost2) * y[6] * y[7]);
		// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[5] = -Delta * rho2rm * rho_6 * gsl_pow_2(y[4] - chi * y[7]) + (rho2rm - r * a2 * sint2) * Delta_1 * rho_2 * gsl_pow_2(y[5]) + 2. * a * sint * lacost * rho_2 * y[5] * y[6] + r * Delta * rho_2 * gsl_pow_2(y[6]) + Delta * r * sint2 * rho_2 * gsl_pow_2(y[7]);
		// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[6] = rho2acost * sint * rho_6 * y[4] * (a * y[4] - 2. * rho2achi * y[7]) - a * sint * lacost * rho_2 * (Delta_1 * gsl_pow_2(y[5]) - gsl_pow_2(y[6])) - 2. * r * rho_2 * y[5] * y[6] - (rho2 * (Delta * chi * lacost - gsl_pow_2(rho2achi) * cost) - 2. * ((r + l2) * a * sint2 + Delta * l * cost) * rho2achi * lacost) * sint * rho_6 * gsl_pow_2(y[7]);
		// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
		dydt[7] = -2. * a * rho2rm * Delta_1 * rho_4 * y[5] * (y[4] - chi * y[7]) + 2. * rho2acost * rho_4 * sint_1 * y[6] * (y[4] - chi * y[7]) - 2. * (1. - a2 * sint2 * Delta_1) * r * rho_2 * y[5] * y[7] - 2. * cost * sint_1 * y[6] * y[7];
		return GSL_SUCCESS;
	}
	int KerrTaubNutTHamiltonianGeodesic(double t, const double y[], double dydt[], void *params) { // TODO:
		return GSL_FAILURE;
	}
	int Jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
		return GSL_SUCCESS;
	}
} // namespace SBody
