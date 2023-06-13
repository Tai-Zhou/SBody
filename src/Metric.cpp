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

#include "Utility.h"

namespace SBody {
	Metric::Metric(metric_mode mode, std::string name) : mode_(mode), name_(name) {}
	std::string &Metric::Name() {
		return name_;
	}
	int Metric::MetricTensor(const double position[], gsl_matrix *metric) {
		gsl_matrix_set_zero(metric);
		gsl_matrix_set(metric, 0, 0, -1.);
		gsl_matrix_set(metric, 1, 1, 1.);
		gsl_matrix_set(metric, 2, 2, gsl_pow_2(position[1]));
		gsl_matrix_set(metric, 3, 3, gsl_pow_2(position[1] * sin(position[2])));
		return GSL_SUCCESS;
	}
	double Metric::DotProduct(const double position[], const double x[], const double y[], const size_t dimension) {
		if (dimension == 3)
			return x[1] * y[1] + gsl_pow_2(position[1]) * x[2] * y[2] + gsl_pow_2(position[1] * sin(position[2])) * x[3] * y[3];
		return -x[0] * y[0] + x[1] * y[1] + gsl_pow_2(position[1]) * x[2] * y[2] + gsl_pow_2(position[1] * sin(position[2])) * x[3] * y[3];
	}
	double Metric::Distance(const double x[], const double y[], const size_t dimension) {
		if (dimension == 3)
			return gsl_pow_2(x[1] - y[1]) + gsl_pow_2(x[1] * (x[2] - y[2])) + gsl_pow_2(x[1] * sin(x[2]) * PhiDifference(x[3] - y[3]));
		return -gsl_pow_2(x[0] - y[0]) + gsl_pow_2(x[1] - y[1]) + gsl_pow_2(x[1] * (x[2] - y[2])) + gsl_pow_2(x[1] * sin(x[2]) * PhiDifference(x[3] - y[3]));
	}
	int Metric::LagrangianToHamiltonian(double y[]) {
		if (mode_ != HAMILTONIAN)
			return 1;
		return 1;
	}
	int Metric::HamiltonianToLagrangian(double y[]) {
		if (mode_ != HAMILTONIAN)
			return 1;
		return 1;
	}
	double Metric::Energy(const double y[]) {
		if (mode_ == T)
			return -y[4];
		else if (mode_ == TAU)
			return 1.; // TODO:
		else
			return 1. - y[4];
	}
	double Metric::AngularMomentum(const double y[]) {
		if (mode_ == T)
			return gsl_pow_2(y[1]) * y[7] / y[4];
		else if (mode_ == TAU)
			return 1.; // TODO:
		else
			return y[7];
	}
	double Metric::CarterConstant(const double y[], const double mu2) {
		if (mode_ == T)
			return gsl_pow_4(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(y[7] * cos(y[2]) * sin(y[2]))) / gsl_pow_2(y[4]);
		else if (mode_ == TAU)
			return 1.; // TODO:
		else
			return gsl_pow_2(y[6]) + gsl_pow_2(y[7] / tan(y[2]));
	}
	int Metric::NormalizeTimelikeGeodesic(double y[]) {
		const double g00 = 1. - 2. / y[1];
		if (g00 <= 0)
			return 1;
		y[4] = sqrt(g00 - (gsl_pow_2(y[5]) / g00 + gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
		return std::isnan(y[4]);
	}
	int Metric::NormalizeNullGeodesic(double y[], double frequency) {
		const double coefficient = GSL_SIGN(y[4]) / sqrt(gsl_pow_2(y[5]) + gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7]));
		y[4] = frequency;
		y[5] *= coefficient;
		y[6] *= coefficient;
		y[7] *= coefficient;
		return 0;
	}
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
	Integrator Metric::GetIntegrator(int coordinate) {
		return Integrator(metric::Schwarzschild::function, metric::Schwarzschild::jacobian, coordinate);
	}

	Newton::Newton(int PN, metric_mode mode, std::string name) : Metric(mode, name), PN_(PN){};
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
			return x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
		return -x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
	}
	double Newton::Distance(const double x[], const double y[], const size_t dimension) {
		if (dimension == 3)
			return gsl_pow_2(x[1] - y[1]) + gsl_pow_2(x[2] - y[2]) + gsl_pow_2(x[3] - y[3]);
		return -gsl_pow_2(x[0] - y[0]) + gsl_pow_2(x[1] - y[1]) + gsl_pow_2(x[2] - y[2]) + gsl_pow_2(x[3] - y[3]);
	}
	int Newton::LagrangianToHamiltonian(double y[]) {
		if (mode_ != HAMILTONIAN)
			return 1;
		return 1;
	}
	int Newton::HamiltonianToLagrangian(double y[]) {
		if (mode_ != HAMILTONIAN)
			return 1;
		return 1;
	}
	double Newton::Energy(const double y[]) {
		const double r = Norm(y + 1), r_1 = 1. / r, v2 = SBody::Dot(y + 5);
		const double m_r = r_1, rdot = SBody::Dot(y + 1, y + 5) * r_1, v4 = gsl_pow_2(v2), v6 = gsl_pow_3(v2), v8 = gsl_pow_2(v4);
		const double m_r2 = gsl_pow_2(m_r), m_r3 = gsl_pow_3(m_r), m_r4 = gsl_pow_4(m_r), rdot2 = gsl_pow_2(rdot);
		double E = 0.5 * v2 - m_r;
		if (PN_ & 1)
			E += 0.5 * m_r2 + 0.375 * v4 + 1.5 * v2 * m_r;
		if (PN_ & 2)
			E += -0.5 * m_r3 + 0.3125 * v6 + 1.75 * m_r2 * v2 + 0.5 * m_r2 * rdot2 + 2.625 * m_r * v4;
		if (PN_ & 8)
			E += 0.375 * m_r4 + 1.25 * m_r3 * v2 + 1.5 * m_r3 * rdot2 + 0.2734375 * v8 + 8.4375 * m_r2 * v4 + 0.75 * m_r2 * v2 * rdot2 + 3.4375 * m_r * v6;
		return E;
	}
	double Newton::AngularMomentum(const double y[]) {
		const double r = Norm(y + 1), r_1 = 1. / r, v2 = SBody::Dot(y + 5);
		const double m_r = r_1, rdot = SBody::Dot(y + 1, y + 5) * r_1, v4 = gsl_pow_2(v2), v6 = gsl_pow_3(v2);
		const double m_r2 = gsl_pow_2(m_r), m_r3 = gsl_pow_3(m_r), rdot2 = gsl_pow_2(rdot);
		double J[3], eff = 0;
		Cross(y, y + 3, J);
		if (PN_ & 1)
			eff += 3. * m_r + 0.5 * v2;
		if (PN_ & 2)
			eff += 3.5 * m_r2 + 0.375 * v4 + 3.5 * m_r * v2;
		if (PN_ & 8)
			eff += 2.5 * m_r3 + 0.3125 * v6 + 11.25 * m_r2 * v2 + 0.5 * m_r2 * rdot2 + 4.125 * m_r * v4;
		for (int i = 0; i < 3; ++i)
			J[i] += J[i] * eff;
		return Norm(J);
	}
	double Newton::CarterConstant(const double y[], const double mu2) { // FIXME: not verified!
		double c[3];
		Cross(y + 1, y + 5, c);
		return SBody::Dot(c) / gsl_pow_2(y[4]);
	}
	int Newton::NormalizeTimelikeGeodesic(double y[]) { // TODO: limit the light speed
		y[4] = 1;
		if (SBody::Dot(y + 5) >= 1)
			return 1;
		return 0;
	}
	int Newton::NormalizeNullGeodesic(double y[], double frequency) {
		y[4] = 1;
		const double v_1 = GSL_SIGN(y[4]) / SBody::Norm(y + 5);
		for (int i = 5; i < 8; ++i)
			y[i] *= v_1;
		return 0;
	}
	Integrator Newton::GetIntegrator(int coordinate) {
		Integrator integrator(metric::Newton::function, metric::Newton::jacobian, coordinate, this);
		return integrator;
	}

	Schwarzschild::Schwarzschild(metric_mode mode, std::string name) : Metric(mode, name) {}
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
	double Schwarzschild::Distance(const double x[], const double y[], const size_t dimension) {
		if (dimension == 3)
			return x[1] * gsl_pow_2(x[1] - y[1]) / (x[1] - 2.) + gsl_pow_2(x[1] * (x[2] - y[2])) + gsl_pow_2(x[1] * sin(x[2]) * PhiDifference(x[3] - y[3]));
		return -(1. - 2. / x[1]) * gsl_pow_2(x[0] - y[0]) + x[1] * gsl_pow_2(x[1] - y[1]) / (x[1] - 2.) + gsl_pow_2(x[1] * (x[2] - y[2])) + gsl_pow_2(x[1] * sin(x[2]) * PhiDifference(x[3] - y[3]));
	}
	int Schwarzschild::LagrangianToHamiltonian(double y[]) {
		if (mode_ != HAMILTONIAN)
			return 1;
		const double dtdtau = 1. / y[4], r2 = gsl_pow_2(y[1]);
		y[4] = (y[4] - 1. + 2. / y[1]) * dtdtau; // 1 + p_t
		y[5] *= y[1] / (y[1] - 2.) * dtdtau;
		y[6] *= r2 * dtdtau;
		y[7] *= r2 * gsl_pow_2(sin(y[2])) * dtdtau;
		return 0;
	}
	int Schwarzschild::HamiltonianToLagrangian(double y[]) {
		if (mode_ != HAMILTONIAN)
			return 1;
		const double r_1 = 1. / y[1], g00 = 1. - 2. * r_1;
		y[4] = -g00 / (y[4] - 1.);
		y[5] *= g00 * y[4];
		y[6] *= gsl_pow_2(r_1) * y[4];
		y[7] *= gsl_pow_2(r_1 / sin(y[2])) * y[4];
		return 0;
	}
	double Schwarzschild::Energy(const double y[]) {
		if (mode_ == T)
			return (y[1] - 2.) / (y[1] * y[4]);
		else if (mode_ == TAU)
			return 1.; // TODO:
		else
			return 1. - y[4];
	}
	double Schwarzschild::AngularMomentum(const double y[]) {
		if (mode_ == T)
			return gsl_pow_2(y[1]) * y[7] / y[4];
		else if (mode_ == TAU)
			return 1.; // TODO:
		else
			return y[7];
	}
	double Schwarzschild::CarterConstant(const double y[], const double mu2) {
		if (mode_ == T)
			return gsl_pow_4(y[1]) * (gsl_pow_2(y[6]) + gsl_pow_2(y[7] * cos(y[2]) * sin(y[2]))) / gsl_pow_2(y[4]);
		else if (mode_ == TAU)
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
		const double coefficient = GSL_SIGN(y[4]) * g00 / sqrt(gsl_pow_2(y[5]) + g00 * (gsl_pow_2(y[1] * y[6]) + gsl_pow_2(y[1] * sin(y[2]) * y[7])));
		y[4] = frequency;
		y[5] *= coefficient;
		y[6] *= coefficient;
		y[7] *= coefficient;
		return 0;
	}
	Integrator Schwarzschild::GetIntegrator(int coordinate) {
		switch (mode_) {
		case T:
			return Integrator(metric::Schwarzschild::function, metric::Schwarzschild::jacobian, coordinate);
		case TAU:
			return Integrator(metric::Schwarzschild::functionTau, metric::Schwarzschild::jacobianTau, coordinate);
		case HAMILTONIAN:
			return Integrator(metric::Schwarzschild::functionHamiltonian, metric::Schwarzschild::jacobianHamiltonian, coordinate);
		case RIAF:
			return Integrator(metric::Schwarzschild::functionRIAF, metric::Schwarzschild::jacobianRIAF, coordinate);
		case HELICAL:
			return Integrator(metric::Schwarzschild::functionHelicalWithFixedRadialSpeed, metric::Schwarzschild::jacobianHelicalWithFixedRadialSpeed, coordinate);
		}
	}

	Kerr::Kerr(double spin, metric_mode mode, std::string name) : Schwarzschild(mode, "Kerr"), a_(spin), a2_(a_ * a_), a4_(a2_ * a2_) {}
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
	double Kerr::Distance(const double x[], const double y[], const size_t dimension) {
		const double r = x[1], r2 = gsl_pow_2(r), d0 = x[0] - y[0], d3 = PhiDifference(x[3] - y[3]);
		const double sint = sin(x[2]), sint2 = gsl_pow_2(sint), cost = cos(x[2]), cost2 = gsl_pow_2(cost);
		const double Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * cost2;
		const double mr_rho2 = 2. * r / rho2;
		if (dimension == 3)
			return (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + ((r2 + a2_) * sint2 + 2. * r * a2_ * gsl_pow_2(sint2) / rho2) * gsl_pow_2(d3);
		return (mr_rho2 - 1.) * gsl_pow_2(d0) - 2. * mr_rho2 * a_ * sint2 * d0 * d3 + (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + ((r2 + a2_) * sint2 + mr_rho2 * a2_ * gsl_pow_2(sint2)) * gsl_pow_2(d3);
	}
	int Kerr::LagrangianToHamiltonian(double y[]) {
		if (mode_ != HAMILTONIAN)
			return 1;
		const double dtdtau = 1. / y[4], r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2]));
		const double Delta = r2 - 2. * y[1] + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double mr_rho2 = 2. * y[1] / rho2;
		y[4] = (y[4] - 1. + mr_rho2 * (1. - a_ * sint2 * y[7])) * dtdtau; // 1 + p_t
		y[5] *= rho2 / Delta * dtdtau;
		y[6] *= rho2 * dtdtau;
		y[7] = (-mr_rho2 * a_ + (a2_ + r2 + mr_rho2 * a2_ * sint2) * y[7]) * sint2 * dtdtau;
		return 0;
	}
	int Kerr::HamiltonianToLagrangian(double y[]) {
		if (mode_ != HAMILTONIAN)
			return 1;
		const double pt = y[4] - 1., r2 = gsl_pow_2(y[1]);
		const double Delta = r2 - 2. * y[1] + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double rho_2 = 1. / rho2, mr_rho2 = 2. * y[1] * rho_2;
		y[4] = -Delta / ((Delta + mr_rho2 * (a2_ + r2)) * pt + mr_rho2 * a_ * y[7]);
		y[5] *= Delta * y[4] * rho_2;
		y[6] *= y[4] * rho_2;
		y[7] = (-mr_rho2 * a_ * pt + (1. - mr_rho2) / gsl_pow_2(sin(y[2])) * y[7]) / Delta * y[4];
		return 0;
	}
	double Kerr::Energy(const double y[]) {
		if (mode_ == T)
			return (1. - 2. * y[1] / (gsl_pow_2(y[1]) + a2_ * gsl_pow_2(cos(y[2]))) * (1. - a_ * gsl_pow_2(sin(y[2])) * y[7])) / y[4];
		else if (mode_ == TAU)
			return 1.; // TODO:
		else
			return 1. - y[4];
	}
	double Kerr::AngularMomentum(const double y[]) {
		const double r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2]));
		const double mr_rho2 = 2. * y[1] / (r2 + a2_ * gsl_pow_2(cos(y[2])));
		if (mode_ == T)
			return (-mr_rho2 * a_ + (a2_ + r2 + mr_rho2 * a2_ * sint2) * y[7]) * sint2 / y[4];
		else if (mode_ == TAU)
			return 1.; // TODO:
		else
			return y[7];
	}
	double Kerr::CarterConstant(const double y[], const double mu2) {
		const double r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2])), cost2 = gsl_pow_2(cos(y[2]));
		const double rho2 = r2 + a2_ * cost2;
		const double mr_rho2 = 2. * y[1] / rho2;
		if (mode_ == T)
			return mu2 * cost2 * a2_ + (gsl_pow_2(rho2 * y[6]) + cost2 * (gsl_pow_2(-mr_rho2 * a_ + (a2_ + r2 + mr_rho2 * a2_ * sint2) * y[7]) * sint2 - a2_ * gsl_pow_2(mr_rho2 * (1. - a_ * sint2 * y[7]) - 1.))) / gsl_pow_2(y[4]);
		else if (mode_ == TAU)
			return 1.; // TODO:
		else
			return gsl_pow_2(y[6]) + gsl_pow_2(cos(y[2])) * (a2_ * (mu2 - gsl_pow_2(1. - y[4])) + gsl_pow_2(y[7] / sin(y[2])));
	}
	int Kerr::NormalizeTimelikeGeodesic(double y[]) {
		const double r = y[1], r2 = gsl_pow_2(r), a2r2 = a2_ + r2;
		const double sint = sin(y[2]), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint);
		const double Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double mr_rho2 = 2. * r / rho2;
		y[7] += 2. * a_ * r / (gsl_pow_2(a2r2) - a2_ * Delta * sint2);
		y[4] = sqrt(1. - mr_rho2 + 2. * mr_rho2 * a_ * sint2 * y[7] - (rho2 / (r2 - 2. * r + a2_) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + (a2r2 * sint2 + mr_rho2 * a2_ * sint4) * gsl_pow_2(y[7])));
		return std::isnan(y[4]);
	}
	int Kerr::NormalizeNullGeodesic(double y[], double frequency) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sint = sin(y[2]), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint);
		const double rho2 = r2 + a2_ * gsl_pow_2(cos(y[2]));
		const double mr_rho2 = 2. * r / rho2;
		const double effa = rho2 / (r2 - 2. * r + a2_) * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + ((a2_ + r2) * sint2 + mr_rho2 * a2_ * sint4) * gsl_pow_2(y[7]);
		const double effb = -2. * mr_rho2 * a_ * sint2 * y[7];
		const double effc = mr_rho2 - 1.;
		const double eff = GSL_SIGN(y[4]) * 0.5 * (-effb + sqrt(gsl_pow_2(effb) - 4. * effa * effc)) / effa;
		y[4] = frequency;
		y[5] *= eff;
		y[6] *= eff;
		y[7] *= eff;
		return 0;
	}
	Integrator Kerr::GetIntegrator(int coordinate) {
		if (mode_ == T)
			return Integrator(metric::Kerr::function, metric::Kerr::jacobian, coordinate, this);
		else if (mode_ == TAU)
			return Integrator(metric::Kerr::functionTau, metric::Kerr::jacobianTau, coordinate, this);
		else
			return Integrator(metric::Kerr::functionHamiltonian, metric::Kerr::jacobianHamiltonian, coordinate, this);
	}

	KerrTaubNUT::KerrTaubNUT(double spin, double charge, double NUT, metric_mode mode, std::string name) : Kerr(spin, mode, "Kerr-Taub-NUT"), e_(charge), e2_(e_ * e_), e4_(e2_ * e2_), l_(NUT), l2_(l_ * l_), l4_(l2_ * l2_) {}
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
	double KerrTaubNUT::Distance(const double x[], const double y[], const size_t dimension) {
		const double r = x[1], sint = sin(x[2]), cost = cos(x[2]), d0 = x[0] - y[0], d3 = PhiDifference(x[3] - y[3]);
		const double r2 = gsl_pow_2(r), sint2 = gsl_pow_2(sint);
		const double Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + gsl_pow_2(l_ + a_ * cost), chi = a_ * sint2 - 2. * l_ * cost;
		const double rho_2 = 1. / rho2;
		if (dimension == 3)
			return (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + (gsl_pow_2(rho2 + a_ * chi) * sint2 - gsl_pow_2(chi) * Delta) * rho_2 * gsl_pow_2(d3);
		return (a2_ * sint2 - Delta) * rho_2 * gsl_pow_2(d0) - 4. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) * rho_2 * d0 * d3 + (rho2 / Delta) * gsl_pow_2(x[1] - y[1]) + rho2 * gsl_pow_2(x[2] - y[2]) + (gsl_pow_2(rho2 + a_ * chi) * sint2 - gsl_pow_2(chi) * Delta) * rho_2 * gsl_pow_2(d3);
	}
	int KerrTaubNUT::LagrangianToHamiltonian(double y[]) {
		if (mode_ != HAMILTONIAN)
			return 1;
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
	int KerrTaubNUT::HamiltonianToLagrangian(double y[]) {
		if (mode_ != HAMILTONIAN)
			return 1;
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
	double KerrTaubNUT::Energy(const double y[]) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sint2 = gsl_pow_2(sin(y[2])), cost = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		if (mode_ == T)
			return (Delta - a2_ * sint2 + 2. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) * y[7]) / ((r2 + gsl_pow_2(l_ + a_ * cost)) * y[4]);
		else if (mode_ == TAU)
			return 1.; // TODO:
		else
			return 1. - y[4];
	}
	double KerrTaubNUT::AngularMomentum(const double y[]) {
		const double r = y[1], r2 = gsl_pow_2(r);
		const double sint2 = gsl_pow_2(sin(y[2])), cost = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		if (mode_ == T)
			return (-2. * ((r + l2_) * a_ * sint2 + Delta * l_ * cost) + (gsl_pow_2(r2 + l2_ + a2_) * sint2 - gsl_pow_2(a_ * sint2 - 2. * l_ * cost) * Delta) * y[7]) / ((r2 + gsl_pow_2(l_ + a_ * cost)) * y[4]);
		else if (mode_ == TAU)
			return 1.; // TODO:
		else
			return y[7];
	}
	double KerrTaubNUT::CarterConstant(const double y[], const double mu2) { // TODO:
		const double r2 = gsl_pow_2(y[1]), sint2 = gsl_pow_2(sin(y[2])), cost2 = gsl_pow_2(cos(y[2]));
		const double rho2 = r2 + a2_ * cost2;
		const double mr_rho2 = 2. * y[1] / rho2;
		if (mode_ == T)
			return mu2 * cost2 * a2_ + (gsl_pow_2(rho2 * y[6]) + cost2 * (gsl_pow_2(-mr_rho2 * a_ + (a2_ + r2 + mr_rho2 * a2_ * sint2) * y[7]) * sint2 - a2_ * gsl_pow_2(mr_rho2 * (1. - a_ * sint2 * y[7]) - 1.))) / gsl_pow_2(y[4]);
		else if (mode_ == TAU)
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
		const double sint = abs(sin(y[2])), sint2 = gsl_pow_2(sint), cost = GSL_SIGN(y[2]) * cos(y[2]);
		const double Delta = r2 - 2. * r - l2_ + a2_;
		const double lacost = l_ + a_ * cost;
		const double rho2 = r2 + gsl_pow_2(lacost), rho_2 = 1. / rho2;
		const double chi = a_ * sint2 - 2. * l_ * cost, rho2achi = r2 + l2_ + a2_;
		const double effa = rho2 / Delta * gsl_pow_2(y[5]) + rho2 * gsl_pow_2(y[6]) + rho_2 * (gsl_pow_2(rho2achi) * sint2 - gsl_pow_2(chi) * Delta) * gsl_pow_2(y[7]);
		const double effb = -4. * rho_2 * ((r + l2_) * chi + l_ * cost * rho2achi) * y[7];
		const double effc = -rho_2 * (Delta - a2_ * sint2);
		const double eff = GSL_SIGN(y[4]) * 0.5 * (-effb + sqrt(gsl_pow_2(effb) - 4. * effa * effc)) / effa;
		y[4] = frequency;
		y[5] *= eff;
		y[6] *= eff;
		y[7] *= eff;
		return 0;
	}
	Integrator KerrTaubNUT::GetIntegrator(int coordinate) {
		if (mode_ == T)
			return Integrator(metric::KerrTaubNUT::function, metric::KerrTaubNUT::jacobian, coordinate, this);
		else if (mode_ == TAU)
			return Integrator(metric::KerrTaubNUT::functionTau, metric::KerrTaubNUT::jacobianTau, coordinate, this);
		else
			return Integrator(metric::KerrTaubNUT::functionHamiltonian, metric::KerrTaubNUT::jacobianHamiltonian, coordinate, this);
	}

	namespace metric {
		namespace Newton {
			int function(double t, const double y[], double dydt[], void *params) { // FIXME: rewrite in spherical coordinates
				class Newton *newton = reinterpret_cast<class Newton *>(params);
				const double r = Norm(y + 1), r_1 = 1. / r, v2 = SBody::Dot(y + 5);
				const double m_r = r_1, rdot = SBody::Dot(y + 1, y + 5) * r_1;
				const double F = m_r * gsl_pow_2(r_1), m_r2 = gsl_pow_2(m_r), rdot2 = gsl_pow_2(rdot);
				dydt[0] = y[4];
				dydt[1] = y[5];
				dydt[2] = y[6];
				dydt[3] = y[7];
				dydt[4] = 0.;
				double A = 0, B = 0;
				if (newton->PN_ & 1) {
					A += (v2 - 4. * m_r);
					B += (-4. * rdot);
					if (params != nullptr) {
						A *= *(static_cast<double *>(params));
						B *= *(static_cast<double *>(params));
					}
				}
				if (newton->PN_ & 2) {
					A += m_r * (-2. * rdot2 + 9. * m_r);
					B += 2. * m_r * rdot;
				}
				if (newton->PN_ & 8) {
					A += -16. * gsl_pow_3(m_r) + rdot2 * m_r2;
					B += -4. * rdot * m_r2;
				}
				dydt[5] = -F * (y[1] + A * y[1] + B * y[5] * r);
				dydt[6] = -F * (y[2] + A * y[2] + B * y[6] * r);
				dydt[7] = -F * (y[3] + A * y[3] + B * y[7] * r);
				return GSL_SUCCESS;
			}
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
		} // namespace Newton
		namespace Schwarzschild {
			int function(double t, const double y[], double dydt[], void *params) {
				dydt[0] = y[4]; // d\tau/dt
				dydt[1] = y[5]; // dr/dt
				dydt[2] = y[6]; // d\theta/dt
				dydt[3] = y[7]; // d\phi/dt
				const double r = y[1], sint = sin(y[2]), cost = cos(y[2]);
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
			}
			int functionTau(double t, const double y[], double dydt[], void *params) {
				dydt[0] = y[4]; // dt/d\tau
				dydt[1] = y[5]; // dr/d\tau
				dydt[2] = y[6]; // d\theta/d\tau
				dydt[3] = y[7]; // d\phi/d\tau
				const double r = y[1], sint = sin(y[2]), cost = cos(y[2]);
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
			}
			int functionHamiltonian(double t, const double y[], double dydt[], void *params) {
				const double r_1 = 1. / y[1], r_2 = gsl_pow_2(r_1), g00 = 1. - 2. * r_1, E = 1. - y[4], L2 = gsl_pow_2(y[7]);
				const double sint_1 = 1. / sin(y[2]), sint_2 = gsl_pow_2(sint_1);
				//[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
				dydt[0] = g00 / E;						 // d\tau/dt
				dydt[1] = g00 * y[5] * dydt[0];			 // dr/dt
				dydt[2] = y[6] * r_2 * dydt[0];			 // d\theta/dt
				dydt[3] = y[7] * sint_2 * r_2 * dydt[0]; // d\phi/dt
				dydt[4] = 0.;
				dydt[5] = (-(gsl_pow_2(y[5]) + gsl_pow_2(E) / gsl_pow_2(g00)) + (gsl_pow_2(y[6]) + L2 * sint_2) * r_1) * r_2 * dydt[0];
				dydt[6] = sint_2 * L2 * cos(y[2]) * sint_1 * r_2 * dydt[0];
				dydt[7] = 0.;
				return GSL_SUCCESS;
			}
			int functionRIAF(double t, const double y[], double dydt[], void *params) {
				return GSL_SUCCESS;
			}
			int functionHelicalWithFixedRadialSpeed(double t, const double y[], double dydt[], void *params) {
				const double g00 = 1. - 2. / y[1], sin_theta = sin(y[2]);
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
			int functionHelicalWithFixedRadialMomentum(double t, const double y[], double dydt[], void *params) {
				// WIP
				const double g00 = 1. - 2. / y[1];
				if (g00 <= 0)
					return GSL_FAILURE;
				dydt[0] = sqrt(g00 - (gsl_pow_2(y[5]) / g00 + gsl_pow_2(y[1] * sin(y[2]) * y[7]))); // d\tau/dt
				dydt[1] = y[5];																		// dr/dt
				dydt[2] = 0.;																		// d\theta/dt = 0.
				dydt[3] = y[7];																		// d\phi/dt
				dydt[4] = 1.;
				// d^2r/dt^2 = 0.
				// g00 * dr/d\tau = C
				// dr/d\tau = C / g00
				// dr/dt = C / g00 * d\tau/dt
				dydt[5] = 1.;
				dydt[6] = 0.;
				dydt[7] = 1.;
				return GSL_SUCCESS;
			}
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianTau(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianRIAF(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianHelicalWithFixedRadialSpeed(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
		} // namespace Schwarzschild
		namespace Kerr {
			int function(double t, const double y[], double dydt[], void *params) {
				class Kerr *kerr = reinterpret_cast<class Kerr *>(params);
				dydt[0] = y[4]; // d\tau/dt
				dydt[1] = y[5]; // dr/dt
				dydt[2] = y[6]; // d\theta/dt
				dydt[3] = y[7]; // d\phi/dt
				const double r = y[1], r2 = gsl_pow_2(r), r4 = gsl_pow_4(r), a = kerr->a_, a2 = kerr->a2_, a4 = kerr->a4_, a2r2 = a2 + r2;
				const double sint = sin(y[2]), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint), cost = cos(y[2]), cost2 = gsl_pow_2(cost), sintcost = sint * cost, cott = cost / sint;
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
			int functionTau(double t, const double y[], double dydt[], void *params) {
				class Kerr *kerr = reinterpret_cast<class Kerr *>(params);
				dydt[0] = y[4]; // dt/d\tau
				dydt[1] = y[5]; // dr/d\tau
				dydt[2] = y[6]; // d\theta/d\tau
				dydt[3] = y[7]; // d\phi/d\tau
				const double r = y[1], r2 = gsl_pow_2(r), a = kerr->a_, a2 = kerr->a2_, a4 = kerr->a4_, a2r2 = a2 + r2;
				const double sint = sin(y[2]), sint2 = gsl_pow_2(sint), sint4 = gsl_pow_4(sint), cost = cos(y[2]), cost2 = gsl_pow_2(cost), sintcost = sint * cost, cott = cost / sint;
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
			int functionHamiltonian(double t, const double y[], double dydt[], void *params) {
				class Kerr *kerr = reinterpret_cast<class Kerr *>(params);
				const double r = y[1], r2 = gsl_pow_2(y[1]), a = kerr->a_, a2 = kerr->a2_, a2r2 = a2 + r2, pr2 = gsl_pow_2(y[5]), ptheta2 = gsl_pow_2(y[6]);
				const double E = 1. - y[4], E2 = gsl_pow_2(E), deltaE2 = (2. - y[4]) * y[4], L2 = gsl_pow_2(y[7]);
				const double sint = sin(y[2]), sint2 = gsl_pow_2(sint), sint_2 = 1. / sint2, sint_4 = gsl_pow_2(sint_2), cost = cos(y[2]), cost2 = gsl_pow_2(cost);
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
			int functionHamiltonianTau(double t, const double y[], double dydt[], void *params) {
				class Kerr *kerr = reinterpret_cast<class Kerr *>(params);
				const double r = y[1], r2 = gsl_pow_2(y[1]), a = kerr->a_, a2 = kerr->a2_, a2r2 = a2 + r2, pr2 = gsl_pow_2(y[5]), ptheta2 = gsl_pow_2(y[6]);
				const double E = 1. - y[4], E2 = gsl_pow_2(E), deltaE2 = (2. - y[4]) * y[4], L2 = gsl_pow_2(y[7]);
				const double sint = sin(y[2]), sint2 = gsl_pow_2(sint), sint_2 = 1. / sint2, sint_4 = gsl_pow_2(sint_2), cost = cos(y[2]), cost2 = gsl_pow_2(cost);
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
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianTau(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
		} // namespace Kerr
		namespace KerrTaubNUT {
			int function(double t, const double y[], double dydt[], void *params) {
				class KerrTaubNUT *kerrtaubnut = reinterpret_cast<class KerrTaubNUT *>(params);
				dydt[0] = y[4]; // d\tau/dt
				dydt[1] = y[5]; // dr/dt
				dydt[2] = y[6]; // d\theta/dt
				dydt[3] = y[7]; // d\phi/dt
				const double r = y[1], r2 = gsl_pow_2(r), a = kerrtaubnut->a_, a2 = kerrtaubnut->a2_, l = kerrtaubnut->l_, l2 = kerrtaubnut->l2_;
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
			int functionTau(double t, const double y[], double dydt[], void *params) {
				class KerrTaubNUT *kerrtaubnut = reinterpret_cast<class KerrTaubNUT *>(params);
				dydt[0] = y[4]; // dt/d\tau
				dydt[1] = y[5]; // dr/d\tau
				dydt[2] = y[6]; // d\theta/d\tau
				dydt[3] = y[7]; // d\phi/d\tau
				const double r = y[1], r2 = gsl_pow_2(r), a = kerrtaubnut->a_, a2 = kerrtaubnut->a2_, l = kerrtaubnut->l_, l2 = kerrtaubnut->l2_;
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
			int functionHamiltonian(double t, const double y[], double dydt[], void *params) { // TODO:
				return GSL_SUCCESS;
			}
			int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianTau(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
			int jacobianHamiltonian(double t, const double y[], double *dfdy, double dfdt[], void *params) {
				return GSL_SUCCESS;
			}
		} // namespace KerrTaubNUT
	}	  // namespace metric
} // namespace SBody
