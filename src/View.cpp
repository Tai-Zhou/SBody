/**
 * @file View.cpp
 * @author Tai Zhou
 * @brief This is the file setting the view of a distant observer.
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "View.h"

#include <cmath>
#include <memory>

#include <fmt/core.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "IO.h"
#include "Metric.h"
#include "Object.h"
#include "Unit.h"
#include "Utility.h"

using namespace std;

namespace SBody {
	View::View(shared_ptr<Metric> metric, double r, double theta, string file_name) : metric_(move(metric)), r_(r), theta_(theta), sin_theta_observer_(sin(theta)), cos_theta_observer_(cos(theta)), t_final_(-r - 2e4 * 1.), output_(make_unique<NumPy>(file_name, vector<int>({9}))) {
		bar_ = make_unique<indicators::BlockProgressBar>(indicators::option::ShowElapsedTime{true},
														 indicators::option::ShowRemainingTime{true},
														 indicators::option::ForegroundColor{indicators::Color(4)},
														 indicators::option::FontStyles{vector<indicators::FontStyle>{indicators::FontStyle::bold}});
	}
	int View::InitializePhoton(double *photon, double alpha, double beta) {
		photon[0] = 0.;
		photon[1] = r_;
		photon[5] = 1.;
		photon[8] = 0.;
		photon[9] = 0.;
		if (sin_theta_observer_ < epsilon) {
			const double k = gsl_hypot(alpha, beta);
			if (theta_ < M_PI_2) {
				photon[2] = 1e-15;
				photon[3] = ModBy2Pi(atan2(alpha, -beta));
				photon[6] = -k / gsl_pow_2(r_);
			} else {
				photon[2] = M_PI - 1e-15;
				photon[3] = ModBy2Pi(atan2(alpha, beta));
				photon[6] = k / gsl_pow_2(r_);
			}
			photon[7] = 0.;
		} else {
			photon[2] = theta_;
			photon[3] = 0.,
			photon[6] = beta / gsl_pow_2(r_);
			photon[7] = -alpha / (gsl_pow_2(r_) * sin_theta_observer_);
		}
		return metric_->NormalizeNullGeodesic(photon, 1.);
	}
	int View::TraceStar(Star &star, int ray_number) { // FIXME:!!!!
		double star_position[8];
		star.Position(star_position);
		const double r_star = star_position[1], sin_theta_star = sin(star_position[2]), cos_theta_star = cos(star_position[2]), sin_phi_star = sin(star_position[3]), cos_phi_star = cos(star_position[3]);
		if (r_star <= 3.)
			fmt::print(stderr, "\033[103m[WRN]\033[0m star radius = {:.6f}\n", r_star);
		const double alpha_coefficient = sin_theta_star * sin_phi_star, beta_coefficient = cos_theta_star * sin_theta_observer_ - sin_theta_star * cos_phi_star * cos_theta_observer_, iteration_coefficient = tanh(0.05 * r_star), cos_observer_star = sin_theta_observer_ * sin_theta_star * cos_phi_star + cos_theta_observer_ * cos_theta_star, sin_observer_star = sqrt(gsl_pow_2(alpha_coefficient) + gsl_pow_2(beta_coefficient)), theta_observer_star = acos(cos_observer_star);
		double alpha0 = GSL_POSINF, alpha1 = 0, beta0 = GSL_POSINF, beta1 = 0, cosph, h;
		if (cos_observer_star <= -1.) {
			fmt::print(stderr, "\033[103m[WRN]\033[0m star behind black hole, cos theta = {:.6f}\n", cos_observer_star);
			alpha1 = 2. * sqrt(r_star);
			beta1 = 0.;
		} else if (cos_observer_star < 1.) {
			if (theta_observer_star < M_PI_2) {
				const double effective_radius = r_star + gsl_pow_3(theta_observer_star / M_PI_2) / sin_observer_star;
				alpha1 = effective_radius * alpha_coefficient;
				beta1 = effective_radius * beta_coefficient;
			} else {
				const double effective_radius = r_star + 0.5 * (1. + sqrt(gsl_pow_2(r_star * sin_observer_star + 1) - 16. * r_star * cos_observer_star)) / sin_observer_star;
				// b=(r*sin(theta)+1+sqrt(gsl_pow_2(r*sin(theta)+1)-16*cos(theta)*r))/2
				alpha1 = effective_radius * alpha_coefficient;
				beta1 = effective_radius * beta_coefficient;
			}
		}
		GslBlock collector;
		gsl_vector *photon = collector.VectorAlloc(10), *last = collector.VectorAlloc(10);
#ifdef RECORD_TRACE
		vector<double> qdq(12);
		IO::NumPy rec("Trace " + to_string(ray_number), 12);
#endif
		Integrator &&integrator = metric_->GetIntegrator(2);
		while (gsl_hypot(alpha1 - alpha0, beta1 - beta0) > epsilon * (1. + gsl_hypot(alpha1, beta1))) {
			alpha0 = alpha1;
			beta0 = beta1;
			InitializePhoton(photon->data, alpha0, beta0);
			metric_->LagrangianToHamiltonian(photon->data);
			int status = 0, fixed = 0;
			h = -1.;
			while (status <= 0) {
				if (status == -1) {
					h *= 0.5;
					fixed = 0;
				}
				gsl_vector_memcpy(last, photon);
				if (fixed)
					status = integrator.ApplyFixedStep(gsl_vector_ptr(photon, 9), h, photon->data);
				else
					status = integrator.Apply(gsl_vector_ptr(photon, 9), t_final_, &h, photon->data);
				cosph = (r_star * sin_theta_star * cos_phi_star - gsl_vector_get(photon, 1) * GSL_SIGN(gsl_vector_get(photon, 2)) * sin(gsl_vector_get(photon, 2)) * cos(gsl_vector_get(photon, 3))) * sin_theta_observer_ + (r_star * cos_theta_star - gsl_vector_get(photon, 1) * GSL_SIGN(gsl_vector_get(photon, 2)) * cos(gsl_vector_get(photon, 2))) * cos_theta_observer_;
				if (cosph > r_star * 1e-12) {
					gsl_vector_memcpy(photon, last);
					h *= 0.3;
					fixed = 1;
					integrator.Reset();
				} else {
					if (gsl_vector_get(photon, 9) < t_final_ * 1e-8) {
						gsl_vector_set(photon, 8, gsl_vector_get(photon, 8) + gsl_vector_get(photon, 9));
						gsl_vector_set(photon, 9, 0);
					}
#ifdef RECORD_TRACE
					copy(photon, photon + 8, qdq.begin());
					metric_->HamiltonianToLagrangian(qdq.data());
					qdq[8] = photon[8] + photon[9];
					qdq[9] = metric::energy(qdq.data());
					qdq[10] = metric::angularMomentum(qdq.data());
					qdq[11] = metric::carter(qdq.data(), 0.);
					rec.save(qdq);
#endif
					if (cosph >= 0) {
						alpha1 += iteration_coefficient * (r_star * sin_theta_star * sin_phi_star - gsl_vector_get(photon, 1) * GSL_SIGN(gsl_vector_get(photon, 2)) * sin(gsl_vector_get(photon, 2)) * sin(gsl_vector_get(photon, 3))); // TODO:OPT!
						beta1 += iteration_coefficient * ((r_star * cos_theta_star - gsl_vector_get(photon, 1) * GSL_SIGN(gsl_vector_get(photon, 2)) * cos(gsl_vector_get(photon, 2))) * sin_theta_observer_ - (r_star * sin_theta_star * cos_phi_star - gsl_vector_get(photon, 1) * GSL_SIGN(gsl_vector_get(photon, 2)) * sin(gsl_vector_get(photon, 2)) * cos(gsl_vector_get(photon, 3))) * cos_theta_observer_);
						metric_->HamiltonianToLagrangian(photon->data);
						break;
					}
				}
			}
			if (status > 0) {
				fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
				return status;
			}
		}
		double photon2[9], cone_record[sample_number][3], local_cone_record[sample_number][3], area_record_initial[sample_number][3], area_record[sample_number][3], center_photon_position[3], center_photon_velocity[3], interval = M_2PI / sample_number;
		auto photon2_position_view = gsl_vector_view_array(photon2, 4), photon2_velocity_view = gsl_vector_view_array(photon2 + 4, 4);
		gsl_vector_set(photon, 0, 0.);
		copy(photon->data, photon->data + 8, photon2);
		photon2[8] = 0.;
		metric_->NormalizeNullGeodesic(photon2, 1.);
		metric_->LagrangianToHamiltonian(photon2);
		h = 1.;
		int status = 0, signum;
		integrator.Reset();
		while (status <= 0 && photon2[8] < 1000.)
			status = integrator.Apply(photon2 + 8, 1000., &h, photon2);
		if (status > 0) {
			fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
			return status;
		}
		metric_->HamiltonianToLagrangian(photon2);
		SphericalToCartesian(photon2);
		copy(photon2 + 1, photon2 + 4, center_photon_position);
		copy(photon2 + 5, photon2 + 8, center_photon_velocity);
		const double photon_norm = Norm(center_photon_velocity);
		for (int i = 0; i < 3; ++i)
			center_photon_velocity[i] /= photon_norm;
		auto gmunu = collector.MatrixAlloc(4, 4), coordinate = collector.MatrixAlloc(4, 4), coordinate_gmunu = collector.MatrixAlloc(4, 4), coordinate_static = collector.MatrixCalloc(4, 4), coordinate_static_gmunu = collector.MatrixAlloc(4, 4);
		auto permutation = collector.PermutationAlloc(4);
		double deltaA[4], Gamma[4][4][4];
		for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 4; ++j)
				for (int k = 0; k < 4; ++k)
					Gamma[i][j][k] = 0.;
		Gamma[0][0][1] = Gamma[0][1][0] = 1. / (gsl_vector_get(photon, 1) * (gsl_vector_get(photon, 1) - 2.));
		Gamma[1][0][0] = (gsl_vector_get(photon, 1) - 2.) / gsl_pow_3(gsl_vector_get(photon, 1));
		Gamma[1][1][1] = Gamma[0][0][1];
		Gamma[1][2][2] = 2. - gsl_vector_get(photon, 1);
		Gamma[1][3][3] = Gamma[1][2][2] * gsl_pow_2(sin(gsl_vector_get(photon, 2)));
		Gamma[2][1][2] = Gamma[2][2][1] = 1. / gsl_vector_get(photon, 1);
		Gamma[2][3][3] = -sin(gsl_vector_get(photon, 2)) * cos(gsl_vector_get(photon, 2));
		Gamma[3][1][3] = Gamma[3][3][1] = Gamma[2][1][2];
		Gamma[3][2][3] = Gamma[3][3][2] = cos(gsl_vector_get(photon, 2)) / sin(gsl_vector_get(photon, 2));
		metric_->MetricTensor(photon->data, gmunu);
		metric_->LocalInertialFrame(photon->data, coordinate);
		gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu, coordinate, 0., coordinate_gmunu);
		gsl_matrix_set(coordinate_static, 0, 0, -sqrt(-1. / gmunu->data[0]));
		gsl_matrix_set(coordinate_static, 1, 1, sqrt(1. / gmunu->data[5]));
		gsl_matrix_set(coordinate_static, 2, 2, sqrt(1. / gmunu->data[10]));
		gsl_matrix_set(coordinate_static, 3, 3, sqrt(1. / gmunu->data[15]));
		gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu, coordinate_static, 0., coordinate_static_gmunu);
		auto photon_transform = collector.VectorAlloc(4), photon_in_static_frame_cartesian = collector.VectorAlloc(4), photon_in_static_frame_spherical = collector.VectorAlloc(4), photon_in_star_frame_cartesian = collector.VectorAlloc(4), photon_in_star_frame_spherical = collector.VectorAlloc(4);
		gsl_vector_set(photon_transform, 0, 1.);
		copy(gsl_vector_ptr(photon, 5), gsl_vector_ptr(photon, 8), gsl_vector_ptr(photon_transform, 1));
		gsl_blas_dgemv(CblasNoTrans, 1., coordinate_gmunu, photon_transform, 0., photon_in_star_frame_cartesian);
		gsl_vector_scale(photon_in_star_frame_cartesian, 1. / gsl_vector_get(photon_in_star_frame_cartesian, 0));
		gsl_blas_dgemv(CblasNoTrans, 1., coordinate_static_gmunu, photon_transform, 0., photon_in_static_frame_cartesian);
		gsl_vector_scale(photon_in_static_frame_cartesian, 1. / gsl_vector_get(photon_in_static_frame_cartesian, 0));
		CartesianToSpherical(gsl_vector_ptr(photon_in_static_frame_cartesian, 1), gsl_vector_ptr(photon_in_static_frame_spherical, 1));
		CartesianToSpherical(gsl_vector_ptr(photon_in_star_frame_cartesian, 1), gsl_vector_ptr(photon_in_star_frame_spherical, 1));
		gsl_linalg_LU_decomp(coordinate_gmunu, permutation, &signum);
		for (int i = 0; i < sample_number; ++i) {
			const double angle = i * interval, sin_angle = sin(angle), cos_angle = cos(angle);
			copy(photon->data, photon->data + 4, photon2);
			photon2[4] = 1.;
			photon2[5] = cos_angle * sin_epsilon;
			photon2[6] = sin_angle * sin_epsilon;
			photon2[7] = cos_epsilon;
			RotateAroundAxis(photon2 + 5, 1, gsl_vector_get(photon_in_star_frame_spherical, 2));
			RotateAroundAxis(photon2 + 5, 2, gsl_vector_get(photon_in_star_frame_spherical, 3));
			gsl_linalg_LU_svx(coordinate_gmunu, permutation, &photon2_velocity_view.vector);
			photon2[8] = 0.;
			metric_->NormalizeNullGeodesic(photon2);
			local_cone_record[i][0] = metric_->DotProduct(photon->data, photon2 + 4, coordinate_static->data + 4, 4);
			local_cone_record[i][1] = metric_->DotProduct(photon->data, photon2 + 4, coordinate_static->data + 8, 4);
			local_cone_record[i][2] = metric_->DotProduct(photon->data, photon2 + 4, coordinate_static->data + 12, 4);
			const double vph_norm = metric_->DotProduct(photon->data, photon2 + 4, coordinate_static->data, 4);
			for (int j = 0; j < 3; ++j)
				local_cone_record[i][j] /= vph_norm;
			metric_->LagrangianToHamiltonian(photon2);
			h = 1.;
			int status = 0;
			integrator.Reset();
			while (status <= 0 && photon2[8] < 1000.)
				status = integrator.Apply(photon2 + 8, 1000., &h, photon2);
			if (status > 0) {
				fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
				return status;
			}
			metric_->HamiltonianToLagrangian(photon2);
			SphericalToCartesian(photon2);
			copy(photon2 + 5, photon2 + 8, cone_record[i]);
			const double rec2_norm = Norm(cone_record[i]);
			for (int j = 0; j < 3; ++j)
				cone_record[i][j] /= rec2_norm;
		}
		/*
		for (int i = 0; i < sample_number; ++i) { // FIXME: https://en.wikipedia.org/wiki/Raychaudhuri_equation
			const double angle = i * interval, sin_angle = sin(angle), cos_angle = cos(angle);
			photon2[0] = 0.;
			photon2[1] = epsilon * cos_angle;
			photon2[2] = epsilon * sin_angle;
			photon2[3] = 0.;
			photon2[4] = 1.;
			RotateAroundAxis(photon2 + 1, 1, gsl_vector_get(photon_in_star_frame_spherical, 2));
			RotateAroundAxis(photon2 + 1, 2, gsl_vector_get(photon_in_star_frame_spherical, 3));
			gsl_linalg_LU_svx(coordinate_gmunu, permutation, &photon2_position_view.vector);
			area_record_initial[i][0] = metric_->DotProduct(photon->data, photon2, coordinate_static->data + 4, 4);
			area_record_initial[i][1] = metric_->DotProduct(photon->data, photon2, coordinate_static->data + 8, 4);
			area_record_initial[i][2] = metric_->DotProduct(photon->data, photon2, coordinate_static->data + 12, 4);
			photon2[5] = gsl_vector_get(photon_in_static_frame_cartesian, 1);
			photon2[6] = gsl_vector_get(photon_in_static_frame_cartesian, 2);
			photon2[7] = gsl_vector_get(photon_in_static_frame_cartesian, 3);
			for (int j = 0; j < 4; ++j) {
				deltaA[j] = gsl_vector_get(photon_transform, j);
				for (int k = 0; k < 4; ++k)
					for (int l = 0; l < 4; ++l)
						deltaA[j] += Gamma[k][j][l] * photon2[4 + k] * photon2[l];
			}
			RotateAroundAxis(photon2 + 5, 2, gsl_vector_get(photon, 2));
			RotateAroundAxis(photon2 + 5, 0, -photon2[3]);
			RotateAroundAxis(photon2 + 5, 2, -gsl_vector_get(photon, 2) - photon2[2]);
			photon2[1] += gsl_vector_get(photon, 1);
			photon2[2] += gsl_vector_get(photon, 2);
			photon2[3] += gsl_vector_get(photon, 3);
			photon2[5] /= sqrt(photon2[1] / (photon2[1] - 2.));
			photon2[6] /= photon2[1];
			photon2[7] /= (photon2[1] * sin(photon2[2]));
			photon2[5] = deltaA[1];
			photon2[6] = deltaA[2];
			photon2[7] = deltaA[3];
			// SphericalToCartesian(area_record_initial[i], 3);
			photon2[8] = 0.;
			metric_->NormalizeNullGeodesic(photon2);
			metric_->LagrangianToHamiltonian(photon2);
			h = 1.;
			int status = 0;
			integrator.Reset();
			while (status <= 0 && photon2[8] < 1000.)
				status = integrator.Apply(photon2 + 8, 1000., &h, photon2);
			if (status > 0) {
				fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
				return status;
			}
			metric_->HamiltonianToLagrangian(photon2);
			copy(photon2 + 1, photon2 + 4, area_record[i]);
			SphericalToCartesian(area_record[i], 3);
		}*/
		double cone_solid_angle = DotCross(center_photon_velocity, cone_record[0], cone_record[sample_number - 1]);
		double cone_local_solid_angle = DotCross(photon_in_static_frame_cartesian->data + 1, local_cone_record[0], local_cone_record[sample_number - 1]);
		// SphericalToCartesian(gsl_vector_ptr(photon, 1), 3);
		double cross_section_area_initial = DotCross(gsl_vector_ptr(photon_in_static_frame_cartesian, 1), area_record_initial[0], area_record_initial[sample_number - 1]);
		double cross_section_area = TriangleArea(center_photon_position, area_record[0], area_record[sample_number - 1]);
		for (int i = 1; i < sample_number; ++i) {
			cone_solid_angle += DotCross(center_photon_velocity, cone_record[i], cone_record[i - 1]);
			cone_local_solid_angle += DotCross(photon_in_static_frame_cartesian->data + 1, local_cone_record[i], local_cone_record[i - 1]);
			cross_section_area_initial += DotCross(gsl_vector_ptr(photon_in_static_frame_cartesian, 1), area_record_initial[i], area_record_initial[i - 1]);
			cross_section_area += TriangleArea(center_photon_position, area_record[i], area_record[i - 1]);
		}
		if (ray_number == 2500) {
			NumPy area_record_npy("area_record", {3}), cone_record_npy("cone_record", {3});
			area_record_npy.Save(center_photon_position, 3);
			cone_record_npy.Save(center_photon_velocity, 3);
			for (int i = 0; i < sample_number; ++i) {
				area_record_npy.Save(area_record[i], 3);
				cone_record_npy.Save(cone_record[i], 3);
			}
		}
#ifdef VIEW_TAU
		output->Save({alpha1, beta1, star.RedshiftTau(photon), (photon[8] + photon[9]) / Unit::s, 0.5 * abs(cone_solid_angle) / epsilon_circle_area, sqrt(-1. / gmunu->data[0]), 0.5 * abs(cone_local_solid_angle) / epsilon_circle_area, cross_section_area / epsilon_circle_area, cross_section_area_initial / epsilon_circle_area});
#else
		output_->Save({alpha1, beta1, star.Redshift(photon->data), (gsl_vector_get(photon, 8) + gsl_vector_get(photon, 9)) / Unit::s, 0.5 * abs(cone_solid_angle) / epsilon_circle_area, sqrt(-1. / gmunu->data[0]), 0.5 * abs(cone_local_solid_angle) / epsilon_circle_area, cross_section_area / epsilon_circle_area, 0.5 * abs(cross_section_area_initial) / epsilon_circle_area});
#endif
		return 0;
	}
	int View::OmegaTest() {
		Integrator &&integrator = metric_->GetIntegrator(2);
		double position_[8] = {0., 3., M_PI_4, 0., 0., 0., 0., 0.};
		metric_->NormalizeTimelikeGeodesic(position_);
		gsl_matrix *coordinate = gsl_matrix_alloc(4, 4), *gmunu = gsl_matrix_alloc(4, 4), *coordinate_gmunu = gsl_matrix_alloc(4, 4);
		gsl_permutation *perm = gsl_permutation_alloc(4);
		metric_->MetricTensor(position_, gmunu);
		gsl_matrix_set_zero(coordinate);
		gsl_matrix_set(coordinate, 0, 0, sqrt(-1. / gmunu->data[0]));
		gsl_matrix_set(coordinate, 1, 1, sqrt(1. / gmunu->data[5]));
		gsl_matrix_set(coordinate, 2, 2, sqrt(1. / gmunu->data[10]));
		gsl_matrix_set(coordinate, 3, 3, sqrt(1. / gmunu->data[15]));
		gsl_matrix_set_zero(coordinate_gmunu);
		gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu, coordinate, 0., coordinate_gmunu);
		int signum;
		gsl_linalg_LU_decomp(coordinate_gmunu, perm, &signum);
		double rec[100][3], center_ph[3];
		double ph[9], area, h, time_limit = 0., interval = M_2PI / sample_number;
		gsl_vector_view ph_view = gsl_vector_view_array(ph + 4, 4);
		NumPy cone_record("Omega_record", {1});
		ProgressBar::bars_[0].set_option(indicators::option::MaxProgress(90));
		for (double angle = 90; angle > 0; angle -= 1) {
			const double sina = sin(angle / 180. * M_PI), cosa = cos(angle / 180. * M_PI);
			copy(position_, position_ + 4, ph);
			ph[4] = 1.;
			ph[5] = sina;
			ph[6] = cosa;
			ph[7] = 0.;
			gsl_linalg_LU_svx(coordinate_gmunu, perm, &ph_view.vector);
			metric_->NormalizeNullGeodesic(ph);
			if (ph[5] < 0) {
				ph[5] = -ph[5];
				ph[6] = -ph[6];
				ph[7] = -ph[7];
			}
			metric_->LagrangianToHamiltonian(ph);
			ph[8] = 0;
			h = 1.;
			int status = 0;
			integrator.Reset();
			while (status <= 0 && ph[1] < 1.e3)
				status = integrator.Apply(&time_limit, -t_final_, &h, ph);
			metric_->HamiltonianToLagrangian(ph);
			const double sin_theta = GSL_SIGN(ph[2]) * sin(ph[2]), cos_theta = GSL_SIGN(ph[2]) * cos(ph[2]), sin_phi = sin(ph[3]), cos_phi = cos(ph[3]);
			center_ph[0] = ph[5] * sin_theta * cos_phi + ph[1] * (cos_theta * cos_phi * ph[6] - sin_theta * sin_phi * ph[7]);
			center_ph[1] = ph[5] * sin_theta * sin_phi + ph[1] * (cos_theta * sin_phi * ph[6] + sin_theta * cos_phi * ph[7]);
			center_ph[2] = ph[5] * cos_theta - ph[1] * sin_theta * ph[6];
			const double vph_norm = Norm(center_ph);
			for (int j = 0; j < 3; ++j)
				center_ph[j] /= vph_norm;
			for (int i = 0; i < 100; ++i) {
				const double angle_i = i * interval, sinai = sin(angle_i), cosai = cos(angle_i);
				copy(position_, position_ + 4, ph);
				ph[4] = 1.;
				ph[5] = sina - sin_epsilon * cosai * cosa;
				ph[6] = cosa + sin_epsilon * cosai * sina;
				ph[7] = sin_epsilon * sinai;
				gsl_linalg_LU_svx(coordinate_gmunu, perm, &ph_view.vector);
				metric_->NormalizeNullGeodesic(ph);
				if (ph[5] < 0) {
					ph[5] = -ph[5];
					ph[6] = -ph[6];
					ph[7] = -ph[7];
				}
				metric_->LagrangianToHamiltonian(ph);
				ph[8] = 0;
				h = 1.;
				status = 0;
				integrator.Reset();
				while (status <= 0 && ph[8] < time_limit)
					status = integrator.Apply(ph + 8, time_limit, &h, ph);
				metric_->HamiltonianToLagrangian(ph);
				const double sin_theta = GSL_SIGN(ph[2]) * sin(ph[2]), cos_theta = GSL_SIGN(ph[2]) * cos(ph[2]), sin_phi = sin(ph[3]), cos_phi = cos(ph[3]);
				rec[i][0] = ph[5] * sin_theta * cos_phi + ph[1] * (cos_theta * cos_phi * ph[6] - sin_theta * sin_phi * ph[7]);
				rec[i][1] = ph[5] * sin_theta * sin_phi + ph[1] * (cos_theta * sin_phi * ph[6] + sin_theta * cos_phi * ph[7]);
				rec[i][2] = ph[5] * cos_theta - ph[1] * sin_theta * ph[6];
				const double vph_norm = Norm(rec[i]);
				for (int j = 0; j < 3; ++j)
					rec[i][j] /= vph_norm;
			}
			area = DotCross(center_ph, rec[0], rec[99]);
			for (int i = 1; i < sample_number; ++i)
				area += DotCross(center_ph, rec[i], rec[i - 1]);
			if (angle == 13) {
				NumPy cone_record("cone_record13", {3});
				cone_record.Save(center_ph, 3);
				for (int i = 0; i < sample_number; ++i)
					cone_record.Save(rec[i], 3);
			}
			cone_record.Save({abs(area) / (M_2PI * gsl_pow_2(epsilon))});
			ProgressBar::bars_[0].tick();
		}
		return 0;
	}
	int View::Shadow() {
		const double interval = M_2PI / sample_number;
		NumPy rec("shadow", {2});
		double h, rin = 2., rout = 10., rmid = 6., photon[10];
		Integrator &&integrator = metric_->GetIntegrator(2);
		bar_->set_option(indicators::option::MaxProgress(sample_number));
		bar_->set_option(indicators::option::PrefixText("? Shadow"));
		int progressBarIndex = ProgressBar::bars_.push_back(*bar_);
		for (int i = 0; i < sample_number; ++i) {
			const double angle = i * interval, sina = sin(angle), cosa = cos(angle);
			int status = 0;
			while (rout - rin > epsilon * (rin + rout)) {
				rmid = 0.5 * (rin + rout);
				InitializePhoton(photon, rmid * cosa, rmid * sina);
				h = -1.;
				while (status <= 0 && photon[8] + photon[9] > t_final_) {
					status = integrator.Apply(photon + 9, t_final_, &h, photon);
					if (photon[9] < t_final_ * 1e-8) {
						photon[8] += photon[9];
						photon[9] = 0.;
					}
					if (photon[4] >= 1e6 || photon[5] <= 0 || photon[1] <= 0)
						break;
				}
				if (status > 0) {
					fmt::print(stderr, "[!] view::shadow status = {}\n", status);
					return status;
				}
				if (photon[5] <= 0)
					rout = rmid;
				else
					rin = rmid;
				integrator.Reset();
			}
			rin -= 2. * interval * rmid;
			rout += 2. * interval * rmid;
			rec.Save({rmid * cosa, rmid * sina});
			if (ProgressBar::display_)
				ProgressBar::bars_[progressBarIndex].tick();
		}
		if (ProgressBar::display_)
			ProgressBar::SetComplete(progressBarIndex, "! Shadow");
		return 0;
	}
	Camera::Camera(shared_ptr<Metric> metric, size_t pixel, double half_angle, double r, double theta, string file_name) : View(move(metric), r, theta, file_name), pixel_(pixel), half_angle_(half_angle) {
		screen_ = vector<vector<double>>(pixel, vector<double>(pixel));
		initials_ = vector<array<double, 10>>(pixel * pixel);
		const double pixel_size = 2. * half_angle * r / pixel, t1 = r + 100.;
		for (size_t i = 0; i < pixel; ++i)
			for (size_t j = 0; j < pixel; ++j)
				InitializePhoton(initials_[i * pixel + j].data(), pixel_size * (i - 0.5 * pixel + 0.5), pixel_size * (j - 0.5 * pixel + 0.5));
#pragma omp parallel for
		for (int p = pixel * pixel - 1; p >= 0; --p) {
			int status = 0;
			initials_[p][9] = -1.;
			Integrator &&integrator = metric_->GetIntegrator(2);
			metric_->NormalizeNullGeodesic(initials_[p].data(), 1.);
			metric_->LagrangianToHamiltonian(initials_[p].data());
			while (status <= 0 && initials_[p][8] > t1 && initials_[p][1] > 100)
				status = integrator.Apply(initials_[p].data() + 8, t1, initials_[p].data() + 9, initials_[p].data());
			if (status > 0)
				fmt::print(stderr, "[!] camera::initialize status = {}\n", status);
		}
	}
	int Camera::TraceStar() {
		const double t1 = -1000.;
#pragma omp parallel for
		for (int p = pixel_ * pixel_ - 1; p >= 0; --p) {
			int i = p / pixel_;
			int j = p - i * pixel_;
			double ph[10], last[10];
			int status = 0;
			copy(initials_[p].begin(), initials_[p].end(), ph);
			Integrator &&integrator = metric_->GetIntegrator(2);
			while (status <= 0 && ph[8] > t1) {
				copy(ph, ph + 10, last);
				status = integrator.Apply(ph + 8, t1, ph + 9, ph);
				for (auto objP : Object::object_list_)
					if (objP->Hit(ph, last))
						screen_[i][j] = objP->Redshift(ph); // FIXME: if multi objects
				if (screen_[i][j] > epsilon)
					break;
			}
			if (status > 0)
				fmt::print(stderr, "[!] camera::traceStar status = {}\n", status);
		}
		return 0;
	}
	int Camera::Lens() {
		const double t1 = -1000. * r_, pixelPerAngle = 0.5 * pixel_ / half_angle_;
		NumPy rec("lens", {2});
		bar_->set_option(indicators::option::MaxProgress(pixel_ * pixel_));
		bar_->set_option(indicators::option::PrefixText("? Lens"));
		int progressBarIndex = ProgressBar::bars_.push_back(*bar_);
		for (size_t i = 0; i < pixel_; ++i)
			for (size_t j = 0; j < pixel_; ++j) {
				double ph[10];
				int status = 0;
				copy(initials_[i * pixel_ + j].begin(), initials_[i * pixel_ + j].end(), ph);
				Integrator &&integrator = metric_->GetIntegrator(2);
				while (status <= 0 && ph[8] > t1 && ph[1] > 3. && ph[1] < 3.e2)
					status = integrator.Apply(ph + 8, t1, ph + 9, ph);
				if (status > 0)
					fmt::print(stderr, "[!] camera::lens status = {}\n", status);
				if (ph[1] <= 3.)
					rec.Save({NAN, NAN});
				else {
					metric_->HamiltonianToLagrangian(ph);
					if (ph[2] < 0)
						ph[2] += M_PI;
					SphericalToCartesian(ph);
					rec.Save({ph[6] * pixelPerAngle, (ph[7] * sin_theta_observer_ - ph[5] * cos_theta_observer_) * pixelPerAngle});
				}
				if (ProgressBar::display_)
					ProgressBar::bars_[progressBarIndex].tick();
			}
		if (ProgressBar::display_)
			ProgressBar::SetComplete(progressBarIndex, "! lens");
		return 0;
	}
	int Camera::Save() {
		for (const vector<double> &line : screen_)
			output_->Save(line);
		return 0;
	}
} // namespace SBody
