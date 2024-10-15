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
	View::View(shared_ptr<Metric> metric, double r, double theta, double iota) : metric_(metric), r_(r), r2_(gsl_pow_2(r)), theta_(theta), sin_theta_(sin(theta)), cos_theta_(cos(theta)), iota_(iota), sin_iota_(sin(iota)), cos_iota_(cos(iota)), t_final_(-2e4) {
		bar_ = make_unique<indicators::BlockProgressBar>(indicators::option::ShowElapsedTime{true},
														 indicators::option::ShowRemainingTime{true},
														 indicators::option::ForegroundColor{indicators::Color(4)},
														 indicators::option::FontStyles{vector<indicators::FontStyle>{indicators::FontStyle::bold}});
	}
	int View::InitializePhoton(double photon[], double alpha, double beta) {
		return metric_->InitializePhoton(photon, alpha, beta, r_, r2_, theta_, sin_theta_);
	}
	int View::Trace(const double position[], TimeSystem object_time, double record[], bool calculate_magnification, bool fast_trace) {
		double photon[9], last_step_record[9], alpha, delta_alpha, beta, delta_beta, h;
		if (fast_trace && metric_->FastTrace(r_, theta_, sin_theta_, cos_theta_, position[1], position[2], position[3], alpha, beta, photon) == GSL_SUCCESS) {
			PhotonInformation(position, object_time, record, photon, alpha, beta);
			return calculate_magnification ? Magnification(position, object_time, record[4], photon, record[2]) : GSL_SUCCESS;
		}
		auto t_start = chrono::steady_clock::now();
		const double r_object = position[1], sin_theta_object = abs(sin(position[2])), cos_theta_object = GSL_SIGN(position[2]) * cos(position[2]), sin_phi_object = sin(position[3]), cos_phi_object = cos(position[3]);
		if (r_object <= 3.) {
			PrintlnWarning("Object orbit radius = {:.6f}", r_object);
			if (r_object < 0)
				return GSL_FAILURE;
		}
		const double alpha_coefficient = sin_theta_object * sin_phi_object, beta_coefficient = cos_theta_object * sin_theta_ - sin_theta_object * cos_phi_object * cos_theta_, cos_observer_object = sin_theta_ * sin_theta_object * cos_phi_object + cos_theta_ * cos_theta_object, sin_observer_object = sqrt(gsl_pow_2(alpha_coefficient) + gsl_pow_2(beta_coefficient)), theta_observer_object = acos(cos_observer_object), iteration_coefficient = tanh((2. + cos_observer_object) * 0.02 * r_object);
		if (cos_observer_object == -1.) {
			PrintlnWarning("Object behind black hole, cos(theta) = {:.6f}\n", cos_observer_object);
			alpha = 2. * sqrt(r_object);
			beta = 0.;
		} else if (cos_observer_object < 1.) { // initial guessing
			if (theta_observer_object < M_PI_2) {
				const double effective_radius = r_object + gsl_pow_3(theta_observer_object / M_PI_2) / sin_observer_object;
				alpha = effective_radius * alpha_coefficient;
				beta = effective_radius * beta_coefficient;
			} else {
				const double effective_radius = 1. / sin_observer_object + 0.5 * (M_PI * r_object - 6. * cos_observer_object) / (M_PI - theta_observer_object - sin_observer_object * cos_observer_object);
				// b-r*sin(theta)=(b-1.)*(2.*theta/pi-1.)+1.+(b-4.)/pi*sin(theta*2.)
				// b=(r*sin(theta)*M_PI+M_2_PI-2.*theta-8.*sin(theta)*cos(theta))/(2.*(M_PI-theta-sin(theta)*cos(theta)))
				alpha = effective_radius * alpha_coefficient;
				beta = effective_radius * beta_coefficient;
			}
		}
		unique_ptr<Integrator> integrator = metric_->GetIntegrator(T, HAMILTONIAN);
		int retry = 0;
		while (true) {
			InitializePhoton(photon, alpha, beta);
			metric_->LagrangianToHamiltonian(photon);
			int status = 0, fixed = 0;
			h = -0.01 * r_;
			while (status <= 0) {
#ifdef GSL_RANGE_CHECK_OFF
				if (auto now = chrono::steady_clock::now(); now - t_start > chrono::seconds(5)) {
					PrintlnWarning("View::Trace() timeout! Object position is:\nr = {}\ntheta = {}\nphi = {}", position[1], position[2], position[3]);
					return GSL_FAILURE;
				}
#endif
				if (status == GSL_FAILURE) {
					h *= 0.5;
					fixed = 0;
				}
				copy(photon, photon + 9, last_step_record);
				if (fixed)
					status = integrator->ApplyFixedStep(photon + 8, h, photon);
				else
					status = integrator->ApplyStep(photon + 8, t_final_, &h, photon);
				if (const double cos_observer_object_photon = (r_object * sin_theta_object * cos_phi_object - photon[1] * abs(sin(photon[2])) * cos(photon[3])) * sin_theta_ + (r_object * cos_theta_object - photon[1] * GSL_SIGN(photon[2]) * cos(photon[2])) * cos_theta_; cos_observer_object_photon > r_object * 1e-13) {
					// photon goes through the plane of the object
					copy(last_step_record, last_step_record + 9, photon);
					h *= 0.3;
					fixed = 1;
					integrator->Reset();
				} else if (cos_observer_object_photon >= 0) {
					// photon in the same plane with the object
					delta_alpha = r_object * sin_theta_object * sin_phi_object - photon[1] * abs(sin(photon[2])) * sin(photon[3]);
					delta_beta = (r_object * cos_theta_object - photon[1] * GSL_SIGN(photon[2]) * cos(photon[2])) * sin_theta_ - (r_object * sin_theta_object * cos_phi_object - photon[1] * abs(sin(photon[2])) * cos(photon[3])) * cos_theta_;
					break;
				}
				// photon not reaches the plane of the object
				if (photon[8] > t_final_)
					continue;
				// photon runs away
				if (++retry < 3) { // retry limit
					delta_alpha = alpha * 0.2;
					delta_beta = beta * 0.2;
					break;
				}
				PrintlnWarning("View::Trace() failed! Object position is:\nr = {}\ntheta = {}\nphi = {}", position[1], position[2], position[3]);
				return GSL_FAILURE;
			}
			if (status > 0) {
				PrintlnWarning("View::Trece() status = {}\n", status);
				return status;
			}
			if (gsl_hypot(delta_alpha, delta_beta) <= EPSILON * (1. + gsl_hypot(alpha, beta))) {
				metric_->HamiltonianToLagrangian(photon);
				break;
			}
			alpha += iteration_coefficient * delta_alpha;
			beta += iteration_coefficient * delta_beta;
			integrator->Reset();
		}
		photon[8] -= r_;
		PhotonInformation(position, object_time, record, photon, alpha, beta);
		return calculate_magnification ? Magnification(position, object_time, record[4], photon, record[2]) : GSL_SUCCESS;
	}
	int View::PhotonInformation(const double position[], TimeSystem object_time, double record[], const double photon[], double alpha, double beta) {
		record[0] = alpha * cos_iota_ + beta * sin_iota_;				 // alpha
		record[1] = beta * cos_iota_ - alpha * sin_iota_;				 // beta
		record[2] = metric_->Redshift(position, photon, object_time, T); // redshift
		record[3] = photon[8] / Unit::s;								 // look back time
		return GSL_SUCCESS;
	}
	int View::Magnification(const double position[], TimeSystem object_time, double &magnification, const double photon[], double redshift) {
		unique_ptr<Integrator> integrator = metric_->GetIntegrator(T, HAMILTONIAN);
		double forward_photon[8], cone_record[SAMPLE_NUMBER][3], local_cone_record[SAMPLE_NUMBER][3], center_photon_velocity[3], h = 1., t = 0.;
		auto forward_photon_velocity_view = gsl_vector_view_array(forward_photon + 4, 4);
		copy(photon, photon + 8, forward_photon);
		metric_->NormalizeNullGeodesic(forward_photon);
		metric_->LagrangianToHamiltonian(forward_photon);
		integrator->Reset();
		if (int status = integrator->Apply(&t, 1000., &h, forward_photon); status > 0) {
			PrintlnError("View::Magnification() status = {}", status);
			return status;
		}
		metric_->HamiltonianToLagrangian(forward_photon);
		SphericalToCartesian(forward_photon);
		copy(forward_photon + 5, forward_photon + 8, center_photon_velocity);
		cblas_dscal(3, 1. / Norm(center_photon_velocity), center_photon_velocity, 1);
		GslBlock collector;
		auto gmunu = collector.MatrixAlloc(4, 4);			 // object local metric tensor
		auto coordinate = collector.MatrixAlloc(4, 4);		 // object local inertial coordinate frame
		auto coordinate_gmunu = collector.MatrixAlloc(4, 4); // object local inertial frame measured by observer
		auto permutation = collector.PermutationAlloc(4);	 // permutation used in the LU decomposition
		auto photon_transform = collector.VectorAlloc(4);	 // photon in TimeSystem TAU
		auto photon_in_object_frame_cartesian = collector.VectorAlloc(4);
		double photon_in_object_frame_spherical[4];
		metric_->MetricTensor(position, gmunu);
		metric_->LocalInertialFrame(position, object_time, coordinate);
		gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu, coordinate, 0., coordinate_gmunu);
		gsl_vector_set(photon_transform, 0, 1.);
		copy(photon + 5, photon + 8, gsl_vector_ptr(photon_transform, 1));
		gsl_blas_dgemv(CblasNoTrans, 1., coordinate_gmunu, photon_transform, 0., photon_in_object_frame_cartesian);
		// gsl_vector_scale(photon_in_object_frame_cartesian, 1. / gsl_vector_get(photon_in_object_frame_cartesian, 0));
		CartesianToSpherical(photon_in_object_frame_cartesian->data, photon_in_object_frame_spherical, 4);
#ifndef GSL_RANGE_CHECK_OFF
		auto coordinate_static = collector.MatrixCalloc(4, 4);		// object local static frame (only dt/d\tau != 0)
		auto coordinate_static_gmunu = collector.MatrixAlloc(4, 4); // object local static frame measured by observer
		auto photon_in_static_frame_cartesian = collector.VectorAlloc(4);
		gsl_matrix_set(coordinate_static, 0, 0, sqrt(-1. / gmunu->data[0]));
		gsl_matrix_set(coordinate_static, 1, 1, sqrt(1. / gmunu->data[5]));
		gsl_matrix_set(coordinate_static, 2, 2, sqrt(1. / gmunu->data[10]));
		gsl_matrix_set(coordinate_static, 3, 3, sqrt(1. / gmunu->data[15]));
		gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu, coordinate_static, 0., coordinate_static_gmunu);
		gsl_blas_dgemv(CblasNoTrans, 1., coordinate_static_gmunu, photon_transform, 0., photon_in_static_frame_cartesian);
		// local_redshift should equal to sqrt(EPSILON_POLYGON_AREA / cone_local_solid_angle),
		// the main error comes from the calculation of the cone_local_solid_angle, due to the limited SAMPLE_NUMBER.
		const double local_redshift = photon_in_object_frame_cartesian->data[0] / photon_in_static_frame_cartesian->data[0];
		gsl_vector_scale(photon_in_static_frame_cartesian, 1. / gsl_vector_get(photon_in_static_frame_cartesian, 0));
#endif
		int signum;
		gsl_linalg_LU_decomp(coordinate_gmunu, permutation, &signum);
		for (int i = 0; i < SAMPLE_NUMBER; ++i) {
			const double angle = i * ANGLE_INTERVAL;
			copy(photon, photon + 4, forward_photon);
			forward_photon[4] = -1.;
			forward_photon[5] = cos(angle) * SIN_EPSILON;
			forward_photon[6] = sin(angle) * SIN_EPSILON;
			forward_photon[7] = COS_EPSILON;
			RotateAroundAxis(forward_photon + 5, Y, photon_in_object_frame_spherical[2]);
			RotateAroundAxis(forward_photon + 5, Z, photon_in_object_frame_spherical[3]);
			gsl_linalg_LU_svx(coordinate_gmunu, permutation, &forward_photon_velocity_view.vector);
			metric_->NormalizeNullGeodesic(forward_photon);
#ifndef GSL_RANGE_CHECK_OFF
			// local_cone_record[i][0] = metric_->DotProduct(photon, forward_photon + 4, coordinate_static->data + 4, 4);
			local_cone_record[i][0] = cblas_ddot(4, forward_photon + 4, 1, coordinate_static_gmunu->data + 4, 1);
			// local_cone_record[i][1] = metric_->DotProduct(photon, forward_photon + 4, coordinate_static->data + 8, 4);
			local_cone_record[i][1] = cblas_ddot(4, forward_photon + 4, 1, coordinate_static_gmunu->data + 8, 1);
			// ocal_cone_record[i][2] = metric_->DotProduct(photon, forward_photon + 4, coordinate_static->data + 12, 4);
			local_cone_record[i][2] = cblas_ddot(4, forward_photon + 4, 1, coordinate_static_gmunu->data + 12, 1);
			cblas_dscal(3, 1. / cblas_ddot(4, forward_photon + 4, 1, coordinate_static_gmunu->data, 1), local_cone_record[i], 1);
#endif
			metric_->LagrangianToHamiltonian(forward_photon);
			h = 1.;
			t = 0.;
			integrator->Reset();
			if (int status = integrator->Apply(&t, 1000., &h, forward_photon); status > 0) {
				PrintlnError("View::Magnification() status = {}", status);
				return status;
			}
			metric_->HamiltonianToLagrangian(forward_photon);
			SphericalToCartesian(forward_photon);
			copy(forward_photon + 5, forward_photon + 8, cone_record[i]);
			cblas_dscal(3, 1. / Norm(cone_record[i]), cone_record[i], 1);
		}
		double cone_solid_angle = TriangleArea(center_photon_velocity, cone_record[0], cone_record[SAMPLE_NUMBER - 1]);
#ifndef GSL_RANGE_CHECK_OFF
		double cone_local_solid_angle = TriangleArea(photon_in_static_frame_cartesian->data + 1, local_cone_record[0], local_cone_record[SAMPLE_NUMBER - 1]);
#endif
		for (int i = 1; i < SAMPLE_NUMBER; ++i) {
			cone_solid_angle += TriangleArea(center_photon_velocity, cone_record[0], cone_record[SAMPLE_NUMBER - 1]);
#ifndef GSL_RANGE_CHECK_OFF
			cone_local_solid_angle += TriangleArea(photon_in_static_frame_cartesian->data + 1, local_cone_record[i], local_cone_record[i - 1]);
#endif
		}
		magnification = EPSILON_POLYGON_AREA / (cone_solid_angle * redshift);
		return GSL_SUCCESS;
	}
	int View::OmegaTest() {
		unique_ptr<Integrator> integrator = metric_->GetIntegrator(T, HAMILTONIAN);
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
		double ph[9], area, h, time_limit = 0.;
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
			integrator->Reset();
			while (status <= 0 && ph[1] < 1.e3)
				status = integrator->Apply(&time_limit, -t_final_, &h, ph);
			metric_->HamiltonianToLagrangian(ph);
			const double sin_theta = GSL_SIGN(ph[2]) * sin(ph[2]), cos_theta = GSL_SIGN(ph[2]) * cos(ph[2]), sin_phi = sin(ph[3]), cos_phi = cos(ph[3]);
			center_ph[0] = ph[5] * sin_theta * cos_phi + ph[1] * (cos_theta * cos_phi * ph[6] - sin_theta * sin_phi * ph[7]);
			center_ph[1] = ph[5] * sin_theta * sin_phi + ph[1] * (cos_theta * sin_phi * ph[6] + sin_theta * cos_phi * ph[7]);
			center_ph[2] = ph[5] * cos_theta - ph[1] * sin_theta * ph[6];
			const double vph_norm = Norm(center_ph);
			for (int j = 0; j < 3; ++j)
				center_ph[j] /= vph_norm;
			for (int i = 0; i < 100; ++i) {
				const double angle_i = i * ANGLE_INTERVAL, sinai = sin(angle_i), cosai = cos(angle_i);
				copy(position_, position_ + 4, ph);
				ph[4] = 1.;
				ph[5] = sina - SIN_EPSILON * cosai * cosa;
				ph[6] = cosa + SIN_EPSILON * cosai * sina;
				ph[7] = SIN_EPSILON * sinai;
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
				integrator->Reset();
				while (status <= 0 && ph[8] < time_limit)
					status = integrator->Apply(ph + 8, time_limit, &h, ph);
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
			for (int i = 1; i < SAMPLE_NUMBER; ++i)
				area += DotCross(center_ph, rec[i], rec[i - 1]);
			if (angle == 13) {
				NumPy cone_record("cone_record13", {3});
				cone_record.Save(center_ph, 3);
				for (int i = 0; i < SAMPLE_NUMBER; ++i)
					cone_record.Save(rec[i], 3);
			}
			cone_record.Save({abs(area) / (M_2PI * gsl_pow_2(EPSILON))});
			ProgressBar::bars_[0].tick();
		}
		return 0;
	}
	int View::Shadow(string file_name) {
		NumPy record(file_name, {2});
		double h, rin = 2., rout = 10., rmid = 6., photon[10];
		unique_ptr<Integrator> integrator = metric_->GetIntegrator(T, HAMILTONIAN);
		bar_->set_option(indicators::option::MaxProgress(SAMPLE_NUMBER));
		bar_->set_option(indicators::option::PrefixText("? Shadow"));
		int progressBarIndex = ProgressBar::bars_.push_back(*bar_);
		for (int i = 0; i < SAMPLE_NUMBER; ++i) {
			const double angle = i * ANGLE_INTERVAL, sin_angle = sin(angle), cos_angle = cos(angle);
			int status = 0;
			while (rout - rin > EPSILON * (rin + rout)) {
				rmid = 0.5 * (rin + rout);
				InitializePhoton(photon, rmid * cos_angle, rmid * sin_angle);
				h = -1.;
				while (status <= 0 && photon[8] + photon[9] > t_final_) {
					status = integrator->Apply(photon + 9, t_final_, &h, photon);
					if (photon[9] < t_final_ * 1e-8) {
						photon[8] += photon[9];
						photon[9] = 0.;
					}
					if (photon[4] >= 1e6 || photon[5] <= 0 || photon[1] <= 0)
						break;
				}
				if (status > 0) {
					PrintlnError("View::Shadow() status = {}", status);
					return status;
				}
				if (photon[5] <= 0)
					rout = rmid;
				else
					rin = rmid;
				integrator->Reset();
			}
			rin -= 2. * ANGLE_INTERVAL * rmid;
			rout += 2. * ANGLE_INTERVAL * rmid;
			record.Save({rmid * (cos_angle * cos_iota_ - sin_angle * sin_iota_), rmid * (sin_angle * cos_iota_ + cos_angle * sin_iota_)});
			if (ProgressBar::display_)
				ProgressBar::bars_[progressBarIndex].tick();
		}
		if (ProgressBar::display_)
			ProgressBar::SetComplete(progressBarIndex, "! Shadow");
		return 0;
	}
	Camera::Camera(shared_ptr<Metric> metric, size_t pixel, double half_angle, double r, double theta, double iota) : View(metric, r, theta, iota), pixel_(pixel), half_angle_(half_angle) {
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
			unique_ptr<Integrator> integrator = metric_->GetIntegrator(T, HAMILTONIAN);
			metric_->NormalizeNullGeodesic(initials_[p].data(), 1.);
			metric_->LagrangianToHamiltonian(initials_[p].data());
			while (status <= 0 && initials_[p][8] > t1 && initials_[p][1] > 100)
				status = integrator->Apply(initials_[p].data() + 8, t1, initials_[p].data() + 9, initials_[p].data());
			if (status > 0)
				PrintlnWarning("Camera::Camera() status = {}", status);
		}
	}
	int Camera::Trace() {
		const double t1 = -1000.;
#pragma omp parallel for
		for (int p = pixel_ * pixel_ - 1; p >= 0; --p) {
			int i = p / pixel_;
			int j = p - i * pixel_;
			double ph[10], last[10];
			int status = 0;
			copy(initials_[p].begin(), initials_[p].end(), ph);
			unique_ptr<Integrator> integrator = metric_->GetIntegrator(T, HAMILTONIAN);
			while (status <= 0 && ph[8] > t1) {
				copy(ph, ph + 10, last);
				status = integrator->Apply(ph + 8, t1, ph + 9, ph);
				for (auto objP : Object::object_list_)
					if (objP->Hit(ph, last))
						screen_[i][j] = objP->Redshift(ph, T); // FIXME: if multi objects
				if (screen_[i][j] > EPSILON)
					break;
			}
			if (status > 0)
				PrintlnWarning("Camera::Trace() status = {}", status);
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
				unique_ptr<Integrator> integrator = metric_->GetIntegrator(T, HAMILTONIAN);
				while (status <= 0 && ph[8] > t1 && ph[1] > 3. && ph[1] < 3.e2)
					status = integrator->Apply(ph + 8, t1, ph + 9, ph);
				if (status > 0)
					PrintlnWarning("Camera::Lens() status = {}", status);
				if (ph[1] <= 3.)
					rec.Save({NAN, NAN});
				else {
					metric_->HamiltonianToLagrangian(ph);
					if (ph[2] < 0)
						ph[2] += M_PI;
					SphericalToCartesian(ph);
					rec.Save({ph[6] * pixelPerAngle, (ph[7] * sin_theta_ - ph[5] * cos_theta_) * pixelPerAngle});
				}
				if (ProgressBar::display_)
					ProgressBar::bars_[progressBarIndex].tick();
			}
		if (ProgressBar::display_)
			ProgressBar::SetComplete(progressBarIndex, "! lens");
		return 0;
	}
	int Camera::Save(string file_name) {
		NumPy record(file_name, {2});
		for (const vector<double> &line : screen_)
			record.Save(line);
		return 0;
	}
} // namespace SBody
