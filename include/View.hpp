/**
 * @file View.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_VIEW_H
#define SBODY_VIEW_H

#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <boost/algorithm/algorithm.hpp>
#include <fmt/core.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "IO.hpp"
#include "Metric.hpp"
#include "Object.hpp"
#include "Utility.hpp"

namespace SBody {
	template <typename Type>
	struct TraceParameters {
		std::shared_ptr<Metric<Type>> metric;
		std::shared_ptr<Integrator> integrator;
		double *photon;
		const Type r, r2;
		const Type theta_obs, sin_theta_obs, cos_theta_obs;
		const Type r_obj, x_obj, y_obj, z_obj;
		const Type t_final;
		TraceParameters(std::shared_ptr<Metric<Type>> metric, std::shared_ptr<Integrator> integrator, Type photon[], Type r, Type r2, Type theta_obs, Type sin_theta_obs, Type cos_theta_obs, Type r_obj, Type sin_theta_obj, Type cos_theta_obj, Type sin_phi_obj, Type cos_phi_obj, Type t_final) : metric(metric), integrator(integrator), photon(photon), r(r), r2(r2), theta_obs(theta_obs), sin_theta_obs(sin_theta_obs), cos_theta_obs(cos_theta_obs), r_obj(r_obj), x_obj(r_obj * sin_theta_obj * cos_phi_obj), y_obj(r_obj * sin_theta_obj * sin_phi_obj), z_obj(r_obj * cos_theta_obj), t_final(t_final) {}
	};
	template <typename Type>
	class View {
	  protected:
		std::shared_ptr<Metric<Type>> metric_;
		/// Distance to the center black hole
		const Type r_;
		/// Square of the distance to the center black hole
		const Type r2_;
		/// Angle between the observer and the \f$z\f$ axis, \f$\theta\f$.
		const Type theta_;
		/// \f$\sin\theta\f$
		const Type sin_theta_;
		/// \f$\cos\theta\f$
		const Type cos_theta_;
		/// Rotational angle of the coordiante of the view, \f$\iota\f$.
		const Type iota_;
		/// \f$\sin\iota\f$
		const Type sin_iota_;
		/// \f$\cos\iota\f$
		const Type cos_iota_;
		/// Position and velocity of the observer.
		std::array<Type, 8> position_;
		/// time limit for the integration of photons.
		const Type t_final_;

	  public:
		/**
		 * @brief Constructor
		 *
		 * @param r Distance to the black hole
		 * @param theta Angle between the observer and the \f$z\f$ axis, \f$\theta\f$.
		 * @param iota Rotational angle of the coordiante of the view, \f$\iota\f$.
		 */
		View(std::shared_ptr<Metric<Type>> metric, Type r, Type theta, Type iota, Type v_alpha = 0.0, Type v_beta = 0.0) : metric_(metric), r_(r), r2_(boost::algorithm::power(r, 2)), theta_(theta), sin_theta_(sin(theta)), cos_theta_(cos(theta)), iota_(iota), sin_iota_(sin(iota)), cos_iota_(cos(iota)), t_final_(-2e4) {
			position_[0] = 0.0;
			position_[1] = r_;
			position_[2] = theta_;
			position_[4] = 1.0;
			position_[5] = 0.0;
			Type v_y = v_alpha * cos_iota_ + v_beta * sin_iota_, v_xz = v_beta * cos_iota_ - v_alpha * sin_iota_;
			if (sin_theta_ < GSL_SQRT_DBL_EPSILON) {
				const Type v = gsl_hypot(v_y, v_xz);
				if (theta_ < M_PI_2) {
					position_[3] = atan2(-v_y, v_xz);
					position_[6] = v / r_;
				} else {
					position_[3] = atan2(-v_y, -v_xz);
					position_[6] = -v / r_;
				}
				position_[7] = 0.0;
			} else {
				position_[3] = 0.0;
				position_[6] = v_xz / r_;
				position_[7] = -v_y / (r * sin_theta_);
			}
			metric_->NormalizeTimelikeGeodesic(position_.data());
			// record[0] = alpha * cos_iota_ - beta * sin_iota_; // alpha
			// record[1] = beta * cos_iota_ + alpha * sin_iota_; // beta
		}

		/**
		 * @brief Initialize the position and velocity of a trace back photon.
		 *
		 * @param photon 9 dimensional vector, photon[8] is used to store the look back time.
		 * @param alpha x position of the target in the observer's view.
		 * @param beta y position of the target in the observer's view.
		 * @return status
		 */
		int
		InitializePhoton(Type photon[], Type alpha, Type beta) {
			return metric_->InitializePhoton(photon, alpha, beta, r_, r2_, theta_, sin_theta_);
		}

		/**
		 * @brief
		 *
		 * @param position 8 dimensional vector
		 * @param object_time time system of the object
		 * @param record position to save the observational values
		 * @param calculate_magnification switch of the magnification calculation
		 * @param fast_trace switch of the fast trace
		 * @return status
		 */
		template <std::size_t N>
		int Trace(const std::array<Type, 8> &position, TimeSystem object_time, std::array<Type, N> &record, bool calculate_magnification, bool fast_trace = true) {
			Type photon[9], alpha, beta;
			if (fast_trace && metric_->FastTrace(r_, theta_, sin_theta_, cos_theta_, position[1], position[2], position[3], alpha, beta, photon) == Status::SUCCESS) {
				PhotonInformation(position.data(), object_time, record.data(), photon, alpha, beta);
				return calculate_magnification ? Magnification(position.data(), object_time, record[4], photon, record[2]) : Status::SUCCESS;
			}
			const Type r_object = position[1], sin_theta_object = abs(sin(position[2])), cos_theta_object = GSL_SIGN(position[2]) * cos(position[2]), sin_phi_object = sin(position[3]), cos_phi_object = cos(position[3]);
			if (r_object <= 3.) {
				PrintlnWarning("Object orbit radius = {:.6f}", r_object);
				if (r_object < 0)
					return GSL_FAILURE;
			}
			const Type cos_observer_object = sin_theta_ * sin_theta_object * cos_phi_object + cos_theta_ * cos_theta_object;
			const Type alpha_coefficient = sin_theta_object * sin_phi_object, beta_coefficient = cos_theta_object * sin_theta_ - sin_theta_object * cos_phi_object * cos_theta_, sin_observer_object = sqrt(Power2(alpha_coefficient) + Power2(beta_coefficient)), theta_observer_object = acos(cos_observer_object);
			GslBlock collector;
			gsl_vector *alpha_beta_initial_value = collector.VectorCalloc(2);
			if (cos_observer_object == -1.) {
				PrintlnWarning("Object behind black hole, cos(theta) = {:.6f}\n", cos_observer_object);
				gsl_vector_set(alpha_beta_initial_value, 0, 2. * sqrt(r_object));
			} else { // initial guessing
				Type effective_radius;
				if (theta_observer_object < M_PI_2)
					effective_radius = r_object + boost::algorithm::power(theta_observer_object / M_PI_2, 3) / sin_observer_object;
				else
					// b-r*sin(theta)=(b-1.)*(2.*theta/pi-1.)+1.+(b-4.)/pi*sin(theta*2.)
					// b=(r*sin(theta)*M_PI+M_2_PI-2.*theta-8.*sin(theta)*cos(theta))/(2.*(M_PI-theta-sin(theta)*cos(theta)))
					effective_radius = 1. / sin_observer_object + 0.5 * (M_PI * r_object - 6. * cos_observer_object) / (M_PI - theta_observer_object - sin_observer_object * cos_observer_object);
				gsl_vector_set(alpha_beta_initial_value, 0, effective_radius * alpha_coefficient);
				gsl_vector_set(alpha_beta_initial_value, 1, effective_radius * beta_coefficient);
			}
			TraceParameters trace_parameters(metric_, metric_->GetIntegrator(T, HAMILTONIAN), photon, r_, r2_, theta_, sin_theta_, cos_theta_, position[1], sin_theta_object, cos_theta_object, sin_phi_object, cos_phi_object, t_final_);
			gsl_multiroot_function alpha_beta_function{TraceToPlane, 2, &trace_parameters};
			int status;
			MultiFunctionSolver alpha_beta_translation_solver(2, gsl_multiroot_fsolver_sbody_dnewton_translation);
			if (status = alpha_beta_translation_solver.Set(&alpha_beta_function, alpha_beta_initial_value, theta_, sin_theta_, cos_theta_, r_object, sin_theta_object, cos_theta_object, position[3], sin_phi_object, cos_phi_object, true); status != Status::SUCCESS)
				return status;
			if (status = alpha_beta_translation_solver.Solve(r_object * GSL_ROOT3_DBL_EPSILON); status == Status::SUCCESS) {
				alpha = gsl_vector_get(alpha_beta_translation_solver.Root(), 0);
				beta = gsl_vector_get(alpha_beta_translation_solver.Root(), 1);
			} else {
				MultiFunctionSolver alpha_beta_rotation_solver(2, gsl_multiroot_fsolver_sbody_dnewton_rotation);
				if (status = alpha_beta_rotation_solver.Set(&alpha_beta_function, alpha_beta_translation_solver.Root(), theta_, sin_theta_, cos_theta_, r_object, sin_theta_object, cos_theta_object, position[3], sin_phi_object, cos_phi_object, true); status != Status::SUCCESS)
					return status;
				if (status = alpha_beta_rotation_solver.Solve(r_object * GSL_ROOT3_DBL_EPSILON); status == Status::SUCCESS) {
					alpha = gsl_vector_get(alpha_beta_rotation_solver.Root(), 0);
					beta = gsl_vector_get(alpha_beta_rotation_solver.Root(), 1);
				} else {
					// PrintlnWarning("Kerr FastTrace() TRANSLATION failed with status = {}", status);
					MultiFunctionSolver alpha_beta_direction_solver(2, gsl_multiroot_fsolver_sbody_direction);
					if (status = alpha_beta_direction_solver.Set(&alpha_beta_function, alpha_beta_rotation_solver.Root()); status != Status::SUCCESS) {
						// PrintlnError("Kerr Trace() set DIRECTION failed with status = {}", status);
						return status;
					}
					if (status = alpha_beta_direction_solver.Solve(r_object * GSL_ROOT3_DBL_EPSILON); status != Status::SUCCESS) {
						// PrintlnError("Kerr Trace() DIRECTION failed with status = {}", status);
						return status;
					}
					alpha = gsl_vector_get(alpha_beta_direction_solver.Root(), 0);
					beta = gsl_vector_get(alpha_beta_direction_solver.Root(), 1);
				}
			}
			photon[8] -= r_;
			PhotonInformation(position.data(), object_time, record.data(), photon, alpha, beta);
			return calculate_magnification ? Magnification(position.data(), object_time, record[4], photon, record[2]) : Status::SUCCESS;
		}

		static int TraceToPlane(const gsl_vector *alpha_beta, void *params, gsl_vector *delta_apparent_alpha_beta) {
			auto *const param = static_cast<TraceParameters<Type> *>(params);
			Type alpha = gsl_vector_get(alpha_beta, 0), beta = gsl_vector_get(alpha_beta, 1);
			auto integrator = param->integrator;
			integrator->Reset();
			if (!isfinite(alpha) || !isfinite(beta))
				return GSL_ERUNAWAY;
			double *photon = param->photon, last_step_record[9];
			if (int status = param->metric->InitializePhoton(photon, alpha, beta, param->r, param->r2, param->theta_obs, param->sin_theta_obs); status != Status::SUCCESS)
				return status;
			param->metric->LagrangianToHamiltonian(photon);
			int status = 0, fixed = 0;
			Type h = -0.01 * param->r, last_h;
#ifdef GSL_RANGE_CHECK_OFF
			auto t_start = std::chrono::steady_clock::now();
#endif
			while (status <= 0) {
				if (status == GSL_FAILURE) {
					h *= 0.5;
					fixed = 0;
				}
				std::copy(photon, photon + 9, last_step_record);
				last_h = h;
				if (fixed)
					status = integrator->ApplyFixedStep(photon + 8, h, photon);
				else
					status = integrator->ApplyStep(photon + 8, param->t_final, &h, photon);
				if (const Type cos_observer_object_photon = (param->x_obj - photon[1] * abs(sin(photon[2])) * cos(photon[3])) * param->sin_theta_obs + (param->z_obj - photon[1] * GSL_SIGN(photon[2]) * cos(photon[2])) * param->cos_theta_obs; cos_observer_object_photon > param->r_obj * GSL_SQRT_DBL_EPSILON) {
					// photon goes through the plane of the object
					std::copy(last_step_record, last_step_record + 9, photon);
					h = last_h * 0.3;
					fixed = 1;
					integrator->Reset();
				} else if (cos_observer_object_photon >= 0) {
					// photon in the same plane with the object
					gsl_vector_set(delta_apparent_alpha_beta, 0, param->y_obj - photon[1] * abs(sin(photon[2])) * sin(photon[3]));
					gsl_vector_set(delta_apparent_alpha_beta, 1, (param->z_obj - photon[1] * GSL_SIGN(photon[2]) * cos(photon[2])) * param->sin_theta_obs - (param->x_obj - photon[1] * abs(sin(photon[2])) * cos(photon[3])) * param->cos_theta_obs);
					return Status::SUCCESS;
				}
				// photon fall into the BH
				if (abs(photon[5]) * GSL_SQRT_DBL_EPSILON > 1.) {
					gsl_vector_set(delta_apparent_alpha_beta, 0, 1.1);
					return GSL_EDOM;
				}
#ifdef GSL_RANGE_CHECK_OFF
				if (auto now = std::chrono::steady_clock::now(); now - t_start > std::chrono::milliseconds(100)) {
					gsl_vector_set(delta_apparent_alpha_beta, 0, 1.1);
					return GSL_EDOM;
				}
#endif
				// photon has not reached the plane of the object
				if (photon[8] > param->t_final)
					continue;
				// photon failed to hit the plane, the impact params need to be larger
				gsl_vector_set(delta_apparent_alpha_beta, 0, 1.01);
				return GSL_EDOM;
			}
			PrintlnWarning("View::Trece() status = {}\n", status);
			return status;
		}

		int PhotonInformation(const Type position[], TimeSystem object_time, Type record[], const Type photon[], Type alpha, Type beta) {
			record[0] = alpha * cos_iota_ - beta * sin_iota_;				 // alpha
			record[1] = beta * cos_iota_ + alpha * sin_iota_;				 // beta
			record[2] = metric_->Redshift(position, photon, object_time, T); // redshift
			record[3] = photon[8];											 // look back time
			return Status::SUCCESS;
		}

		/**
		 * @brief Calculate the magnification of the object.
		 *
		 * @param position 8 dimensional vector
		 * @param object_time time system of the object
		 * @param magnification position to save the magnification
		 * @param photon 8 dimensional vector of the photon traced to the object
		 * @param redshift redshift of the photon
		 * @return int
		 */
		int Magnification(const Type position[], TimeSystem object_time, Type &magnification, const Type photon[], Type redshift) {
			std::unique_ptr<Integrator> integrator = metric_->GetIntegrator(T, HAMILTONIAN);
			Type forward_photon[8], h = 1., t = 0.;
			std::array<Type, 3> cone_record[SAMPLE_NUMBER], local_cone_record[SAMPLE_NUMBER], center_photon_velocity;
			auto forward_photon_velocity_view = gsl_vector_view_array(forward_photon + 4, 4);
			std::copy(photon, photon + 8, forward_photon);
			metric_->NormalizeNullGeodesic(forward_photon);
			metric_->LagrangianToHamiltonian(forward_photon);
			integrator->Reset();
			if (int status = integrator->Apply(&t, 1000., &h, forward_photon); status > 0) {
				PrintlnError("View::Magnification() status = {}", status);
				return status;
			}
			metric_->HamiltonianToLagrangian(forward_photon);
			SphericalToCartesian(forward_photon);
			std::copy(forward_photon + 5, forward_photon + 8, center_photon_velocity.begin());
			cblas_dscal(3, 1. / Norm(center_photon_velocity), center_photon_velocity.data(), 1);
			GslBlock collector;
			auto gmunu = collector.MatrixAlloc(4, 4);			 // object local metric tensor
			auto coordinate = collector.MatrixAlloc(4, 4);		 // object local inertial coordinate frame
			auto coordinate_gmunu = collector.MatrixAlloc(4, 4); // object local inertial frame measured by observer
			auto permutation = collector.PermutationAlloc(4);	 // permutation used in the LU decomposition
			auto photon_transform = collector.VectorAlloc(4);	 // photon in TimeSystem TAU
			auto photon_in_object_frame_cartesian = collector.VectorAlloc(4);
			Type photon_in_object_frame_spherical[4];
			metric_->MetricTensor(position, gmunu);
			metric_->LocalInertialFrame(position, object_time, coordinate);
			gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu, coordinate, 0., coordinate_gmunu);
			gsl_vector_set(photon_transform, 0, 1.);
			std::copy(photon + 5, photon + 8, gsl_vector_ptr(photon_transform, 1));
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
			const Type local_redshift = photon_in_object_frame_cartesian->data[0] / photon_in_static_frame_cartesian->data[0];
			gsl_vector_scale(photon_in_static_frame_cartesian, 1. / gsl_vector_get(photon_in_static_frame_cartesian, 0));
#endif
			int signum;
			gsl_linalg_LU_decomp(coordinate_gmunu, permutation, &signum);
			for (int i = 0; i < SAMPLE_NUMBER; ++i) {
				const Type angle = i * ANGLE_INTERVAL;
				std::copy(photon, photon + 4, forward_photon);
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
				cblas_dscal(3, 1. / cblas_ddot(4, forward_photon + 4, 1, coordinate_static_gmunu->data, 1), local_cone_record[i].data(), 1);
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
				std::copy(forward_photon + 5, forward_photon + 8, cone_record[i].begin());
				cblas_dscal(3, 1. / Norm(cone_record[i]), cone_record[i].data(), 1);
			}
			Type cone_solid_angle = TriangleArea(center_photon_velocity, cone_record[0], cone_record[SAMPLE_NUMBER - 1]);
			std::array<Type, 3> photon_in_static_frame_cartesian_position;
#ifndef GSL_RANGE_CHECK_OFF
			std::copy(photon_in_static_frame_cartesian->data + 1, photon_in_static_frame_cartesian->data + 4, photon_in_static_frame_cartesian_position.begin());
			Type cone_local_solid_angle = TriangleArea(photon_in_static_frame_cartesian_position, local_cone_record[0], local_cone_record[SAMPLE_NUMBER - 1]);
#endif
			for (int i = 1; i < SAMPLE_NUMBER; ++i) {
				cone_solid_angle += TriangleArea(center_photon_velocity, cone_record[0], cone_record[SAMPLE_NUMBER - 1]);
#ifndef GSL_RANGE_CHECK_OFF
				cone_local_solid_angle += TriangleArea(photon_in_static_frame_cartesian_position, local_cone_record[i], local_cone_record[i - 1]);
#endif
			}
			magnification = EPSILON_POLYGON_AREA / (cone_solid_angle * redshift);
			return Status::SUCCESS;
		}

		/**
		 * @brief Calculate the radius of the black hole shadow at different direction, saved in the file.
		 *
		 * @param file_name file name.
		 * @return status
		 */
		int Shadow(std::string file_name, std::optional<ProgressBar> &bars) {
			NumPy record(file_name, {2});
			Type h, rin = 2., rout = 10., rmid = 6., photon[10];
			std::unique_ptr<Integrator> integrator = metric_->GetIntegrator(T, HAMILTONIAN);
			indicators::BlockProgressBar bar{
				indicators::option::ShowElapsedTime{true},
				indicators::option::ShowRemainingTime{true},
				indicators::option::ForegroundColor{indicators::Color(4)},
				indicators::option::FontStyles{
					std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}};
			bar.set_option(indicators::option::MaxProgress(SAMPLE_NUMBER));
			bar.set_option(indicators::option::PrefixText("? Shadow"));
			int progress_bar_index = bars.has_value() ? bars->push_back(bar) : -1;
			for (int i = 0; i < SAMPLE_NUMBER; ++i) {
				const Type angle = i * ANGLE_INTERVAL, sin_angle = sin(angle), cos_angle = cos(angle);
				int status = 0;
				while (rout - rin > GSL_SQRT_DBL_EPSILON * (rin + rout)) {
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
				if (progress_bar_index >= 0)
					(*bars)[progress_bar_index].tick();
			}
			if (progress_bar_index >= 0)
				bars->SetComplete(progress_bar_index, "! Shadow");
			return Status::SUCCESS;
		}

		int OmegaTest(std::optional<ProgressBar> &bars = std::nullopt) {
			std::unique_ptr<Integrator> integrator = metric_->GetIntegrator(T, HAMILTONIAN);
			Type position_[8] = {0., 3., M_PI_4, 0., 0., 0., 0., 0.};
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
			Type rec[100][3], center_ph[3];
			Type ph[9], area, h, time_limit = 0.;
			gsl_vector_view ph_view = gsl_vector_view_array(ph + 4, 4);
			NumPy cone_record("Omega_record", {1});
			bars.value()[0].set_option(indicators::option::MaxProgress(90));
			for (Type angle = 90; angle > 0; angle -= 1) {
				const Type sina = sin(angle / 180. * M_PI), cosa = cos(angle / 180. * M_PI);
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
				const Type sin_theta = GSL_SIGN(ph[2]) * sin(ph[2]), cos_theta = GSL_SIGN(ph[2]) * cos(ph[2]), sin_phi = sin(ph[3]), cos_phi = cos(ph[3]);
				center_ph[0] = ph[5] * sin_theta * cos_phi + ph[1] * (cos_theta * cos_phi * ph[6] - sin_theta * sin_phi * ph[7]);
				center_ph[1] = ph[5] * sin_theta * sin_phi + ph[1] * (cos_theta * sin_phi * ph[6] + sin_theta * cos_phi * ph[7]);
				center_ph[2] = ph[5] * cos_theta - ph[1] * sin_theta * ph[6];
				const Type vph_norm = Norm(center_ph);
				for (int j = 0; j < 3; ++j)
					center_ph[j] /= vph_norm;
				for (int i = 0; i < 100; ++i) {
					const Type angle_i = i * ANGLE_INTERVAL, sinai = sin(angle_i), cosai = cos(angle_i);
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
					const Type sin_theta = GSL_SIGN(ph[2]) * sin(ph[2]), cos_theta = GSL_SIGN(ph[2]) * cos(ph[2]), sin_phi = sin(ph[3]), cos_phi = cos(ph[3]);
					rec[i][0] = ph[5] * sin_theta * cos_phi + ph[1] * (cos_theta * cos_phi * ph[6] - sin_theta * sin_phi * ph[7]);
					rec[i][1] = ph[5] * sin_theta * sin_phi + ph[1] * (cos_theta * sin_phi * ph[6] + sin_theta * cos_phi * ph[7]);
					rec[i][2] = ph[5] * cos_theta - ph[1] * sin_theta * ph[6];
					const Type vph_norm = Norm(rec[i]);
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
				cone_record.Save({abs(area) / (M_2PI * Power2(GSL_SQRT_DBL_EPSILON))});
				bars.value()[0].tick();
			}
			return Status::SUCCESS;
		}
	};

	template <typename Type>
	class Camera : public View<Type> {
	  protected:
		const size_t pixel_;
		const Type half_angle_;
		std::vector<std::array<Type, 10>> initials_;
		std::vector<std::vector<Type>> screen_;

	  public:
		/**
		 * @brief Construct a new camera object
		 *
		 * @param pixel
		 * @param half_angle
		 * @param r
		 * @param theta
		 * @param file_name
		 */
		Camera(std::shared_ptr<Metric<Type>> metric, size_t pixel, Type half_angle, Type r, Type theta, Type iota) : View<Type>(metric, r, theta, iota), pixel_(pixel), half_angle_(half_angle) {
			screen_ = std::vector<std::vector<Type>>(pixel, std::vector<Type>(pixel));
			initials_ = std::vector<std::array<Type, 10>>(pixel * pixel);
			const Type pixel_size = 2. * half_angle * r / pixel, t1 = r + 100.;
			for (size_t i = 0; i < pixel; ++i)
				for (size_t j = 0; j < pixel; ++j)
					this->InitializePhoton(initials_[i * pixel + j].data(), pixel_size * (i - 0.5 * pixel + 0.5), pixel_size * (j - 0.5 * pixel + 0.5));
#pragma omp parallel for
			for (int p = pixel * pixel - 1; p >= 0; --p) {
				int status = 0;
				initials_[p][9] = -1.;
				std::unique_ptr<Integrator> integrator = this->metric_->GetIntegrator(T, HAMILTONIAN);
				this->metric_->NormalizeNullGeodesic(initials_[p].data(), 1.);
				this->metric_->LagrangianToHamiltonian(initials_[p].data());
				while (status <= 0 && initials_[p][8] > t1 && initials_[p][1] > 100)
					status = integrator->Apply(initials_[p].data() + 8, t1, initials_[p].data() + 9, initials_[p].data());
				if (status > 0)
					PrintlnWarning("Camera::Camera() status = {}", status);
			}
		}

		/**
		 * @brief
		 *
		 * @return int
		 */
		int Trace() {
			const Type t1 = -1000.;
#pragma omp parallel for
			for (int p = pixel_ * pixel_ - 1; p >= 0; --p) {
				int i = p / pixel_;
				int j = p - i * pixel_;
				Type ph[10], last[10];
				int status = 0;
				std::copy(initials_[p].begin(), initials_[p].end(), ph);
				std::unique_ptr<Integrator> integrator = this->metric_->GetIntegrator(T, HAMILTONIAN);
				while (status <= 0 && ph[8] > t1) {
					std::copy(ph, ph + 10, last);
					status = integrator->Apply(ph + 8, t1, ph + 9, ph);
					for (auto objP : Object<Type>::object_list_)
						if (objP->Hit(ph, last))
							screen_[i][j] = objP->Redshift(ph, T); // FIXME: if multi objects
					if (screen_[i][j] > GSL_SQRT_DBL_EPSILON)
						break;
				}
				if (status > 0)
					PrintlnWarning("Camera::Trace() status = {}", status);
			}
			return Status::SUCCESS;
		}

		/**
		 * @brief
		 *
		 * @return int
		 */
		int Lens(std::optional<ProgressBar> &bars) {
			const Type t1 = -1000. * this->r_, pixelPerAngle = 0.5 * pixel_ / half_angle_;
			NumPy rec("lens", {2});
			indicators::BlockProgressBar bar{
				indicators::option::ForegroundColor{indicators::Color(4)},
				indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}},
				indicators::option::MaxProgress(pixel_ * pixel_),
				indicators::option::PrefixText("? Lens"),
				indicators::option::ShowElapsedTime{true},
				indicators::option::ShowRemainingTime{true}};
			int progress_bar_index = bars.has_value() ? bars->push_back(bar) : -1;
			for (size_t i = 0; i < pixel_; ++i)
				for (size_t j = 0; j < pixel_; ++j) {
					Type ph[10];
					int status = 0;
					std::copy(initials_[i * pixel_ + j].begin(), initials_[i * pixel_ + j].end(), ph);
					std::unique_ptr<Integrator> integrator = this->metric_->GetIntegrator(T, HAMILTONIAN);
					while (status <= 0 && ph[8] > t1 && ph[1] > 3. && ph[1] < 3.e2)
						status = integrator->Apply(ph + 8, t1, ph + 9, ph);
					if (status > 0)
						PrintlnWarning("Camera::Lens() status = {}", status);
					if (ph[1] <= 3.)
						rec.Save({NAN, NAN});
					else {
						this->metric_->HamiltonianToLagrangian(ph);
						if (ph[2] < 0)
							ph[2] += M_PI;
						SphericalToCartesian(ph);
						rec.Save({ph[6] * pixelPerAngle, (ph[7] * this->sin_theta_ - ph[5] * this->cos_theta_) * pixelPerAngle});
					}
					if (progress_bar_index >= 0)
						(*bars)[progress_bar_index].tick();
				}
			if (progress_bar_index >= 0)
				bars->SetComplete(progress_bar_index, "! lens");
			return Status::SUCCESS;
		}

		/**
		 * @brief
		 *
		 * @return int
		 */
		int Save(std::string file_name) {
			NumPy record(file_name, {2});
			for (const std::vector<Type> &line : screen_)
				record.Save(line);
			return Status::SUCCESS;
		}
	};
} // namespace SBody

#endif
