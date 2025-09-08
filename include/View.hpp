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

#ifndef SBODY_VIEW_HPP
#define SBODY_VIEW_HPP

#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <boost/algorithm/algorithm.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fmt/core.h>

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
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> integration_system;
		boost::numeric::ublas::bounded_vector<Type, 8> &photon;
		Type &photon_time;
		const Type r, r2;
		const Type theta_obs, sin_theta_obs, cos_theta_obs;
		const Type r_obj, x_obj, y_obj, z_obj;
		const Type t_final;
		TraceParameters(std::shared_ptr<Metric<Type>> metric, std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> integration_system, boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time, Type r, Type r2, Type theta_obs, Type sin_theta_obs, Type cos_theta_obs, Type r_obj, Type sin_theta_obj, Type cos_theta_obj, Type sin_phi_obj, Type cos_phi_obj, Type t_final) : metric(metric), integration_system(integration_system), photon(photon), photon_time(photon_time), r(r), r2(r2), theta_obs(theta_obs), sin_theta_obs(sin_theta_obs), cos_theta_obs(cos_theta_obs), r_obj(r_obj), x_obj(r_obj * sin_theta_obj * cos_phi_obj), y_obj(r_obj * sin_theta_obj * sin_phi_obj), z_obj(r_obj * cos_theta_obj), t_final(t_final) {}
	};
	template <typename Type>
	class View {
	  protected:
		std::shared_ptr<Metric<Type>> metric_;
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> integration_system_;
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
		boost::numeric::ublas::bounded_vector<Type, 8> position_;
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
		View(std::shared_ptr<Metric<Type>> metric, Type r, Type theta, Type iota, Type v_alpha = 0.0, Type v_beta = 0.0) : metric_(metric), integration_system_(metric->GetIntegrationSystem(T, HAMILTONIAN)), r_(r), r2_(r * r), theta_(theta), sin_theta_(std::sin(theta)), cos_theta_(std::cos(theta)), iota_(iota), sin_iota_(std::sin(iota)), cos_iota_(std::cos(iota)), t_final_(-2e4) {
			position_[0] = 0.0;
			position_[1] = r_;
			position_[2] = theta_;
			position_[4] = 1.0;
			position_[5] = 0.0;
			Type v_y = v_alpha * cos_iota_ + v_beta * sin_iota_, v_xz = v_beta * cos_iota_ - v_alpha * sin_iota_;
			if (sin_theta_ < boost::math::tools::root_epsilon<Type>()) {
				const Type v = std::hypot(v_y, v_xz);
				if (theta_ < boost::math::constants::half_pi<Type>()) {
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
			metric_->NormalizeTimelikeGeodesic(position_);
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
		int InitializePhoton(boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time, Type alpha, Type beta) {
			return metric_->InitializePhoton(photon, photon_time, alpha, beta, r_, r2_, theta_, sin_theta_);
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
		int Trace(const boost::numeric::ublas::bounded_vector<Type, 8> &position, TimeSystem object_time, boost::numeric::ublas::bounded_vector<Type, 5> &record, bool calculate_magnification, bool fast_trace = true) {
			namespace ublas = boost::numeric::ublas;
			Type alpha, beta;
			ublas::bounded_vector<Type, 8> photon;
			Type photon_time;
			if (fast_trace && metric_->FastTrace(r_, theta_, sin_theta_, cos_theta_, position[1], position[2], position[3], alpha, beta, photon, photon_time) == Status::SUCCESS) {
				PhotonInformation(position, object_time, record, photon, photon_time, alpha, beta);
				return calculate_magnification ? Magnification(position, object_time, record[4], photon, record[2]) : Status::SUCCESS;
			}
			const Type r_object = position[1], sin_theta_object = SinTheta(position[2]), cos_theta_object = CosTheta(position[2]), sin_phi_object = std::sin(position[3]), cos_phi_object = std::cos(position[3]);
			if (r_object <= 3.) {
				PrintlnWarning("Object orbit radius = {:.6f}", r_object);
				if (r_object < 0)
					return Status::FALLS_INTO_BLACK_HOLE;
			}
			const Type cos_observer_object = sin_theta_ * sin_theta_object * cos_phi_object + cos_theta_ * cos_theta_object;
			const Type alpha_coefficient = sin_theta_object * sin_phi_object, beta_coefficient = cos_theta_object * sin_theta_ - sin_theta_object * cos_phi_object * cos_theta_, sin_observer_object = std::hypot(alpha_coefficient, beta_coefficient), theta_observer_object = std::acos(cos_observer_object);
			ublas::bounded_vector<Type, 2> alpha_beta_initial_value;
			if (cos_observer_object == -1.) {
				PrintlnWarning("Object behind black hole, cos(theta) = {:.6f}\n", cos_observer_object);
				alpha_beta_initial_value(0) = 2. * std::sqrt(r_object);
			} else { // initial guessing
				Type effective_radius;
				if (theta_observer_object < boost::math::constants::half_pi<Type>())
					effective_radius = r_object + Power3(theta_observer_object / boost::math::constants::half_pi<Type>()) / sin_observer_object;
				else
					// b-r*sin(theta)=(b-1.)*(2.*theta/pi-1.)+1.+(b-4.)/pi*sin(theta*2.)
					// b=(r*sin(theta)*boost::math::constants::pi<Type>()+M_2_PI-2.*theta-8.*sin(theta)*cos(theta))/(2.*(boost::math::constants::pi<Type>()t::math::constants::pi<Type>()-theta-sin(theta)*cos(theta)))
					effective_radius = 1. / sin_observer_object + (boost::math::constants::half_pi<Type>() * r_object - 3. * cos_observer_object) / (boost::math::constants::pi<Type>() - theta_observer_object - sin_observer_object * cos_observer_object);
				alpha_beta_initial_value(0) = effective_radius * alpha_coefficient;
				alpha_beta_initial_value(1) = effective_radius * beta_coefficient;
			}
			TraceParameters trace_parameters(metric_, metric_->GetIntegrationSystem(T, HAMILTONIAN), photon, photon_time, r_, r2_, theta_, sin_theta_, cos_theta_, position[1], sin_theta_object, cos_theta_object, sin_phi_object, cos_phi_object, t_final_);
			int status;
			DNewtonTranslationMultiFunctionSolver<Type, 2, TraceParameters<Type>> alpha_beta_translation_solver(TraceToPlane, trace_parameters);
			if (status = alpha_beta_translation_solver.Set(alpha_beta_initial_value, theta_, sin_theta_, cos_theta_, r_object, sin_theta_object, cos_theta_object, position[3], sin_phi_object, cos_phi_object, true); status == Status::SUCCESS)
				if (status = alpha_beta_translation_solver.Solve(r_object * boost::math::tools::root_epsilon<Type>()); status == Status::SUCCESS) {
					alpha = alpha_beta_translation_solver.Root()(0);
					beta = alpha_beta_translation_solver.Root()(1);
					photon_time -= r_;
					PhotonInformation(position, object_time, record, photon, photon_time, alpha, beta);
					return calculate_magnification ? Magnification(position, object_time, record[4], photon, record[2]) : Status::SUCCESS;
				}
			DNewtonRotationMultiFunctionSolver<Type, 2, TraceParameters<Type>> alpha_beta_rotation_solver(TraceToPlane, trace_parameters);
			if (status = alpha_beta_rotation_solver.Set(alpha_beta_translation_solver.Root(), theta_, sin_theta_, cos_theta_, r_object, sin_theta_object, cos_theta_object, position[3], sin_phi_object, cos_phi_object, true); status == Status::SUCCESS)
				if (status = alpha_beta_rotation_solver.Solve(r_object * boost::math::tools::root_epsilon<Type>()); status == Status::SUCCESS) {
					alpha = alpha_beta_rotation_solver.Root()(0);
					beta = alpha_beta_rotation_solver.Root()(1);
					photon_time -= r_;
					PhotonInformation(position, object_time, record, photon, photon_time, alpha, beta);
					return calculate_magnification ? Magnification(position, object_time, record[4], photon, record[2]) : Status::SUCCESS;
				}
			// PrintlnWarning("Kerr FastTrace() TRANSLATION failed with status = {}", status);
			DirectionMultiFunctionSolver<Type, 2, TraceParameters<Type>> alpha_beta_direction_solver(TraceToPlane, trace_parameters);
			if (status = alpha_beta_direction_solver.Set(alpha_beta_rotation_solver.Root()); status == Status::SUCCESS)
				// PrintlnError("Kerr Trace() set DIRECTION failed with status = {}", status);
				if (status = alpha_beta_direction_solver.Solve(r_object * boost::math::tools::root_epsilon<Type>()); status == Status::SUCCESS) {
					// PrintlnError("Kerr Trace() DIRECTION failed with status = {}", status);
					alpha = alpha_beta_direction_solver.Root()(0);
					beta = alpha_beta_direction_solver.Root()(1);
					photon_time -= r_;
					PhotonInformation(position, object_time, record, photon, photon_time, alpha, beta);
					return calculate_magnification ? Magnification(position, object_time, record[4], photon, record[2]) : Status::SUCCESS;
				}
			return status;
		}

		static int TraceToPlane(const boost::numeric::ublas::bounded_vector<Type, 2> &alpha_beta, boost::numeric::ublas::bounded_vector<Type, 2> &delta_apparent_alpha_beta, TraceParameters<Type> &params) {
			Type alpha = alpha_beta(0), beta = alpha_beta(1);
			if (!std::isfinite(alpha) || !std::isfinite(beta))
				return Status::NUMERIC_ERROR;
			Type photon_time;
			if (int status = params.metric->InitializePhoton(params.photon, photon_time, alpha, beta, params.r, params.r2, params.theta_obs, params.sin_theta_obs); status != Status::SUCCESS)
				return status;
			params.metric->LagrangianToHamiltonian(params.photon);
			Type h = -0.01 * params.r, last_h;
#ifdef SBODY_RELEASE
			auto wall_clock_start_time = std::chrono::steady_clock::now();
#endif
			auto integration_stepper = boost::numeric::odeint::make_dense_output(absolute_accuracy, relative_accuracy, boost::numeric::odeint::runge_kutta_dopri5<boost::numeric::ublas::bounded_vector<Type, 8>>());
			integration_stepper.initialize(params.photon, photon_time, h);
			Type photon_last_theta = params.photon[2], photon_hit_time, photon_hit_time_lower_bound, photon_hit_time_upper_bound;
			while (true) {
				try {
					integration_stepper.do_step(params.integration_system);
				} catch (const std::exception &e) {
					return Status::FAILURE;
				}
				const boost::numeric::ublas::bounded_vector<Type, 8> &photon_state = integration_stepper.current_state();
				if (const Type cos_observer_object_photon = (params.x_obj - photon_state[1] * SinTheta(photon_state[2]) * std::cos(photon_state[3])) * params.sin_theta_obs + (params.z_obj - photon_state[1] * CosTheta(photon_state[2])) * params.cos_theta_obs; cos_observer_object_photon >= 0.) {
					// photon goes through the plane of the object
					photon_hit_time_lower_bound = integration_stepper.current_time();
					boost::numeric::ublas::bounded_vector<Type, 8> photon_hit_state;
					while (photon_hit_time_upper_bound - photon_hit_time_lower_bound > boost::math::tools::root_epsilon<Type>()) {
						photon_hit_time = 0.5 * (photon_hit_time_upper_bound + photon_hit_time_lower_bound);
						integration_stepper.calc_state(photon_hit_time, photon_hit_state);
						if ((params.x_obj - photon_hit_state[1] * SinTheta(photon_hit_state[2]) * std::cos(photon_hit_state[3])) * params.sin_theta_obs + (params.z_obj - photon_hit_state[1] * CosTheta(photon_hit_state[2])) * params.cos_theta_obs > 0.)
							photon_hit_time_lower_bound = photon_hit_time;
						else
							photon_hit_time_upper_bound = photon_hit_time;
					}
					// photon in the same plane with the object
					params.photon_time = photon_hit_time;
					delta_apparent_alpha_beta(0) = params.y_obj - photon_hit_state[1] * SinTheta(photon_hit_state[2]) * std::sin(photon_hit_state[3]);
					delta_apparent_alpha_beta(1) = (params.z_obj - photon_hit_state[1] * CosTheta(photon_hit_state[2])) * params.sin_theta_obs - (params.x_obj - photon_hit_state[1] * SinTheta(photon_hit_state[2]) * std::cos(photon_hit_state[3])) * params.cos_theta_obs;
					return Status::SUCCESS;
				}
				// photon fall into the BH
				if (std::abs(photon_state[5]) * boost::math::tools::root_epsilon<Type>() > 1.) {
					delta_apparent_alpha_beta(0) = 1.1;
					return Status::SCALING_REQUIRED;
				}
#ifdef SBODY_RELEASE
				if (std::chrono::steady_clock::now() - wall_clock_start_time > std::chrono::milliseconds(100)) {
					delta_apparent_alpha_beta(0) = 1.1;
					return Status::COMPUTATION_TIMEOUT;
				}
#endif
				// photon has not reached the plane of the object
				if (integration_stepper.current_time() > params.t_final) {
					photon_hit_time_upper_bound = integration_stepper.current_time();
					DenseStepperMapTheta(integration_stepper, photon_last_theta);
					continue;
				}
				// photon failed to hit the plane, the impact params need to be larger
				delta_apparent_alpha_beta(0) = 1.01;
				return Status::SCALING_REQUIRED;
			}
		}

		int PhotonInformation(const boost::numeric::ublas::bounded_vector<Type, 8> &position, TimeSystem object_time, boost::numeric::ublas::bounded_vector<Type, 5> &record, const boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type photon_time, Type alpha, Type beta) {
			record(0) = alpha * cos_iota_ - beta * sin_iota_;				 // alpha
			record(1) = beta * cos_iota_ + alpha * sin_iota_;				 // beta
			record(2) = metric_->Redshift(position, photon, object_time, T); // redshift
			record(3) = photon_time;										 // look back time
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
		int Magnification(const boost::numeric::ublas::bounded_vector<Type, 8> &position, TimeSystem object_time, Type &magnification, const boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type redshift) {
			namespace ublas = boost::numeric::ublas;
			boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_dopri5<ublas::bounded_vector<Type, 8>>> integration_stepper;
			ublas::bounded_vector<Type, 8> forward_photon;
			ublas::bounded_vector<Type, 3> cone_record[SAMPLE_NUMBER], local_cone_record[SAMPLE_NUMBER], center_photon_velocity;
			auto forward_photon_velocity = ublas::vector_range(forward_photon, ublas::range(4, 8));
			forward_photon = photon;
			metric_->NormalizeNullGeodesic(forward_photon);
			metric_->LagrangianToHamiltonian(forward_photon);
			try {
				boost::numeric::odeint::integrate_adaptive(integration_stepper, integration_system_, forward_photon, 0.0, 1000., 1.0);
			} catch (const std::exception &e) {
				return Status::FAILURE;
			}
			metric_->HamiltonianToLagrangian(forward_photon);
			SphericalToCartesian(forward_photon);
			std::copy(forward_photon.begin() + 5, forward_photon.end(), center_photon_velocity.begin());
			center_photon_velocity /= ublas::norm_2(center_photon_velocity);
			ublas::bounded_matrix<Type, 4, 4> gmunu;			   // object local metric tensor
			ublas::bounded_matrix<Type, 4, 4> coordinate;		   // object local inertial coordinate frame
			ublas::bounded_matrix<Type, 4, 4> coordinate_gmunu;	   // object local inertial coordinate frame
			ublas::permutation_matrix<std::size_t> permutation(4); // permutation used in the LU decomposition
			ublas::bounded_vector<Type, 4> photon_transform;	   // photon in TimeSystem TAU
			ublas::bounded_vector<Type, 4> photon_in_object_frame_cartesian, photon_in_object_frame_spherical;
			metric_->MetricTensor(position, gmunu);
			metric_->LocalInertialFrame(position, object_time, coordinate);
			coordinate_gmunu = ublas::prod(gmunu, coordinate);									// gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu_gsl, coordinate_gsl, 0., coordinate_gmunu_gsl);
			std::copy(photon.begin() + 5, photon.end(), photon_transform.begin() + 1);			// std::copy(photon.begin() + 5, photon.end(), gsl_vector_ptr(photon_transform_gsl, 1));
			photon_transform(0) = 1.;															// gsl_vector_set(photon_transform_gsl, 0, 1.);
			photon_in_object_frame_cartesian = ublas::prod(coordinate_gmunu, photon_transform); // gsl_blas_dgemv(CblasNoTrans, 1., coordinate_gmunu_gsl, photon_transform_gsl, 0., photon_in_object_frame_cartesian_gsl);
			photon_in_object_frame_cartesian /= photon_in_object_frame_cartesian(0);			// gsl_vector_scale(photon_in_object_frame_cartesian, 1. / gsl_vector_get(photon_in_object_frame_cartesian, 0));
			CartesianToSpherical(photon_in_object_frame_cartesian, photon_in_object_frame_spherical);
#ifndef SBODY_RELEASE
			ublas::bounded_matrix<Type, 4, 4> coordinate_static;	   // object local static frame (only dt/d\tau != 0)
			ublas::bounded_matrix<Type, 4, 4> coordinate_static_gmunu; // object local static frame measured by observer
			ublas::bounded_vector<Type, 4> photon_in_static_frame_cartesian;
			coordinate_static(0, 0) = std::sqrt(-1. / gmunu(0, 0));
			coordinate_static(1, 1) = std::sqrt(1. / gmunu(1, 1));
			coordinate_static(2, 2) = std::sqrt(1. / gmunu(2, 2));
			coordinate_static(3, 3) = std::sqrt(1. / gmunu(3, 3));
			coordinate_static_gmunu = ublas::prod(gmunu, coordinate_static);
			photon_in_static_frame_cartesian = ublas::prod(coordinate_static_gmunu, photon_transform);
			// gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu_gsl, coordinate_static, 0., coordinate_static_gmunu);
			// gsl_blas_dgemv(CblasNoTrans, 1., coordinate_static_gmunu, photon_transform, 0., photon_in_static_frame_cartesian);
			// local_redshift should equal to sqrt(EPSILON_POLYGON_AREA / cone_local_solid_angle),
			// the main error comes from the calculation of the cone_local_solid_angle, due to the limited SAMPLE_NUMBER.
			const Type local_redshift = photon_in_object_frame_cartesian(0) / photon_in_static_frame_cartesian(0);
			photon_in_static_frame_cartesian /= photon_in_static_frame_cartesian(0);
#endif
			int signum;
			ublas::lu_factorize(coordinate_gmunu, permutation);
			for (int i = 0; i < SAMPLE_NUMBER; ++i) {
				const Type angle = i * ANGLE_INTERVAL;
				std::copy(photon.begin(), photon.begin() + 4, forward_photon.begin());
				forward_photon[4] = -1.;
				forward_photon[5] = std::cos(angle) * SIN_EPSILON;
				forward_photon[6] = std::sin(angle) * SIN_EPSILON;
				forward_photon[7] = COS_EPSILON;
				RotateAroundAxis(forward_photon, true, Y, photon_in_object_frame_spherical[2]);
				RotateAroundAxis(forward_photon, true, Z, photon_in_object_frame_spherical[3]);
				ublas::lu_substitute(coordinate_gmunu, permutation, forward_photon_velocity);
				metric_->NormalizeNullGeodesic(forward_photon);
#ifndef SBODY_RELEASE
				// local_cone_record[i][0] = metric_->DotProduct(photon, forward_photon + 4, coordinate_static->data + 4, 4);
				local_cone_record[i][0] = ublas::inner_prod(forward_photon_velocity, ublas::matrix_row(coordinate_static_gmunu, 1));
				// local_cone_record[i][1] = metric_->DotProduct(photon, forward_photon + 4, coordinate_static->data + 8, 4);
				local_cone_record[i][1] = ublas::inner_prod(forward_photon_velocity, ublas::matrix_row(coordinate_static_gmunu, 2));
				// local_cone_record[i][2] = metric_->DotProduct(photon, forward_photon + 4, coordinate_static->data + 12, 4);
				local_cone_record[i][2] = ublas::inner_prod(forward_photon_velocity, ublas::matrix_row(coordinate_static_gmunu, 3));
				local_cone_record[i] /= ublas::inner_prod(forward_photon_velocity, ublas::matrix_row(coordinate_static_gmunu, 0));
#endif
				metric_->LagrangianToHamiltonian(forward_photon);
				try {
					boost::numeric::odeint::integrate_adaptive(integration_stepper, integration_system_, forward_photon, 0.0, 1000., 1.0);
				} catch (const std::exception &e) {
					return Status::FAILURE;
				}
				metric_->HamiltonianToLagrangian(forward_photon);
				SphericalToCartesian(forward_photon);
				std::copy(forward_photon.begin() + 5, forward_photon.end(), cone_record[i].begin());
				cone_record[i] /= ublas::norm_2(cone_record[i]);
			}
			Type cone_solid_angle = TriangleArea(center_photon_velocity, cone_record[0], cone_record[SAMPLE_NUMBER - 1]);
#ifndef SBODY_RELEASE
			boost::numeric::ublas::bounded_vector<Type, 3> photon_in_static_frame_cartesian_position;
			photon_in_static_frame_cartesian_position = ublas::vector_range(photon_in_static_frame_cartesian, ublas::range(1, 4));
			Type cone_local_solid_angle = TriangleArea(photon_in_static_frame_cartesian_position, local_cone_record[0], local_cone_record[SAMPLE_NUMBER - 1]);
#endif
			for (int i = 1; i < SAMPLE_NUMBER; ++i) {
				cone_solid_angle += TriangleArea(center_photon_velocity, cone_record[0], cone_record[SAMPLE_NUMBER - 1]);
#ifndef SBODY_RELEASE
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
			boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_dopri5<boost::numeric::ublas::bounded_vector<Type, 8>>> integration_stepper;
			NumPy record(file_name, {2});
			Type h, rin = 2., rout = 10., rmid = 6.;
			boost::numeric::ublas::bounded_vector<Type, 8> photon;
			Type photon_time;
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
				const Type angle = i * ANGLE_INTERVAL, sin_angle = std::sin(angle), cos_angle = std::cos(angle);
				while (rout - rin > boost::math::tools::root_epsilon<Type>() * (rin + rout)) {
					rmid = 0.5 * (rin + rout);
					InitializePhoton(photon, photon_time, rmid * cos_angle, rmid * sin_angle);
					h = -1.;
					while (photon_time > t_final_) {
						integration_stepper.try_step(integration_system_, photon, photon_time, h);
						if (photon[4] >= 1e6 || photon[5] <= 0 || photon[1] <= 0)
							break;
					}
					if (photon[5] <= 0)
						rout = rmid;
					else
						rin = rmid;
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

		// int OmegaTest(std::optional<ProgressBar> &bars = std::nullopt) {
		// 	boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_dopri5<boost::numeric::ublas::bounded_vector<Type, 8>>> integration_stepper;
		// 	Type position_[8] = {0., 3., M_PI_4, 0., 0., 0., 0., 0.};
		// 	metric_->NormalizeTimelikeGeodesic(position_);
		// 	gsl_matrix *coordinate = gsl_matrix_alloc(4, 4), *gmunu_gsl = gsl_matrix_alloc(4, 4), *coordinate_gmunu = gsl_matrix_alloc(4, 4);
		// 	gsl_permutation *perm = gsl_permutation_alloc(4);
		// 	metric_->MetricTensor(position_, gmunu_gsl);
		// 	gsl_matrix_set_zero(coordinate);
		// 	gsl_matrix_set(coordinate, 0, 0, std::sqrt(-1. / gmunu_gsl->data[0]));
		// 	gsl_matrix_set(coordinate, 1, 1, std::sqrt(1. / gmunu_gsl->data[5]));
		// 	gsl_matrix_set(coordinate, 2, 2, std::sqrt(1. / gmunu_gsl->data[10]));
		// 	gsl_matrix_set(coordinate, 3, 3, std::sqrt(1. / gmunu_gsl->data[15]));
		// 	gsl_matrix_set_zero(coordinate_gmunu);
		// 	gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu_gsl, coordinate, 0., coordinate_gmunu);
		// 	int signum;
		// 	gsl_linalg_LU_decomp(coordinate_gmunu, perm, &signum);
		// 	Type rec[100][3], center_photon[3];
		// 	Type photon_obs_time, area, h, time_limit = 0.;
		// 	boost::numeric::ublas::bounded_vector<Type, 8> photon;
		// 	gsl_vector_view ph_view = gsl_vector_view_array(photon.data() + 4, 4);
		// 	NumPy cone_record("Omega_record", {1});
		// 	bars.value()[0].set_option(indicators::option::MaxProgress(90));
		// 	for (Type angle = 90; angle > 0; angle -= 1) {
		// 		const Type sina = std::sin(angle / 180. * boost::math::constants::pi<Type>()), cosa = std::cos(angle / 180. * boost::math::constants::pi<Type>());
		// 		copy(position_, position_ + 4, photon.data());
		// 		photon(4) = 1.;
		// 		photon(5) = sina;
		// 		photon(6) = cosa;
		// 		photon(7) = 0.;
		// 		gsl_linalg_LU_svx(coordinate_gmunu, perm, &ph_view.vector);
		// 		metric_->NormalizeNullGeodesic(photon);
		// 		if (photon(5) < 0) {
		// 			photon(5) = -photon(5);
		// 			photon(6) = -photon(6);
		// 			photon(7) = -photon(7);
		// 		}
		// 		metric_->LagrangianToHamiltonian(photon);
		// 		photon_obs_time = 0;
		// 		h = 1.;
		// 		int status = 0;
		// 		boost::numeric::odeint::integrate_adaptive(integration_stepper, integration_system_, photon, time_limit, -t_final_, 1.0);
		// 		metric_->HamiltonianToLagrangian(photon);
		// 		const Type sin_theta = SinTheta(photon[2]), cos_theta = CosTheta(photon[2]), sin_phi = std::sin(photon[3]), cos_phi = std::cos(photon[3]);
		// 		center_photon[0] = photon[5] * sin_theta * cos_phi + photon[1] * (cos_theta * cos_phi * photon[6] - sin_theta * sin_phi * photon[7]);
		// 		center_photon[1] = photon[5] * sin_theta * sin_phi + photon[1] * (cos_theta * sin_phi * photon[6] + sin_theta * cos_phi * photon[7]);
		// 		center_photon[2] = photon[5] * cos_theta - photon[1] * sin_theta * photon[6];
		// 		const Type vph_norm = Norm(center_photon);
		// 		for (int j = 0; j < 3; ++j)
		// 			center_photon[j] /= vph_norm;
		// 		for (int i = 0; i < 100; ++i) {
		// 			const Type angle_i = i * ANGLE_INTERVAL, sinai = std::sin(angle_i), cosai = std::cos(angle_i);
		// 			copy(position_, position_ + 4, photon.begin());
		// 			photon[4] = 1.;
		// 			photon[5] = sina - SIN_EPSILON * cosai * cosa;
		// 			photon[6] = cosa + SIN_EPSILON * cosai * sina;
		// 			photon[7] = SIN_EPSILON * sinai;
		// 			gsl_linalg_LU_svx(coordinate_gmunu, perm, &ph_view.vector);
		// 			metric_->NormalizeNullGeodesic(photon);
		// 			if (photon[5] < 0) {
		// 				photon[5] = -photon[5];
		// 				photon[6] = -photon[6];
		// 				photon[7] = -photon[7];
		// 			}
		// 			metric_->LagrangianToHamiltonian(photon);
		// 			photon_obs_time = 0;
		// 			h = 1.;
		// 			status = 0;
		// 			boost::numeric::odeint::integrate_adaptive(integration_stepper, integration_system_, photon, photon_obs_time, time_limit, 1.0);
		// 			metric_->HamiltonianToLagrangian(photon);
		// 			const Type sin_theta = SinTheta(photon[2]), cos_theta = CosTheta(photon[2]), sin_phi = std::sin(photon[3]), cos_phi = std::cos(photon[3]);
		// 			rec[i][0] = photon[5] * sin_theta * cos_phi + photon[1] * (cos_theta * cos_phi * photon[6] - sin_theta * sin_phi * photon[7]);
		// 			rec[i][1] = photon[5] * sin_theta * sin_phi + photon[1] * (cos_theta * sin_phi * photon[6] + sin_theta * cos_phi * photon[7]);
		// 			rec[i][2] = photon[5] * cos_theta - photon[1] * sin_theta * photon[6];
		// 			const Type vph_norm = Norm(rec[i]);
		// 			for (int j = 0; j < 3; ++j)
		// 				rec[i][j] /= vph_norm;
		// 		}
		// 		area = DotCross(center_photon, rec[0], rec[99]);
		// 		for (int i = 1; i < SAMPLE_NUMBER; ++i)
		// 			area += DotCross(center_photon, rec[i], rec[i - 1]);
		// 		if (angle == 13) {
		// 			NumPy cone_record("cone_record13", {3});
		// 			cone_record.Save(center_photon, 3);
		// 			for (int i = 0; i < SAMPLE_NUMBER; ++i)
		// 				cone_record.Save(rec[i], 3);
		// 		}
		// 		cone_record.Save({std::abs(area) / (boost::math::constants::two_pi<Type>() * boost::math::tools::epsilon<Type>())});
		// 		bars.value()[0].tick();
		// 	}
		// 	return Status::SUCCESS;
		// }
	};

	template <typename Type>
	class Camera : public View<Type> {
	  protected:
		const std::size_t pixel_;
		const Type half_angle_;
		std::vector<boost::numeric::ublas::bounded_vector<Type, 8>> initial_positions_;
		std::vector<Type> initial_times_;
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
		Camera(std::shared_ptr<Metric<Type>> metric, std::size_t pixel, Type half_angle, Type r, Type theta, Type iota) : View<Type>(metric, r, theta, iota), pixel_(pixel), half_angle_(half_angle) {
			screen_ = std::vector<std::vector<Type>>(pixel, std::vector<Type>(pixel));
			initial_positions_ = std::vector<boost::numeric::ublas::bounded_vector<Type, 8>>(pixel * pixel);
			const Type pixel_size = 2. * half_angle * r / pixel, t1 = r + 100.;
			boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_dopri5<boost::numeric::ublas::bounded_vector<Type, 8>>> integration_stepper;
			for (std::size_t i = 0; i < pixel; ++i)
				for (std::size_t j = 0; j < pixel; ++j)
					this->InitializePhoton(initial_positions_[i * pixel + j], initial_times_[i * pixel + j], pixel_size * (i - 0.5 * pixel + 0.5), pixel_size * (j - 0.5 * pixel + 0.5));
#pragma omp parallel for
			for (int p = pixel * pixel - 1; p >= 0; --p) {
				this->metric_->NormalizeNullGeodesic(initial_positions_[p], 1.);
				this->metric_->LagrangianToHamiltonian(initial_positions_[p]);
				boost::numeric::odeint::integrate_adaptive(integration_stepper, this->integration_system_, initial_positions_[p], initial_times_[p], t1, -1.0);
			}
		}

		/**
		 * @brief
		 *
		 * @return int
		 */
		int Trace(std::vector<Object<Type> *> &object_list) {
			boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_dopri5<boost::numeric::ublas::bounded_vector<Type, 8>>> integration_stepper;
			const Type t1 = -1000.;
#pragma omp parallel for
			for (int p = pixel_ * pixel_ - 1; p >= 0; --p) {
				int i = p / pixel_;
				int j = p - i * pixel_;
				boost::numeric::ublas::bounded_vector<Type, 8> photon, last;
				Type t = initial_times_[p], h = 1.0;
				int status = 0;
				std::copy(initial_positions_[p].begin(), initial_positions_[p].end(), photon.begin());
				while (status <= 0 && t > t1) {
					std::copy(photon.begin(), photon.end(), last.begin());
					integration_stepper.try_step(this->integration_system_, photon, t, h);
					for (auto object_pointer : object_list)
						if (object_pointer->Hit(photon, last))
							screen_[i][j] = object_pointer->Redshift(photon, T); // FIXME: if multi objects
					if (screen_[i][j] > boost::math::tools::root_epsilon<Type>())
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
			boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_dopri5<boost::numeric::ublas::bounded_vector<Type, 8>>> integration_stepper;
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
			for (std::size_t i = 0; i < pixel_; ++i)
				for (std::size_t j = 0; j < pixel_; ++j) {
					std::size_t idx = i * pixel_ + j;
					boost::numeric::ublas::bounded_vector<Type, 8> photon;
					Type t = initial_times_[idx], h = -1;
					int status = 0;
					std::copy(initial_positions_[idx].begin(), initial_positions_[idx].end(), photon.begin());
					boost::numeric::odeint::integrate_adaptive(integration_stepper, this->integration_system_, photon, t, t1, -1.0);
					if (status > 0)
						PrintlnWarning("Camera::Lens() status = {}", status);
					if (photon[1] <= 3.)
						rec.Save({NAN, NAN});
					else {
						this->metric_->HamiltonianToLagrangian(photon);
						if (photon[2] < 0)
							photon[2] += boost::math::constants::pi<Type>();
						SphericalToCartesian(photon);
						rec.Save({photon[6] * pixelPerAngle, (photon[7] * this->sin_theta_ - photon[5] * this->cos_theta_) * pixelPerAngle});
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
