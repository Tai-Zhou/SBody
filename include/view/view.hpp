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

#pragma once

#include <cmath>
#include <string>
#include <vector>

#include <boost/algorithm/algorithm.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fmt/core.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "file/file.hpp"
#include "infra/infra_context.hpp"
#include "infra/progress.hpp"
#include "infra/status.hpp"
#include "math/constants.hpp"
#include "math/geometry.hpp"
#include "math/solver/system/gauss_newton.hpp"
#include "metric/metric.hpp"

namespace sbody {
template <typename T>
struct TraceParameters {
    const Metric<T> &metric_;
    std::function<void(const boost::numeric::ublas::bounded_vector<T, 8> &, boost::numeric::ublas::bounded_vector<T, 8> &,
                       const T)>
        integration_system_;
    boost::numeric::ublas::bounded_vector<T, 8> &photon_;
    T &photon_time_;
    const T r_, r2_;
    const T theta_obs_, sin_theta_obs_, cos_theta_obs_;
    const T r_obj_, x_obj_, y_obj_, z_obj_;
    const T t_final_;
    const T absolute_accuracy_;
    const T relative_accuracy_;
    TraceParameters(const Metric<T> &metric, boost::numeric::ublas::bounded_vector<T, 8> &photon, T &photon_time, T r, T r2,
                    T theta_obs, T sin_theta_obs, T cos_theta_obs, T r_obj, T sin_theta_obj, T cos_theta_obj, T sin_phi_obj,
                    T cos_phi_obj, T t_final, T absolute_accuracy, T relative_accuracy)
        : metric_(metric), photon_(photon), photon_time_(photon_time), r_(r), r2_(r2), theta_obs_(theta_obs),
          sin_theta_obs_(sin_theta_obs), cos_theta_obs_(cos_theta_obs), r_obj_(r_obj),
          x_obj_(r_obj * sin_theta_obj * cos_phi_obj), y_obj_(r_obj * sin_theta_obj * sin_phi_obj),
          z_obj_(r_obj * cos_theta_obj), t_final_(t_final), absolute_accuracy_(absolute_accuracy),
          relative_accuracy_(relative_accuracy) {
        metric.getIntegrationSystem(integration_system_, TimeSystem::T, DynamicalSystem::HAMILTONIAN);
    }
};
template <typename T>
class View {
  protected:
    const Metric<T> &metric_;
    std::function<void(const boost::numeric::ublas::bounded_vector<T, 8> &, boost::numeric::ublas::bounded_vector<T, 8> &,
                       const T)>
        integration_system_;
    /// Distance to the center black hole
    const T r_;
    /// Square of the distance to the center black hole
    const T r2_;
    /// Angle between the observer and the \f$z\f$ axis, \f$\theta\f$.
    const T theta_;
    /// \f$\sin\theta\f$
    const T sin_theta_;
    /// \f$\cos\theta\f$
    const T cos_theta_;
    /// Rotational angle of the coordiante of the view, \f$\iota\f$.
    const T iota_;
    /// \f$\sin\iota\f$
    const T sin_iota_;
    /// \f$\cos\iota\f$
    const T cos_iota_;
    /// Position and velocity of the observer.
    boost::numeric::ublas::bounded_vector<T, 8> position_;
    /// time limit for the integration of photons.
    const T t_final_;

  public:
    /**
     * @brief Constructor
     *
     * @param r Distance to the black hole
     * @param theta Angle between the observer and the \f$z\f$ axis, \f$\theta\f$.
     * @param iota Rotational angle of the coordiante of the view, \f$\iota\f$.
     */
    View(const Metric<T> &metric, T r, T theta, T iota, T v_alpha = 0.0, T v_beta = 0.0)
        : metric_(metric), r_(r), r2_(r * r), theta_(theta), sin_theta_(std::sin(theta)), cos_theta_(std::cos(theta)),
          iota_(iota), sin_iota_(std::sin(iota)), cos_iota_(std::cos(iota)), t_final_(-2e4) {
        position_[0] = 0.0;
        position_[1] = r_;
        position_[2] = theta_;
        position_[4] = 1.0;
        position_[5] = 0.0;
        T v_y = v_alpha * cos_iota_ + v_beta * sin_iota_, v_xz = v_beta * cos_iota_ - v_alpha * sin_iota_;
        if (sin_theta_ < boost::math::tools::root_epsilon<T>()) {
            const T v = std::hypot(v_y, v_xz);
            if (theta_ < boost::math::constants::half_pi<T>()) {
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
        metric_.normalizeTimelikeGeodesic(position_);
        metric_.getIntegrationSystem(integration_system_, TimeSystem::T, DynamicalSystem::HAMILTONIAN);
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
    int InitializePhoton(boost::numeric::ublas::bounded_vector<T, 8> &photon, T &photon_time, T alpha, T beta) const {
        return metric_.initializePhoton(photon, photon_time, alpha, beta, r_, r2_, theta_, sin_theta_);
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
    int Trace(const boost::numeric::ublas::bounded_vector<T, 8> &position, TimeSystem object_time,
              boost::numeric::ublas::bounded_vector<T, 5> &record, T absolute_accuracy, T relative_accuracy,
              bool calculate_magnification, bool fast_trace_only, const InfraContext &infra_context) const {
        namespace ublas = boost::numeric::ublas;
        T alpha, beta;
        ublas::bounded_vector<T, 8> photon;
        T photon_time;
        int status = StatusCode::SUCCESS;
        try {
            status = metric_.fastTrace(r_, theta_, sin_theta_, cos_theta_, position[1], position[2], position[3], alpha, beta,
                                       photon, photon_time, absolute_accuracy, relative_accuracy, infra_context);
        } catch (...) {
            status = StatusCode::FAILURE;
        }
        if (status == StatusCode::SUCCESS) {
            photonInformation(position, object_time, record, photon, photon_time, alpha, beta);
            return calculate_magnification ? Magnification(position, object_time, record[4], photon, record[2])
                                           : StatusCode::SUCCESS;
        }
        if (fast_trace_only)
            return status;
        const T r_object = position[1], sin_theta_object = sinTheta(position[2]), cos_theta_object = cosTheta(position[2]),
                sin_phi_object = std::sin(position[3]), cos_phi_object = std::cos(position[3]);
        if (r_object <= 3.) {
            infra_context.logger_.warning("Object orbit radius = {:.6f}", r_object);
            if (r_object < 0)
                return StatusCode::FALL_INTO_BLACK_HOLE;
        }
        // [RYT-01] formulas/ryt/ray.py
        const T cos_observer_object = sin_theta_ * sin_theta_object * cos_phi_object + cos_theta_ * cos_theta_object;
        // [RYT-02] formulas/ryt/ray.py
        const T alpha_coefficient = sin_theta_object * sin_phi_object,
                beta_coefficient = cos_theta_object * sin_theta_ - sin_theta_object * cos_phi_object * cos_theta_,
                sin_observer_object = std::hypot(alpha_coefficient, beta_coefficient),
                theta_observer_object = std::atan2(sin_observer_object, cos_observer_object);
        ublas::bounded_vector<T, 2> alpha_beta_initial_value;
        if (cos_observer_object == -1.) {
            infra_context.logger_.warning("Object behind black hole, cos(theta) = {:.6f}\n", cos_observer_object);
            alpha_beta_initial_value(0) = 2. * std::sqrt(r_object);
        } else { // initial guessing
            T effective_radius;
            if (theta_observer_object < boost::math::constants::half_pi<T>())
                effective_radius =
                    r_object + power3(theta_observer_object / boost::math::constants::half_pi<T>()) / sin_observer_object;
            else
                // b-r*sin(theta)=(b-1.)*(2.*theta/pi-1.)+1.+(b-4.)/pi*sin(theta*2.)
                // b=(r*sin(theta)*boost::math::constants::pi<T>()+M_2_PI-2.*theta-8.*sin(theta)*cos(theta))/(2.*(boost::math::constants::pi<T>()t::math::constants::pi<T>()-theta-sin(theta)*cos(theta)))
                effective_radius =
                    1. / sin_observer_object +
                    (boost::math::constants::half_pi<T>() * r_object - 3. * cos_observer_object) /
                        (boost::math::constants::pi<T>() - theta_observer_object - sin_observer_object * cos_observer_object);
            alpha_beta_initial_value(0) = effective_radius * alpha_coefficient;
            alpha_beta_initial_value(1) = effective_radius * beta_coefficient;
        }
        TraceParameters trace_parameters(metric_, photon, photon_time, r_, r2_, theta_, sin_theta_, cos_theta_, position[1],
                                         sin_theta_object, cos_theta_object, sin_phi_object, cos_phi_object, t_final_,
                                         absolute_accuracy, relative_accuracy);
        GaussNewtonSolver<T, 2, TraceParameters<T>> alpha_beta_solver(TraceToPlane, trace_parameters);
        if (status = alpha_beta_solver.Set(alpha_beta_initial_value); status != StatusCode::SUCCESS)
            return status;
        if (status = alpha_beta_solver.Solve(r_object * boost::math::tools::root_epsilon<double>());
            status != StatusCode::SUCCESS)
            return status;
        // DNewtonRotationMultiFunctionSolver<T, 2, TraceParameters<T>> alpha_beta_rotation_solver(TraceToPlane,
        // trace_parameters); if (alpha_beta_rotation_solver.Set(alpha_beta_initial_value, theta_, sin_theta_, cos_theta_,
        // r_object, sin_theta_object, cos_theta_object, position[3], sin_phi_object, cos_phi_object, SPIN_CENTER_OFFSET *
        // metric_->a_ * sin_theta_object / r_object, 0., true) != StatusCode::SUCCESS ||
        // alpha_beta_rotation_solver.Solve(r_object * boost::math::tools::root_epsilon<double>()) != StatusCode::SUCCESS) {
        // 	DNewtonTranslationMultiFunctionSolver<T, 2, TraceParameters<T>> alpha_beta_translation_solver(TraceToPlane,
        // trace_parameters); 	if (alpha_beta_translation_solver.Set(alpha_beta_rotation_solver.Root(), theta_, sin_theta_,
        // cos_theta_, r_object, sin_theta_object, cos_theta_object, position[3], sin_phi_object, cos_phi_object, true) !=
        // StatusCode::SUCCESS || alpha_beta_translation_solver.Solve(r_object * boost::math::tools::root_epsilon<double>()) !=
        // StatusCode::SUCCESS) { 		DirectionMultiFunctionSolver<T, 2, TraceParameters<T>>
        // alpha_beta_direction_solver(TraceToPlane, trace_parameters); 		if (status =
        // alpha_beta_direction_solver.Set(alpha_beta_translation_solver.Root()); status != StatusCode::SUCCESS)
        // return status; 		if (status = alpha_beta_direction_solver.Solve(r_object *
        // boost::math::tools::root_epsilon<double>()); status != StatusCode::SUCCESS) 			return status;
        // alpha = alpha_beta_direction_solver.Root()(0); 		beta = alpha_beta_direction_solver.Root()(1); 	} else {
        // alpha = alpha_beta_translation_solver.Root()(0); 		beta = alpha_beta_translation_solver.Root()(1);
        // 	}
        // } else {
        // 	alpha = alpha_beta_rotation_solver.Root()(0);
        // 	beta = alpha_beta_rotation_solver.Root()(1);
        // }
        alpha = alpha_beta_solver.Root()(0);
        beta = alpha_beta_solver.Root()(1);
        photonInformation(position, object_time, record, photon, photon_time, alpha, beta);
        return calculate_magnification ? Magnification(position, object_time, record[4], photon, record[2])
                                       : StatusCode::SUCCESS;
    }

    static int TraceToPlane(const boost::numeric::ublas::bounded_vector<T, 2> &alpha_beta,
                            boost::numeric::ublas::bounded_vector<T, 2> &delta_apparent_alpha_beta,
                            TraceParameters<T> &params) {
        // [RYT-03] formulas/ryt/ray.py
        T alpha = alpha_beta(0), beta = alpha_beta(1);
        if (!std::isfinite(alpha) || !std::isfinite(beta))
            return StatusCode::NUMERIC_ERROR;
        T photon_time;
        if (int status = params.metric_.initializePhoton(params.photon_, photon_time, alpha, beta, params.r_, params.r2_,
                                                         params.theta_obs_, params.sin_theta_obs_);
            status != StatusCode::SUCCESS)
            return status;
        params.metric_.lagrangianToHamiltonian(params.photon_);
        T h = -0.01 * params.r_, last_h;
#ifdef SBODY_RELEASE
        auto wall_clock_start_time = std::chrono::steady_clock::now();
#endif
        auto integration_stepper = boost::numeric::odeint::make_dense_output(
            params.absolute_accuracy_, params.relative_accuracy_,
            boost::numeric::odeint::runge_kutta_dopri5<boost::numeric::ublas::bounded_vector<T, 8>>());
        integration_stepper.initialize(params.photon_, photon_time, h);
        T photon_last_theta = params.photon_[2], photon_hit_time, photon_hit_time_lower_bound,
          photon_hit_time_upper_bound = photon_time;
        while (photon_hit_time_upper_bound > params.t_final_) {
            try {
                integration_stepper.do_step(params.integration_system_);
            } catch (const std::exception &e) {
                return StatusCode::FAILURE;
            }
            const boost::numeric::ublas::bounded_vector<T, 8> &photon_state = integration_stepper.current_state();
            if (const T cos_observer_object_photon =
                    (params.x_obj_ - photon_state[1] * sinTheta(photon_state[2]) * std::cos(photon_state[3])) *
                        params.sin_theta_obs_ +
                    (params.z_obj_ - photon_state[1] * cosTheta(photon_state[2])) * params.cos_theta_obs_;
                cos_observer_object_photon >= 0.) {
                // photon goes through the plane of the object
                photon_hit_time_lower_bound = integration_stepper.current_time();
                boost::numeric::ublas::bounded_vector<T, 8> photon_hit_state;
                while (photon_hit_time_upper_bound - photon_hit_time_lower_bound > boost::math::tools::root_epsilon<T>()) {
                    photon_hit_time = 0.5 * (photon_hit_time_upper_bound + photon_hit_time_lower_bound);
                    integration_stepper.calc_state(photon_hit_time, photon_hit_state);
                    if ((params.x_obj_ - photon_hit_state[1] * sinTheta(photon_hit_state[2]) * std::cos(photon_hit_state[3])) *
                                params.sin_theta_obs_ +
                            (params.z_obj_ - photon_hit_state[1] * cosTheta(photon_hit_state[2])) * params.cos_theta_obs_ >
                        0.)
                        photon_hit_time_lower_bound = photon_hit_time;
                    else
                        photon_hit_time_upper_bound = photon_hit_time;
                }
                // photon in the same plane with the object
                params.photon_time_ = photon_hit_time - params.r_;
                delta_apparent_alpha_beta(0) =
                    params.y_obj_ - photon_hit_state[1] * sinTheta(photon_hit_state[2]) * std::sin(photon_hit_state[3]);
                delta_apparent_alpha_beta(1) =
                    (params.z_obj_ - photon_hit_state[1] * cosTheta(photon_hit_state[2])) * params.sin_theta_obs_ -
                    (params.x_obj_ - photon_hit_state[1] * sinTheta(photon_hit_state[2]) * std::cos(photon_hit_state[3])) *
                        params.cos_theta_obs_;
                return StatusCode::SUCCESS;
            }
            // photon fall into the BH
            if (std::abs(photon_state[5]) * boost::math::tools::root_epsilon<T>() > 1.) {
                delta_apparent_alpha_beta(0) = 1.1;
                return StatusCode::SCALING_REQUIRED;
            }
#ifdef SBODY_RELEASE
            if (std::chrono::steady_clock::now() - wall_clock_start_time > std::chrono::milliseconds(100)) {
                delta_apparent_alpha_beta(0) = 1.1;
                return StatusCode::COMPUTATION_TIMEOUT;
            }
#endif
            // photon has not reached the plane of the object
            photon_hit_time_upper_bound = integration_stepper.current_time();
            denseStepperMapTheta(integration_stepper, photon_last_theta);
        }
        // photon failed to hit the plane, the impact params need to be larger
        delta_apparent_alpha_beta(0) = 1.01;
        return StatusCode::SCALING_REQUIRED;
    }

    // [OBS-01] formulas/obs/obs.py
    int photonInformation(const boost::numeric::ublas::bounded_vector<T, 8> &position, TimeSystem object_time,
                          boost::numeric::ublas::bounded_vector<T, 5> &record,
                          const boost::numeric::ublas::bounded_vector<T, 8> &photon, T photon_time, T alpha, T beta) const {
        record(0) = alpha * cos_iota_ - beta * sin_iota_;                                // alpha
        record(1) = beta * cos_iota_ + alpha * sin_iota_;                                // beta
        record(2) = metric_.dopplerFactor(position, photon, object_time, TimeSystem::T); // Doppler factor g=1/(1+z)
        record(3) = photon_time;                                                         // look back time
        return StatusCode::SUCCESS;
    }

    /**
     * @brief Calculate the magnification of the object.
     *
     * @param position 8 dimensional vector
     * @param object_time time system of the object
     * @param magnification position to save the magnification
     * @param photon 8 dimensional vector of the photon traced to the object
     * @param doppler_factor Doppler factor \f$g=1/(1+z)\f$ of the photon
     * @return int
     */
    // [OBS-02] formulas/obs/obs.py
    int Magnification(const boost::numeric::ublas::bounded_vector<T, 8> &position, TimeSystem object_time, T &magnification,
                      const boost::numeric::ublas::bounded_vector<T, 8> &photon, T doppler_factor) const {
        namespace ublas = boost::numeric::ublas;
        boost::numeric::odeint::controlled_runge_kutta<
            boost::numeric::odeint::runge_kutta_fehlberg78<ublas::bounded_vector<T, 8>>>
            integration_stepper;
        ublas::bounded_matrix<T, 4, 4> gmunu;         // object local metric tensor
        ublas::bounded_matrix<T, 4, 4> coordinate;    // object local inertial coordinate frame
        ublas::bounded_vector<T, 4> photon_transform; // photon in TimeSystem TAU
        ublas::bounded_vector<T, 4> photon_in_object_frame_cartesian, photon_in_object_frame_spherical;
        metric_.metricTensor(position, gmunu);
        metric_.localInertialFrame(position, object_time, coordinate);
        ublas::bounded_matrix<T, 4, 4> coordinate_gmunu =
            ublas::prod(coordinate, gmunu); // object local inertial coordinate frame
        std::copy(photon.begin() + 5, photon.end(), photon_transform.begin() + 1);
        photon_transform(0) = 1.;
        photon_in_object_frame_cartesian = ublas::prod(coordinate_gmunu, photon_transform);
        photon_in_object_frame_cartesian /= T{-photon_in_object_frame_cartesian(0)};
        cartesianToSpherical(photon_in_object_frame_cartesian, photon_in_object_frame_spherical);
#ifndef SBODY_RELEASE
        ublas::bounded_matrix<T, 4, 4> coordinate_static;       // object local static frame (only dt/d\tau != 0)
        ublas::bounded_matrix<T, 4, 4> coordinate_static_gmunu; // object local static frame measured by observer
        ublas::bounded_vector<T, 4> photon_in_static_frame_cartesian;
        coordinate_static(0, 0) = std::sqrt(-1. / gmunu(0, 0));
        coordinate_static(1, 1) = std::sqrt(1. / gmunu(1, 1));
        coordinate_static(2, 2) = std::sqrt(1. / gmunu(2, 2));
        coordinate_static(3, 3) = std::sqrt(1. / gmunu(3, 3));
        coordinate_static_gmunu = ublas::prod(coordinate_static, gmunu);
        photon_in_static_frame_cartesian = ublas::prod(coordinate_static_gmunu, photon_transform);
        // gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu_gsl, coordinate_static, 0., coordinate_static_gmunu);
        // gsl_blas_dgemv(CblasNoTrans, 1., coordinate_static_gmunu, photon_transform, 0., photon_in_static_frame_cartesian);
        // local_redshift should equal to sqrt(EPSILON_POLYGON_AREA / cone_local_solid_angle),
        // the main error comes from the calculation of the cone_local_solid_angle, due to the limited SAMPLE_NUMBER.
        const T local_redshift = photon_in_object_frame_cartesian(0) / photon_in_static_frame_cartesian(0);
        photon_in_static_frame_cartesian /= T{photon_in_static_frame_cartesian(0)};
#endif
        ublas::bounded_vector<T, 8> forward_photon = photon;
        ublas::bounded_vector<T, 3> cone_record[SAMPLE_NUMBER], local_cone_record[SAMPLE_NUMBER], center_photon_velocity;
        metric_.normalizeNullGeodesic(forward_photon);
        metric_.lagrangianToHamiltonian(forward_photon);
        try {
            boost::numeric::odeint::integrate_adaptive(integration_stepper, integration_system_, forward_photon, 0.0, 1000.,
                                                       1.0);
        } catch (const std::exception &e) {
            return StatusCode::FAILURE;
        }
        metric_.hamiltonianToLagrangian(forward_photon);
        sphericalToCartesian(forward_photon);
        std::copy(forward_photon.begin() + 5, forward_photon.end(), center_photon_velocity.begin());
        center_photon_velocity /= ublas::norm_2(center_photon_velocity);
        ublas::permutation_matrix<std::size_t> permutation(4); // permutation used in the LU decomposition
        ublas::lu_factorize(coordinate_gmunu, permutation);
        auto forward_photon_velocity = ublas::vector_range(forward_photon, ublas::range(4, 8));
        for (int i = 0; i < SAMPLE_NUMBER; ++i) {
            const T angle = i * ANGLE_INTERVAL;
            std::copy(photon.begin(), photon.begin() + 4, forward_photon.begin());
            forward_photon[4] = -1.;
            forward_photon[5] = std::cos(angle) * T{SIN_EPSILON};
            forward_photon[6] = std::sin(angle) * T{SIN_EPSILON};
            forward_photon[7] = COS_EPSILON;
            rotateAroundAxis(forward_photon, true, Axis::Y, photon_in_object_frame_spherical[2]);
            rotateAroundAxis(forward_photon, true, Axis::Z, photon_in_object_frame_spherical[3]);
            ublas::lu_substitute(coordinate_gmunu, permutation, forward_photon_velocity);
            metric_.normalizeNullGeodesic(forward_photon);
#ifndef SBODY_RELEASE
            // local_cone_record[i][0] = metric_->DotProduct(photon, forward_photon + 4, coordinate_static->data + 4, 4);
            local_cone_record[i][0] = ublas::inner_prod(forward_photon_velocity, ublas::matrix_row(coordinate_static_gmunu, 1));
            // local_cone_record[i][1] = metric_->DotProduct(photon, forward_photon + 4, coordinate_static->data + 8, 4);
            local_cone_record[i][1] = ublas::inner_prod(forward_photon_velocity, ublas::matrix_row(coordinate_static_gmunu, 2));
            // local_cone_record[i][2] = metric_->DotProduct(photon, forward_photon + 4, coordinate_static->data + 12, 4);
            local_cone_record[i][2] = ublas::inner_prod(forward_photon_velocity, ublas::matrix_row(coordinate_static_gmunu, 3));
            local_cone_record[i] /= ublas::inner_prod(forward_photon_velocity, ublas::matrix_row(coordinate_static_gmunu, 0));
#endif
            metric_.lagrangianToHamiltonian(forward_photon);
            try {
                boost::numeric::odeint::integrate_adaptive(integration_stepper, integration_system_, forward_photon, 0.0, 1000.,
                                                           1.0);
            } catch (const std::exception &e) {
                return StatusCode::FAILURE;
            }
            metric_.hamiltonianToLagrangian(forward_photon);
            sphericalToCartesian(forward_photon);
            std::copy(forward_photon.begin() + 5, forward_photon.end(), cone_record[i].begin());
            cone_record[i] /= ublas::norm_2(cone_record[i]);
        }
        T cone_solid_angle = triangleArea(center_photon_velocity, cone_record[0], cone_record[SAMPLE_NUMBER - 1]);
#ifndef SBODY_RELEASE
        boost::numeric::ublas::bounded_vector<T, 3> photon_in_static_frame_cartesian_position;
        photon_in_static_frame_cartesian_position = ublas::vector_range(photon_in_static_frame_cartesian, ublas::range(1, 4));
        T cone_local_solid_angle =
            triangleArea(photon_in_static_frame_cartesian_position, local_cone_record[0], local_cone_record[SAMPLE_NUMBER - 1]);
#endif
        for (int i = 1; i < SAMPLE_NUMBER; ++i) {
            cone_solid_angle += triangleArea(center_photon_velocity, cone_record[i], cone_record[i - 1]);
#ifndef SBODY_RELEASE
            cone_local_solid_angle +=
                triangleArea(photon_in_static_frame_cartesian_position, local_cone_record[i], local_cone_record[i - 1]);
#endif
        }
        // The documented magnification is Omega_emit * g / cone_solid_angle.
        magnification = EPSILON_POLYGON_AREA * doppler_factor / cone_solid_angle;
        return StatusCode::SUCCESS;
    }

    /**
     * @brief Calculate the radius of the black hole shadow at different direction, saved in the file.
     *
     * @param file_name file name.
     * @return status
     */
    int Shadow(std::string file_name, const InfraContext &infra_context) const {
        boost::numeric::odeint::controlled_runge_kutta<
            boost::numeric::odeint::runge_kutta_fehlberg78<boost::numeric::ublas::bounded_vector<T, 8>>>
            integration_stepper;
        File record(FileFormat::NUMPY, file_name, {2});
        T h, rin = 2., rout = 10., rmid = 6.;
        boost::numeric::ublas::bounded_vector<T, 8> photon;
        T photon_time;
        int progress_bar_index = -1;
        if (infra_context.progress_bar_ != nullptr) {
            indicators::BlockProgressBar bar{
                indicators::option::ShowElapsedTime{true}, indicators::option::ShowRemainingTime{true},
                indicators::option::ForegroundColor{indicators::Color(4)},
                indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}};
            bar.set_option(indicators::option::MaxProgress(SAMPLE_NUMBER));
            bar.set_option(indicators::option::PrefixText("? Shadow"));
            progress_bar_index = infra_context.progress_bar_->push_back(bar);
        }
        for (int i = 0; i < SAMPLE_NUMBER; ++i) {
            const T angle = i * ANGLE_INTERVAL, sin_angle = std::sin(angle), cos_angle = std::cos(angle);
            while (rout - rin > boost::math::tools::root_epsilon<T>() * (rin + rout)) {
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
            record.save({rmid * (cos_angle * cos_iota_ - sin_angle * sin_iota_),
                         rmid * (sin_angle * cos_iota_ + cos_angle * sin_iota_)});
            if (progress_bar_index >= 0)
                infra_context.progress_bar_->tick(progress_bar_index);
        }
        if (progress_bar_index >= 0)
            infra_context.progress_bar_->setComplete(progress_bar_index, "! Shadow");
        return StatusCode::SUCCESS;
    }

    // int OmegaTest(std::optional<ProgressBar> &bars = std::nullopt) {
    // 	boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_fehlberg78<boost::numeric::ublas::bounded_vector<T,
    // 8>>> integration_stepper; 	T position_[8] = {0., 3., M_PI_4, 0., 0., 0., 0., 0.};
    // 	metric_->NormalizeTimelikeGeodesic(position_);
    // 	gsl_matrix *coordinate = gsl_matrix_alloc(4, 4), *gmunu_gsl = gsl_matrix_alloc(4, 4), *coordinate_gmunu =
    // gsl_matrix_alloc(4, 4); 	gsl_permutation *perm = gsl_permutation_alloc(4); 	metric_->MetricTensor(position_,
    // gmunu_gsl); 	gsl_matrix_set_zero(coordinate); 	gsl_matrix_set(coordinate, 0, 0, std::sqrt(-1. /
    // gmunu_gsl->data[0])); 	gsl_matrix_set(coordinate, 1, 1, std::sqrt(1. / gmunu_gsl->data[5]));
    // gsl_matrix_set(coordinate, 2, 2, std::sqrt(1. / gmunu_gsl->data[10])); 	gsl_matrix_set(coordinate, 3, 3, std::sqrt(1. /
    // gmunu_gsl->data[15])); 	gsl_matrix_set_zero(coordinate_gmunu); 	gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu_gsl,
    // coordinate, 0., coordinate_gmunu); 	int signum; 	gsl_linalg_LU_decomp(coordinate_gmunu, perm, &signum); 	T
    // rec[100][3], center_photon[3]; 	T photon_obs_time, area, h, time_limit = 0.; 	boost::numeric::ublas::bounded_vector<T,
    // 8> photon; 	gsl_vector_view ph_view = gsl_vector_view_array(photon.data() + 4, 4); 	NumPy
    // cone_record("Omega_record", {1}); 	bars.value()[0].set_option(indicators::option::MaxProgress(90)); 	for (T
    // angle = 90; angle > 0; angle -= 1) { 		const T sina = std::sin(angle / 180. * boost::math::constants::pi<T>()),
    // cosa = std::cos(angle / 180. * boost::math::constants::pi<T>()); 		copy(position_, position_ + 4,
    // photon.data()); photon(4) = 1.; photon(5) = sina; 		photon(6) = cosa; 		photon(7) = 0.;
    // gsl_linalg_LU_svx(coordinate_gmunu, perm, &ph_view.vector); 		metric_->NormalizeNullGeodesic(photon);
    // if (photon(5) < 0) { 			photon(5) = -photon(5); 			photon(6) = -photon(6);
    // photon(7) = -photon(7);
    // 		}
    // 		metric_->LagrangianToHamiltonian(photon);
    // 		photon_obs_time = 0;
    // 		h = 1.;
    // 		int status = 0;
    // 		boost::numeric::odeint::integrate_adaptive(integration_stepper, integration_system_, photon, time_limit,
    // -t_final_, 1.0); 		metric_->HamiltonianToLagrangian(photon); 		const T sin_theta
    // =sinTheta(photon[2]), cos_theta =cosTheta(photon[2]), sin_phi = std::sin(photon[3]), cos_phi = std::cos(photon[3]);
    // center_photon[0] = photon[5] * sin_theta * cos_phi + photon[1] * (cos_theta * cos_phi * photon[6] - sin_theta * sin_phi *
    // photon[7]); 		center_photon[1] = photon[5] * sin_theta * sin_phi + photon[1] * (cos_theta * sin_phi *
    // photon[6] + sin_theta * cos_phi * photon[7]); 		center_photon[2] = photon[5] * cos_theta - photon[1] * sin_theta
    // * photon[6]; 		const T vph_norm = Norm(center_photon); 		for (int j = 0; j < 3; ++j)
    // center_photon[j] /= vph_norm; 		for (int i = 0; i < 100; ++i) { 			const T angle_i = i *
    // ANGLE_INTERVAL, sinai = std::sin(angle_i), cosai = std::cos(angle_i); 			copy(position_, position_ + 4,
    // photon.begin()); 			photon[4] = 1.; 			photon[5] = sina - SIN_EPSILON * cosai *
    // cosa; 			photon[6] = cosa + SIN_EPSILON * cosai * sina; 			photon[7] = SIN_EPSILON * sinai;
    // gsl_linalg_LU_svx(coordinate_gmunu, perm, &ph_view.vector); 			metric_->NormalizeNullGeodesic(photon);
    // if (photon[5] < 0) { 				photon[5] = -photon[5]; 				photon[6] =
    // -photon[6]; 				photon[7] = -photon[7];
    // 			}
    // 			metric_->LagrangianToHamiltonian(photon);
    // 			photon_obs_time = 0;
    // 			h = 1.;
    // 			status = 0;
    // 			boost::numeric::odeint::integrate_adaptive(integration_stepper, integration_system_, photon,
    // photon_obs_time, time_limit, 1.0); 			metric_->HamiltonianToLagrangian(photon);
    // const T sin_theta =sinTheta(photon[2]), cos_theta =cosTheta(photon[2]), sin_phi = std::sin(photon[3]), cos_phi =
    // std::cos(photon[3]); 			rec[i][0] = photon[5] * sin_theta * cos_phi + photon[1] * (cos_theta * cos_phi *
    // photon[6] - sin_theta * sin_phi * photon[7]); 			rec[i][1] = photon[5] * sin_theta * sin_phi + photon[1]
    // * (cos_theta * sin_phi * photon[6] + sin_theta * cos_phi * photon[7]); 			rec[i][2] = photon[5] *
    // cos_theta - photon[1] * sin_theta * photon[6]; 			const T vph_norm = Norm(rec[i]);
    // for (int j = 0; j < 3; ++j) 				rec[i][j] /= vph_norm;
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
    // 		cone_record.Save({std::abs(area) / (boost::math::constants::two_pi<T>() * boost::math::tools::epsilon<T>())});
    // 		bars.value()[0].tick();
    // 	}
    // 	return StatusCode::SUCCESS;
    // }
};
} // namespace sbody
