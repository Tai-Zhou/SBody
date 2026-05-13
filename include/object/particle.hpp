/**
 * @file particle.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <boost/numeric/odeint.hpp>

#include "math/geometry.hpp"
#include "object/object.hpp"

namespace sbody {
/**
 * @brief
 *
 */
template <typename T>
class Particle : public Object<T> {
  protected:
    TimeSystem time_;
    DynamicalSystem coordinate_;
    /// 8 dimensional information
    boost::numeric::ublas::bounded_vector<T, 8> position_;
    /// if the position fixed
    const bool fixed_;
    /// integration
    std::function<void(const boost::numeric::ublas::bounded_vector<T, 8> &, boost::numeric::ublas::bounded_vector<T, 8> &, T)>
        integration_system_;
    boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_fehlberg78<
        boost::numeric::ublas::bounded_vector<T, 8>, T, boost::numeric::ublas::bounded_vector<T, 8>, T>>
        integration_stepper_;

  public:
    /**
     * @brief Construct a new star object
     *
     * @param metric
     * @param radius radius
     * @param fixed whether the position of the star is fixed
     */
    Particle(const Metric<T> &metric, TimeSystem time, DynamicalSystem coordinate, T absolute_accuracy, T relative_accuracy,
             bool fixed = false)
        : Object<T>(metric), time_(time), coordinate_(coordinate), fixed_(fixed),
          integration_stepper_(
              boost::numeric::odeint::make_controlled<boost::numeric::odeint::runge_kutta_fehlberg78<
                  boost::numeric::ublas::bounded_vector<T, 8>, T, boost::numeric::ublas::bounded_vector<T, 8>, T>>(
                  absolute_accuracy, relative_accuracy)) {}
    T position(std::size_t idx) { return position_(idx); }
    void position(T *position) { std::copy(position_.begin(), position_.end(), position); }
    void position(boost::numeric::ublas::bounded_vector<T, 8> &position) { position = position_; }
    int initializeKeplerian(T a, T e, T inclination, T periapsis, T ascending_node, T true_anomaly, T observer_inclination = 0.,
                            T observer_rotation = 0.) {
        if (int status = this->metric_.getIntegrationSystem(integration_system_, time_, coordinate_, MotionMode::GEODESIC);
            status != StatusCode::SUCCESS)
            return status;
        position_[0] = 0.;
        periapsis += true_anomaly;
        T r = a * (1. - e * e) / (1. + e * std::cos(true_anomaly));
        const T tp1 = -r * std::cos(periapsis), tp2 = -r * std::sin(periapsis) * std::cos(inclination);
        const T xp1 = tp1 * std::cos(ascending_node) - tp2 * std::sin(ascending_node),
                xp2 = tp2 * std::cos(ascending_node) + tp1 * std::sin(ascending_node),
                xp3 = -r * std::sin(periapsis) * std::sin(inclination);
        position_[1] =
            (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::cos(observer_inclination) +
            xp3 * std::sin(observer_inclination);
        position_[2] = xp2 * std::cos(observer_rotation) - xp1 * std::sin(observer_rotation);
        position_[3] = xp3 * std::cos(observer_inclination) -
                       (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::sin(observer_inclination);
        if (fixed_) {
            position_[5] = 0.;
            position_[6] = 0.;
            position_[7] = 0.;
        } else {
            const T vphi = std::sqrt((1 - e * e) * a) / r,
                    vr = std::copysign(squareRoot(2. / r - 1. / a - vphi * vphi), phiDifference(true_anomaly));
            const T tp5 = vphi * std::sin(periapsis) - vr * std::cos(periapsis),
                    tp6 = -(vphi * std::cos(periapsis) + vr * std::sin(periapsis)) * std::cos(inclination);
            const T xp5 = tp5 * std::cos(ascending_node) - tp6 * std::sin(ascending_node),
                    xp6 = tp5 * std::sin(ascending_node) + tp6 * std::cos(ascending_node),
                    xp7 = -(vphi * std::cos(periapsis) + vr * std::sin(periapsis)) * std::sin(inclination);
            position_[5] =
                (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::cos(observer_inclination) +
                xp7 * std::sin(observer_inclination);
            position_[6] = xp6 * std::cos(observer_rotation) - xp5 * std::sin(observer_rotation);
            position_[7] =
                xp7 * std::cos(observer_inclination) -
                (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::sin(observer_inclination);
        }
        cartesianToSpherical(position_);
        return this->metric_.normalizeTimelikeGeodesic(position_);
    }
    int initializeGeodesic(T orbital_radius, T inclination, T periapsis, T ascending_node, T v_r, T v_phi,
                           T observer_inclination = 0., T observer_rotation = 0.) {
        if (int status = this->metric_.getIntegrationSystem(integration_system_, time_, coordinate_, MotionMode::GEODESIC);
            status != StatusCode::SUCCESS)
            return status;
        position_[0] = 0.;
        const T tp1 = -orbital_radius * std::cos(periapsis),
                tp2 = -orbital_radius * std::sin(periapsis) * std::cos(inclination);
        const T xp1 = tp1 * std::cos(ascending_node) - tp2 * std::sin(ascending_node),
                xp2 = tp2 * std::cos(ascending_node) + tp1 * std::sin(ascending_node),
                xp3 = -orbital_radius * std::sin(periapsis) * std::sin(inclination);
        position_[1] =
            (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::cos(observer_inclination) +
            xp3 * std::sin(observer_inclination);
        position_[2] = xp2 * std::cos(observer_rotation) - xp1 * std::sin(observer_rotation);
        position_[3] = xp3 * std::cos(observer_inclination) -
                       (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::sin(observer_inclination);
        if (fixed_) {
            position_[5] = 0.;
            position_[6] = 0.;
            position_[7] = 0.;
        } else {
            const T tp5 = v_phi * std::sin(periapsis) - v_r * std::cos(periapsis),
                    tp6 = -(v_phi * std::cos(periapsis) + v_r * std::sin(periapsis)) * std::cos(inclination);
            const T xp5 = tp5 * std::cos(ascending_node) - tp6 * std::sin(ascending_node),
                    xp6 = tp5 * std::sin(ascending_node) + tp6 * std::cos(ascending_node),
                    xp7 = -(v_phi * std::cos(periapsis) + v_r * std::sin(periapsis)) * std::sin(inclination);
            position_[5] =
                (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::cos(observer_inclination) +
                xp7 * std::sin(observer_inclination);
            position_[6] = xp6 * std::cos(observer_rotation) - xp5 * std::sin(observer_rotation);
            position_[7] =
                xp7 * std::cos(observer_inclination) -
                (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::sin(observer_inclination);
        }
        cartesianToSpherical(position_);
        return this->metric_.normalizeTimelikeGeodesic(position_);
    }
    int initializeSchwarzschildKeplerianPericenter(T a, T e, T inclination, T periapsis, T ascending_node,
                                                   T observer_inclination, T observer_rotation) {
        if (int status = this->metric_.getIntegrationSystem(integration_system_, time_, coordinate_, MotionMode::GEODESIC);
            status != StatusCode::SUCCESS)
            return status;
        position_[0] = 0.;
        T r = a * (1. - e); // pericenter
        const T tp1 = -r * std::cos(periapsis), tp2 = -r * std::sin(periapsis) * std::cos(inclination);
        const T xp1 = tp1 * std::cos(ascending_node) - tp2 * std::sin(ascending_node),
                xp2 = tp2 * std::cos(ascending_node) + tp1 * std::sin(ascending_node),
                xp3 = -r * std::sin(periapsis) * std::sin(inclination);
        position_[1] =
            (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::cos(observer_inclination) +
            xp3 * std::sin(observer_inclination);
        position_[2] = xp2 * std::cos(observer_rotation) - xp1 * std::sin(observer_rotation);
        position_[3] = xp3 * std::cos(observer_inclination) -
                       (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::sin(observer_inclination);
        if (fixed_) {
            position_[5] = 0.;
            position_[6] = 0.;
            position_[7] = 0.;
        } else {
            const T E = (a - 2.) / std::sqrt(a * (a - 3.));
            const T vphi = std::sqrt(E * E * r * (r - 2.) - power2(r - 2.)) / (r * E);
            const T tp5 = vphi * std::sin(periapsis), tp6 = -vphi * std::cos(periapsis) * std::cos(inclination);
            const T xp5 = tp5 * std::cos(ascending_node) - tp6 * std::sin(ascending_node),
                    xp6 = tp5 * std::sin(ascending_node) + tp6 * std::cos(ascending_node),
                    xp7 = -vphi * std::cos(periapsis) * std::sin(inclination);
            position_[5] =
                (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::cos(observer_inclination) +
                xp7 * std::sin(observer_inclination);
            position_[6] = xp6 * std::cos(observer_rotation) - xp5 * std::sin(observer_rotation);
            position_[7] =
                xp7 * std::cos(observer_inclination) -
                (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::sin(observer_inclination);
        }
        cartesianToSpherical(position_);
        return this->metric_.normalizeTimelikeGeodesic(position_);
    }
    int initializeSchwarzschildKeplerianApocenter(T a, T e, T inclination, T periapsis, T ascending_node,
                                                  T observer_inclination, T observer_rotation) {
        if (int status = this->metric_.getIntegrationSystem(integration_system_, time_, coordinate_, MotionMode::GEODESIC);
            status != StatusCode::SUCCESS)
            return status;
        position_[0] = 0.;
        T r = a * (1. + e); // apocenter
        const T tp1 = r * std::cos(periapsis), tp2 = r * std::sin(periapsis) * std::cos(inclination);
        const T xp1 = tp1 * std::cos(ascending_node) - tp2 * std::sin(ascending_node),
                xp2 = tp2 * std::cos(ascending_node) + tp1 * std::sin(ascending_node),
                xp3 = r * std::sin(periapsis) * std::sin(inclination);
        position_[1] =
            (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::cos(observer_inclination) +
            xp3 * std::sin(observer_inclination);
        position_[2] = xp2 * std::cos(observer_rotation) - xp1 * std::sin(observer_rotation);
        position_[3] = xp3 * std::cos(observer_inclination) -
                       (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::sin(observer_inclination);
        if (fixed_) {
            position_[5] = 0.;
            position_[6] = 0.;
            position_[7] = 0.;
        } else {
            const T E = (a - 2.) / std::sqrt(a * (a - 3.));
            const T vphi = std::sqrt(E * E * r * (r - 2.) - power2(r - 2.)) / (r * E);
            const T tp5 = -vphi * std::sin(periapsis), tp6 = vphi * std::cos(periapsis) * std::cos(inclination);
            const T xp5 = tp5 * std::cos(ascending_node) - tp6 * std::sin(ascending_node),
                    xp6 = tp5 * std::sin(ascending_node) + tp6 * std::cos(ascending_node),
                    xp7 = vphi * std::cos(periapsis) * std::sin(inclination);
            position_[5] =
                (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::cos(observer_inclination) +
                xp7 * std::sin(observer_inclination);
            position_[6] = xp6 * std::cos(observer_rotation) - xp5 * std::sin(observer_rotation);
            position_[7] =
                xp7 * std::cos(observer_inclination) -
                (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::sin(observer_inclination);
        }
        cartesianToSpherical(position_);
        return this->metric_.normalizeTimelikeGeodesic(position_);
    }
    int initializeKeplerianHarmonic(T a, T e, T inclination, T periapsis, T ascending_node, T true_anomaly,
                                    T observer_inclination, T observer_rotation) {
        initializeKeplerian(a, e, inclination, periapsis, ascending_node, true_anomaly, observer_inclination,
                            observer_rotation);
        position_[1] += 1.;
        return this->metric_.normalizeTimelikeGeodesic(position_);
    }
    int initializeCircular(T r, T phi, T v_phi_ratio) {
        if (int status = this->metric_.getIntegrationSystem(integration_system_, time_, coordinate_, MotionMode::CIRCULAR);
            status != StatusCode::SUCCESS)
            return status;
        position_[0] = 0.;
        position_[1] = r;
        position_[2] = boost::math::constants::half_pi<T>();
        position_[3] = phi;
        if (fixed_) {
            position_[5] = 0.;
            position_[6] = 0.;
            position_[7] = 0.;
        } else {
            position_[5] = 0.;
            position_[6] = 0.;
            position_[7] = v_phi_ratio / (r * std::sqrt(r));
        }
        return this->metric_.normalizeTimelikeGeodesic(position_);
    }
    int initializeConical(T r, T theta, T phi, T v_r, T v_phi, bool fixed_energy) {
        if (int status = this->metric_.getIntegrationSystem(
                integration_system_, time_, coordinate_, fixed_energy ? MotionMode::CONICAL_FIXED_ENERGY : MotionMode::CONICAL);
            status != StatusCode::SUCCESS)
            return status;
        position_[0] = 0.;
        position_[1] = r;
        position_[2] = theta;
        position_[3] = phi;
        if (fixed_) {
            position_[5] = 0.;
            position_[6] = 0.;
            position_[7] = 0.;
        } else {
            position_[5] = v_r;
            position_[6] = 0.;
            position_[7] = v_phi;
        }
        return this->metric_.normalizeTimelikeGeodesic(position_);
    }
    int initializeCylindrical(T r, T phi, T z, T v_z, T v_phi, bool fixed_energy) {
        if (int status = this->metric_.getIntegrationSystem(integration_system_, time_, coordinate_,
                                                            fixed_energy ? MotionMode::CYLINDRICAL_FIXED_ENERGY
                                                                         : MotionMode::CYLINDRICAL);
            status != StatusCode::SUCCESS)
            return status;
        position_[0] = 0.;
        position_[1] = std::hypot(r, z);
        position_[2] = std::atan2(r, z);
        position_[3] = phi;
        if (fixed_) {
            position_[5] = 0.;
            position_[6] = 0.;
            position_[7] = 0.;
        } else {
            position_[5] = v_z * std::cos(position_[2]);
            position_[6] = -v_z * std::sin(position_[2]) / position_[1];
            position_[7] = v_phi;
        }
        return this->metric_.normalizeTimelikeGeodesic(position_);
    }
    int integratorApply(T &t, T t1, T &h) {
        if (fixed_) {
            t = t1;
            return StatusCode::SUCCESS;
        }
        if (h == 0.)
            h = std::copysign(std::max(std::sqrt(std::abs(t1 - t)), T{1.}), t1 - t);
        T last_theta = position_[2];
        if (h > 0) {
            while (t1 - t > 0.) {
                if (h > t1 - t)
                    h = t1 - t;
                try {
                    integration_stepper_.try_step(integration_system_, position_, t, h);
                } catch (const std::exception &e) {
                    return StatusCode::FAILURE;
                }
                mapTheta(position_, last_theta);
                if (h * 10000 < t1 - t)
                    return StatusCode::STEP_SIZE_TOO_SMALL;
            }
        } else {
            while (t1 - t < 0.) {
                if (h < t1 - t)
                    h = t1 - t;
                try {
                    integration_stepper_.try_step(integration_system_, position_, t, h);
                } catch (const std::exception &e) {
                    return StatusCode::FAILURE;
                }
                mapTheta(position_, last_theta);
                if (h * 10000 > t1 - t)
                    return StatusCode::STEP_SIZE_TOO_SMALL;
            }
        }
        return StatusCode::SUCCESS;
    }
    int integratorApplyStep(T &t, const T h) {
        if (fixed_) {
            t += h;
            return StatusCode::SUCCESS;
        }
        T start_time = t;
        try {
            integration_stepper_.try_step(integration_system_, position_, t, h);
        } catch (const std::exception &e) {
            return StatusCode::FAILURE;
        }
        if (t == start_time)
            return StatusCode::FAILURE;
        return StatusCode::SUCCESS;
    }
    int normalize() { return this->metric_.normalizeTimelikeGeodesic(position_); }

    /**
     * @brief Check if the star hit by the photon
     *
     * @param current
     * @param last
     * @return bool
     */
    bool hit(const boost::numeric::ublas::bounded_vector<T, 8> &current,
             const boost::numeric::ublas::bounded_vector<T, 8> &last) const override {
        return false;
    }

    /**
     * @brief return frequency
     *
     * @param photon
     * @param photon_time
     * @return T
     */
    T redshift(const boost::numeric::ublas::bounded_vector<T, 8> &photon, TimeSystem photon_time) const override {
        return this->metric_.redshift(position_, photon, time_, photon_time);
    }

    int metricTensor(boost::numeric::ublas::bounded_matrix<T, 4, 4> &metric) const {
        return this->metric_.metricTensor(position_, metric);
    }
    template <typename VecExp1, typename VecExp2>
    T dotProduct(const VecExp1 &x, const VecExp2 &y, bool include_time) const {
        return this->metric_.dotProduct(position_, x, y, include_time);
    }
    int localInertialFrame(boost::numeric::ublas::bounded_matrix<T, 4, 4> &coordinate) const {
        return this->metric_.localInertialFrame(position_, time_, coordinate);
    }
    T energy() const { return this->metric_.energy(position_, time_, coordinate_); }
    T angularMomentum() const { return this->metric_.angularMomentum(position_, time_, coordinate_); }
    T carterConstant() const { return this->metric_.carterConstant(position_, 1., time_, coordinate_); }
};

template <typename T>
class Star : public Particle<T> {
  protected:
    /// radius
    const T radius_;
    /// radius square
    const T radius_square_;

  public:
    Star(const Metric<T> &metric, TimeSystem time, DynamicalSystem coordinate, T absolute_accuracy, T relative_accuracy,
         T radius = 0, bool fixed = false)
        : Particle<T>(metric, time, coordinate, absolute_accuracy, relative_accuracy, fixed), radius_(radius),
          radius_square_(radius * radius) {}
    bool hit(const boost::numeric::ublas::bounded_vector<T, 8> &current,
             const boost::numeric::ublas::bounded_vector<T, 8> &last) const override {
        T a2 = this->metric_.distanceSquare(this->position_, current, 3);
        if (a2 <= radius_square_)
            return true;
        T b2 = this->metric_.distanceSquare(this->position_, last, 3), c2 = this->metric_.distanceSquare(current, last, 3);
        if (a2 + c2 > b2 && b2 + c2 > a2 && 2 * a2 * c2 - power2(a2 - b2 + c2) <= 4 * c2 * radius_square_)
            return true; // if min distance between current and last < radius, return 1;
        return false;
    }
};

template <typename T>
class Hotspot : public Particle<T> {
  protected:
    const T spectral_index_;
    const T luminosity_;
    const T luminosity_mu_;
    const T luminosity_sigma_;

  public:
    Hotspot(const Metric<T> &metric, TimeSystem time, DynamicalSystem coordinate, T absolute_accuracy, T relative_accuracy,
            T spectral_index, T luminosity, T luminosity_mu, T luminosity_sigma, bool fixed)
        : Particle<T>(metric, time, coordinate, absolute_accuracy, relative_accuracy, fixed), spectral_index_(spectral_index),
          luminosity_(luminosity), luminosity_mu_(luminosity_mu), luminosity_sigma_(luminosity_sigma) {}
    T luminosity(T t) const {
        T intrinsic_time = this->time_ == TimeSystem::T ? this->position_[0] : t;
        return luminosity_ * exp(-0.5 * power2((intrinsic_time - luminosity_mu_) / luminosity_sigma_));
    }
    T spectralDensity(T t, T redshift) const {
        T intrinsic_time = this->time_ == TimeSystem::T ? this->position_[0] : t;
        return luminosity_ * pow(redshift, spectral_index_) *
               exp(-0.5 * power2((intrinsic_time - luminosity_mu_) / luminosity_sigma_));
    }
};
} // namespace sbody
