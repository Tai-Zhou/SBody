/**
 * @file schwarzschild.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-06
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
// #include <boost/math/special_functions/next.hpp>
#include <boost/math/tools/roots.hpp>

#include "infra/infra_context.hpp"
#include "math/constants.hpp"
#include "math/geometry.hpp"
#include "math/solver/polynomial.hpp"
#include "math/special_functions.hpp"
#include "metric/metric.hpp"

namespace sbody {
template <typename T>
class Schwarzschild : public Metric<T> {
  public:
    Schwarzschild() = default;
    std::string name() const override { return "Schwarzschild"; }
    // [MET-08] formulas/met/metric.py
    int angularMomentumCarterConstantToAlphaBeta(T l, T q2, T cos_theta, T sin_theta, T &alpha, T &abs_beta) const override {
        if (sin_theta == 0.)
            return StatusCode::FAILURE;
        alpha = -l / sin_theta;
        const T beta2 = q2 - power2(cos_theta * alpha);
        if (beta2 < 0.)
            return StatusCode::DOMAIN_ERROR;
        abs_beta = std::sqrt(beta2);
        return StatusCode::SUCCESS;
    }
    // [MET-03] formulas/met/metric.py
    int metricTensor(const boost::numeric::ublas::bounded_vector<T, 8> &position,
                     boost::numeric::ublas::bounded_matrix<T, 4, 4> &metric) const override {
        metric.clear();
        metric(0, 0) = -(1. - 2. / position[1]);
        metric(1, 1) = position[1] / (position[1] - 2.);
        metric(2, 2) = position[1] * position[1];
        metric(3, 3) = power2(position[1] * std::sin(position[2]));
        return position[1] == 2. ? StatusCode::DIVISION_BY_ZERO : StatusCode::SUCCESS;
    }
    // [MET-04] formulas/met/metric.py
    T dotProduct3D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x1, T x2, T x3, T y1, T y2,
                   T y3) const override {
        return position[1] * x1 * y1 / (position[1] - 2.) + position[1] * position[1] * x2 * y2 +
               power2(position[1] * std::sin(position[2])) * x3 * y3;
    }
    T dotProduct4D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x0, T x1, T x2, T x3, T y0, T y1, T y2,
                   T y3) const override {
        return -(1. - 2. / position[1]) * x0 * y0 + position[1] * x1 * y1 / (position[1] - 2.) +
               position[1] * position[1] * x2 * y2 + power2(position[1] * std::sin(position[2])) * x3 * y3;
    }
    // [MET-04] formulas/met/metric.py
    T distanceSquare(const boost::numeric::ublas::bounded_vector<T, 8> &x, const boost::numeric::ublas::bounded_vector<T, 8> &y,
                     bool include_time) const override {
        if (include_time)
            return -(1. - 2. / x[1]) * power2(x[0] - y[0]) + x[1] * power2(x[1] - y[1]) / (x[1] - 2.) +
                   power2(x[1] * (x[2] - y[2])) + power2(x[1] * std::sin(x[2]) * phiDifference(x[3] - y[3]));
        return x[1] * power2(x[1] - y[1]) / (x[1] - 2.) + power2(x[1] * (x[2] - y[2])) +
               power2(x[1] * std::sin(x[2]) * phiDifference(x[3] - y[3]));
    }
    // [MET-06] formulas/met/metric.py
    int lagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T dt_dtau = 1. / y[4], r2 = power2(y[1]);
        y[4] = (y[4] - 1. + 2. / y[1]) * dt_dtau; // 1 + p_t
        y[5] *= y[1] / (y[1] - 2.) * dt_dtau;
        y[6] *= r2 * dt_dtau;
        y[7] *= r2 * power2(std::sin(y[2])) * dt_dtau;
        return StatusCode::SUCCESS;
    }
    // [MET-06] formulas/met/metric.py
    int hamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T r_1 = 1. / y[1], g11_1 = 1. - 2. * r_1;
        y[4] = -g11_1 / (y[4] - 1.);
        y[5] *= g11_1 * y[4];
        y[6] *= r_1 * r_1 * y[4];
        y[7] *= power2(r_1 / std::sin(y[2])) * y[4];
        return StatusCode::SUCCESS;
    }
    int fastTrace(const T r_observer, const T theta_observer, const T sin_theta_observer, const T cos_theta_observer,
                  const T r_object, const T theta_object, const T phi_object, T &alpha, T &beta,
                  boost::numeric::ublas::bounded_vector<T, 8> &photon, T &photon_time, T absolute_accuracy, T relative_accuracy,
                  const InfraContext &infra_context) const override {
        // [RYT-04] formulas/ryt/ray.py
        const T sin_theta_object = sinTheta(theta_object), cos_theta_object = cosTheta(theta_object),
                sin_phi_object = std::sin(phi_object), cos_phi_object = std::cos(phi_object);
        const T sin_observer_object =
            std::hypot(sin_theta_observer * cos_theta_object * cos_phi_object - cos_theta_observer * sin_theta_object,
                       sin_theta_observer * sin_phi_object);
        const T cos_observer_object =
            sin_theta_observer * sin_theta_object * cos_phi_object + cos_theta_observer * cos_theta_object;
        const T delta_phi = std::atan2(sin_observer_object, cos_observer_object);
        const T u0 = 1. / r_observer, u1 = 1. / r_object, u12 = u1 * u1;
        const T g11_1 = 1. - 2. * u1;
        if (delta_phi == 0.) {
            alpha = 0.;
            beta = 0.;
            photon(0) = (r_object - r_observer + 4. * log(r_object * u0) + 4. * (u1 - u0)) * r_observer / (r_observer - 2.);
            photon(1) = r_object;
            photon(2) = theta_object;
            photon(3) = phi_object;
            photon(4) = r_object * (r_observer - 2.) / (r_observer * (r_object - 2.));
            photon(5) = (r_object - 2.) * u1;
            photon(6) = 0.;
            photon(7) = 0.;
            photon_time = r_object - r_observer + 2. * log(r_object * u0);
            return StatusCode::SUCCESS;
        }
        T impact_parameter_upper_limit = r_object / std::sqrt(g11_1), turning_phi;
        if (T x02[2];
            u1 * 3. < 1. && polySolveQuadratic(T{2.}, -g11_1, -g11_1 * u1, x02) == 2) { // remove (x-u1) as cubic_roots[1] = u1
            turning_phi = boost::math::constants::half_root_two<T>() *
                          ellipticIntegral(0, u0, u1, T{0.}, T{0.}, u1, T{-1.}, -x02[0], T{1.}, x02[1], T{-1.});
        } else { // there is no turning point on the trajectory
            impact_parameter_upper_limit = T{M_SQRT27} - boost::math::tools::epsilon<T>();
            turning_phi = boost::math::constants::pi<T>(); // make delta_phi <= turning_phi
        }
        std::pair<T, T> impact_root;
        std::array<T, 3> cubic_roots;
        try {
            boost::uintmax_t max_iteration = 64;
            if (delta_phi > turning_phi) {
                impact_root = boost::math::tools::toms748_solve(
                    [u0, u1, delta_phi, &cubic_roots](T l) -> T {
                        if (l == 0.)
                            return -delta_phi;
                        T half_l_2 = 0.5 / (l * l);
                        if (polySolveCubic(T{-0.5}, T{0.}, half_l_2, cubic_roots.begin()) == 3) {
                            if (cubic_roots[1] < u1) { // cubic_roots[1] < u1 < cubic_roots[1] + EPSILON
                                return boost::math::constants::half_root_two<T>() *
                                           ellipticIntegral(0, u0, u1, T{0.}, T{1.}, u1, T{-1.}, -cubic_roots[0], T{1.},
                                                            cubic_roots[2], T{-1.}) -
                                       delta_phi;
                            } else { // u1 -> turning point -> u1 -> u0
                                return boost::math::constants::root_two<T>() *
                                           ellipticIntegral(0., u1, cubic_roots[1], T{0.}, T{1.}, cubic_roots[1], T{-1.},
                                                            -cubic_roots[0], T{1.}, cubic_roots[2], T{-1.}) +
                                       boost::math::constants::half_root_two<T>() *
                                           ellipticIntegral(0, u0, u1, T{0.}, T{1.}, cubic_roots[1], T{-1.}, -cubic_roots[0],
                                                            T{1.}, cubic_roots[2], T{-1.}) -
                                       delta_phi;
                            }
                        }
                        // impact_parameter < std::sqrt(27)
                        cubic_roots[1] = std::nan("");
                        return boost::math::constants::half_root_two<T>() *
                                   ellipticIntegral2Complex(0, u0, u1, T{0.}, T{1.}, -half_l_2 / cubic_roots[0],
                                                            cubic_roots[0] - 0.5, T{1.}, -cubic_roots[0], T{1.}) -
                               delta_phi;
                    },
                    T{M_SQRT27} + boost::math::tools::root_epsilon<T>(), impact_parameter_upper_limit,
                    boost::math::tools::eps_tolerance<double>(), max_iteration);
            } else {
                impact_root = boost::math::tools::toms748_solve(
                    [u0, u1, delta_phi, &cubic_roots](T l) -> T {
                        if (l == 0.)
                            return -delta_phi;
                        T half_l_2 = 0.5 / (l * l);
                        if (polySolveCubic(T{-0.5}, T{0.}, half_l_2, cubic_roots.begin()) == 3) {
                            if (cubic_roots[1] < u1) { // cubic_roots[1] < u1 < cubic_roots[1] + EPSILON
                                return boost::math::constants::half_root_two<T>() *
                                           ellipticIntegral(0, u0, u1, T{0.}, T{1.}, u1, T{-1.}, -cubic_roots[0], T{1.},
                                                            cubic_roots[2], T{-1.}) -
                                       delta_phi;
                            } else { // u1 -> u0
                                return boost::math::constants::half_root_two<T>() *
                                           ellipticIntegral(0, u0, u1, T{0.}, T{1.}, cubic_roots[1], T{-1.}, -cubic_roots[0],
                                                            T{1.}, cubic_roots[2], T{-1.}) -
                                       delta_phi;
                            }
                        }
                        // impact_parameter < std::sqrt(27)
                        cubic_roots[1] = std::nan("");
                        return boost::math::constants::half_root_two<T>() *
                                   ellipticIntegral2Complex(0, u0, u1, T{0.}, T{1.}, -half_l_2 / cubic_roots[0],
                                                            cubic_roots[0] - 0.5, T{1.}, -cubic_roots[0], T{1.}) -
                               delta_phi;
                    },
                    T{0.}, impact_parameter_upper_limit, boost::math::tools::eps_tolerance<double>(), max_iteration);
            }
        } catch (const std::exception &e) {
            infra_context.logger_.warning("Failed to find impact root: {}, delta_phi {} turning_phi", e.what(),
                                          delta_phi < turning_phi ? "<" : ">");
            return StatusCode::FAILURE;
        }
        T impact_root_value = 0.5 * (impact_root.first + impact_root.second),
          impact_root_value2 = impact_root_value * impact_root_value;
        alpha = impact_root_value / sin_observer_object * sin_theta_object * sin_phi_object;
        beta = impact_root_value / sin_observer_object *
               (cos_theta_object * sin_theta_observer - sin_theta_object * cos_phi_object * cos_theta_observer);
        if (!std::isnan(cubic_roots[1])) {
            if (cubic_roots[1] < u1) { // cubic_roots[1] < u1 < cubic_roots[1] + EPSILON
                const T ellip_int_4 =
                    ellipticIntegral(-4, u0, u1, T{0.}, T{1.}, u1, T{-1.}, -cubic_roots[0], T{1.}, cubic_roots[2], T{-1.});
                photon(0) = -boost::math::constants::half_root_two<T>() * ellip_int_4 / (impact_root_value * (1. - 2. * u0));
                photon_time = -boost::math::constants::half_root_two<T>() *
                              (ellip_int_4 +
                               2. * ellipticIntegral(-2, u0, u1, T{0.}, T{1.}, u1, T{-1.}, -cubic_roots[0], T{1.},
                                                     cubic_roots[2], T{-1.}) +
                               4. * ellipticIntegral(-2, u0, u1, T{1.}, T{-2.}, u1, T{-1.}, -cubic_roots[0], T{1.},
                                                     cubic_roots[2], T{-1.})) /
                              impact_root_value;
            } else if (delta_phi > turning_phi) { // u1 -> turning point -> u1 -> u0
                const T ellip_int_4 = 2. * ellipticIntegral(-4, u1, cubic_roots[1], T{0.}, T{1.}, cubic_roots[1], T{-1.},
                                                            -cubic_roots[0], T{1.}, cubic_roots[2], T{-1.}) +
                                      ellipticIntegral(-4, u0, u1, T{0.}, T{1.}, cubic_roots[1], T{-1.}, -cubic_roots[0], T{1.},
                                                       cubic_roots[2], T{-1.});
                photon(0) = -boost::math::constants::half_root_two<T>() * ellip_int_4 / (impact_root_value * (1. - 2. * u0));
                photon_time = -boost::math::constants::half_root_two<T>() *
                              (ellip_int_4 +
                               4. * ellipticIntegral(-2, u1, cubic_roots[1], T{0.}, T{1.}, cubic_roots[1], T{-1.},
                                                     -cubic_roots[0], T{1.}, cubic_roots[2], T{-1.}) +
                               2. * ellipticIntegral(-2, u0, u1, T{0.}, T{1.}, cubic_roots[1], T{-1.}, -cubic_roots[0], T{1.},
                                                     cubic_roots[2], T{-1.}) +
                               8. * ellipticIntegral(-2, u1, cubic_roots[1], T{1.}, T{-2.}, cubic_roots[1], T{-1.},
                                                     -cubic_roots[0], T{1.}, cubic_roots[2], T{-1.}) +
                               4. * ellipticIntegral(-2, u0, u1, T{1.}, T{-2.}, cubic_roots[1], T{-1.}, -cubic_roots[0], T{1.},
                                                     cubic_roots[2], T{-1.})) /
                              impact_root_value;
            } else { // u1 -> u0
                const T ellip_int_4 = ellipticIntegral(-4, u0, u1, T{0.}, T{1.}, cubic_roots[1], T{-1.}, -cubic_roots[0], T{1.},
                                                       cubic_roots[2], T{-1.});
                photon(0) = -boost::math::constants::half_root_two<T>() * ellip_int_4 / (impact_root_value * (1. - 2. * u0));
                photon_time = -boost::math::constants::half_root_two<T>() *
                              (ellip_int_4 +
                               2. * ellipticIntegral(-2, u0, u1, T{0.}, T{1.}, cubic_roots[1], T{-1.}, -cubic_roots[0], T{1.},
                                                     cubic_roots[2], T{-1.}) +
                               4. * ellipticIntegral(-2, u0, u1, T{1.}, T{-2.}, cubic_roots[1], T{-1.}, -cubic_roots[0], T{1.},
                                                     cubic_roots[2], T{-1.})) /
                              impact_root_value;
            }
        } else { // impact_parameter < std::sqrt(27)
            const T ellip_int_4 =
                ellipticIntegral2Complex(-4, u0, u1, T{0.}, T{1.}, -0.5 / (impact_root_value2 * cubic_roots[0]),
                                         cubic_roots[0] - 0.5, T{1.}, -cubic_roots[0], T{1.});
            photon(0) = -boost::math::constants::half_root_two<T>() * ellip_int_4 / (impact_root_value * (1. - 2. * u0));
            photon_time =
                -boost::math::constants::half_root_two<T>() *
                (ellip_int_4 +
                 2. * ellipticIntegral2Complex(-2, u0, u1, T{0.}, T{1.}, -0.5 / (impact_root_value2 * cubic_roots[0]),
                                               cubic_roots[0] - 0.5, T{1.}, -cubic_roots[0], T{1.}) +
                 4. * ellipticIntegral2Complex(-2, u0, u1, T{1.}, T{-2.}, -0.5 / (impact_root_value2 * cubic_roots[0]),
                                               cubic_roots[0] - 0.5, T{1.}, -cubic_roots[0], T{1.})) /
                impact_root_value;
        }
        photon(1) = r_object;
        photon(2) = theta_object;
        photon(3) = phi_object;
        photon(4) = g11_1 / (1. - 2. * u0); // dt/d\tau = 1. when r = r_observer
        if (delta_phi > turning_phi) {
            photon(5) = -g11_1 * squareRoot(1. - g11_1 * u12 * impact_root_value2);
        } else {
            photon(5) = g11_1 * squareRoot(1. - g11_1 * u12 * impact_root_value2);
        }
        // the direction of the angular momentum is [alpha * cos_theta_observer, beta, -alpha * sin_theta_observer],
        // the y component should have the same sign as the component perpendicular to [std::cos(phi), std::sin(phi), 0],
        // which is also the component along the [-std::sin(phi), std::cos(phi), 0] direction.
        if (alpha == 0.) {
            photon(6) = std::copysign(g11_1 * u12 * impact_root_value,
                                      beta * cos_phi_object - alpha * cos_theta_observer * sin_phi_object);
            photon(7) = 0.;
        } else {
            photon(6) = std::copysign(
                g11_1 * u12 * squareRoot(impact_root_value2 - power2(alpha * sin_theta_observer / sin_theta_object)),
                beta * cos_phi_object - alpha * cos_theta_observer * sin_phi_object);
            photon(7) = -g11_1 * alpha * sin_theta_observer * u12 / (sin_theta_object * sin_theta_object);
        }
        return StatusCode::SUCCESS;
    }
    int fastShadow([[maybe_unused]] const T r_observer, [[maybe_unused]] const T theta_observer,
                   [[maybe_unused]] const T sin_theta_observer, [[maybe_unused]] const T cos_theta_observer,
                   [[maybe_unused]] const T alpha, [[maybe_unused]] const T beta,
                   [[maybe_unused]] const T r_min) const override {
        // [RYT-04] formulas/ryt/ray.py
        if (alpha * alpha + beta * beta <= 27.)
            return 0;
        return 1;
    }
    // [ORB-07] formulas/orb/orbit.py
    T energy(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time, DynamicalSystem dynamics) const override {
        if (dynamics == DynamicalSystem::LAGRANGIAN) {
            if (time == TimeSystem::T)
                return (y[1] - 2.) / (y[1] * y[4]);
            // time == TimeSystem::TAU
            return (1. - 2. / y[1]) * y[4];
        } // dynamics == DynamicalSystem::HAMILTONIAN
        return 1. - y[4];
    }
    // [ORB-07] formulas/orb/orbit.py
    T angularMomentum(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time,
                      DynamicalSystem dynamics) const override {
        if (dynamics == DynamicalSystem::LAGRANGIAN) {
            if (time == TimeSystem::T)
                return power2(y[1] * std::sin(y[2])) * y[7] / y[4];
            // time == TimeSystem::TAU
            return power2(y[1] * std::sin(y[2])) * y[7];
        } // dynamics == DynamicalSystem::HAMILTONIAN
        return y[7];
    }
    // [ORB-07] formulas/orb/orbit.py
    T carterConstant(const boost::numeric::ublas::bounded_vector<T, 8> &y, [[maybe_unused]] const T mu2, TimeSystem time,
                     DynamicalSystem dynamics) const override {
        if (dynamics == DynamicalSystem::LAGRANGIAN) {
            if (time == TimeSystem::T)
                return power4(y[1]) * (power2(y[6]) + power2(y[7] * std::cos(y[2]) * std::sin(y[2]))) / power2(y[4]);
            // time == TimeSystem::TAU
            return power4(y[1]) * (power2(y[6]) + power2(y[7] * std::cos(y[2]) * std::sin(y[2])));
        } // dynamics == DynamicalSystem::HAMILTONIAN
        return power2(y[6]) + power2(y[7] / tan(y[2]));
    }
    // [ORB-07] formulas/orb/orbit.py
    int normalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T g11_1 = 1. - 2. / y[1];
        if (g11_1 <= 0)
            return 1;
        y[4] = std::sqrt(g11_1 - (power2(y[5]) / g11_1 + power2(y[1] * y[6]) + power2(y[1] * std::sin(y[2]) * y[7])));
        return std::isnan(y[4]) ? StatusCode::DOMAIN_ERROR : StatusCode::SUCCESS;
    }
    // [MET-07] formulas/met/metric.py; [ORB-07] formulas/orb/orbit.py
    int normalizeNullGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y, T frequency = 1.) const override {
        if (frequency <= 0.)
            return StatusCode::DOMAIN_ERROR;
        const T g11_1 = 1. - 2. / y[1];
        if (g11_1 <= 0)
            return 1;
        const T coefficient =
            g11_1 / std::sqrt(y[5] * y[5] + g11_1 * (power2(y[1] * y[6]) + power2(y[1] * std::sin(y[2]) * y[7])));
        if (std::isnan(coefficient))
            return StatusCode::DOMAIN_ERROR;
        y[4] = frequency;
        y[5] *= coefficient;
        y[6] *= coefficient;
        y[7] *= coefficient;
        return StatusCode::SUCCESS;
    }
    // [ORB-15] formulas/met/schwarzschild.py
    int getIntegrationSystem(std::function<void(const boost::numeric::ublas::bounded_vector<T, 8> &,
                                                boost::numeric::ublas::bounded_vector<T, 8> &, T)> &system,
                             TimeSystem time, DynamicalSystem dynamics,
                             MotionMode motion = MotionMode::GEODESIC) const override {
        if (time == TimeSystem::T) {
            if (dynamics == DynamicalSystem::LAGRANGIAN) {
                if (motion == MotionMode::GEODESIC) {
                    system = [](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                                boost::numeric::ublas::bounded_vector<T, 8> &dydt, T) {
                        dydt[0] = y[4]; // d\tau/dt
                        dydt[1] = y[5]; // dr/dt
                        dydt[2] = y[6]; // d\theta/dt
                        dydt[3] = y[7]; // d\phi/dt
                        const T r = y[1], sin_theta = sinTheta(y[2]), cos_theta = cosTheta(y[2]);
                        const T rm2 = r - 2., rm3 = r - 3.;
                        const T rm2r_1 = 1. / (rm2 * r);
                        // d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
                        dydt[4] = 2. * y[5] * rm2r_1 * y[4];
                        // d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                        dydt[5] =
                            -rm2 / power3(r) + 3. * rm2r_1 * power2(y[5]) + rm2 * (power2(y[6]) + power2(sin_theta * y[7]));
                        // d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                        dydt[6] = -2. * rm3 * rm2r_1 * y[5] * y[6] + sin_theta * cos_theta * power2(y[7]);
                        // d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                        if (y[7] == 0.)
                            dydt[7] = T{0.};
                        else if (sin_theta == 0.)
                            throw std::domain_error("Non-zero angular velocity with sin(theta) == 0.");
                        else
                            dydt[7] = -2. * (rm3 * rm2r_1 * y[5] + cos_theta / sin_theta * y[6]) * y[7];
                    };
                    return StatusCode::SUCCESS;
                } else if (motion == MotionMode::CIRCULAR) {
                    system = [](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                                boost::numeric::ublas::bounded_vector<T, 8> &dydt, T) {
                        dydt[0] = y[4]; // d\tau/dt
                        dydt[1] = 0.;   // dr/dt
                        dydt[2] = 0.;   // d\theta/dt
                        dydt[3] = y[7]; // d\phi/dt
                        dydt[4] = 0.;
                        dydt[5] = 0.;
                        dydt[6] = 0.;
                        dydt[7] = 0.;
                    };
                    return StatusCode::SUCCESS;
                } else if (motion == MotionMode::CONICAL) {
                    system = [](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                                boost::numeric::ublas::bounded_vector<T, 8> &dydt, T) {
                        const T r = y[1];
                        if (r <= 2.)
                            throw std::domain_error("Schwarzschild orbital radius r <= 2");
                        const T sin2_theta = power2(std::sin(y[2]));
                        const T g00 = -1. + 2. / r, g11 = -1. / g00, g33 = r * r * sin2_theta;
                        dydt[0] = y[4]; // d\tau/dt
                        dydt[1] = y[5]; // dr/dt
                        dydt[2] = 0.;   // d\theta/dt = 0.
                        dydt[3] = y[7]; // d\phi/dt
                        // d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt
                        dydt[4] = sin2_theta * y[5] * (r * (-g00 - g11 * y[5] * y[5] - y[4] * y[4]) + 1. + power2(g11 * y[5])) /
                                  (y[4] * (g33 + power2(g33 * y[7] / y[4])));
                        // d^2r/dt^2 = 0.
                        dydt[5] = 0.;
                        // d^2\theta/dt^2 = 0.
                        dydt[6] = 0.;
                        // d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
                        dydt[7] = y[7] * (-2. / r * y[5] + dydt[4] / y[4]);
                    };
                    return StatusCode::SUCCESS;
                } else if (motion == MotionMode::CYLINDRICAL) {
                    system = [](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                                boost::numeric::ublas::bounded_vector<T, 8> &dydt, T) {
                        const T r = y[1];
                        if (r <= 2.)
                            throw std::domain_error("Schwarzschild orbital radius r <= 2");
                        const T r_1 = 1. / r;
                        const T g11 = 1. / (1. - 2. * r_1);
                        dydt[0] = y[4]; // d\tau/dt
                        dydt[1] = y[5]; // dr/dt
                        dydt[2] = y[6]; // d\theta/dt.
                        dydt[3] = y[7]; // d\phi/dt
                        // d^2r/dt^2 = r * (d\theta/dt)^2
                        dydt[5] = r * y[6] * y[6];
                        // d^2\theta/dt^2 =
                        dydt[6] = -2. * r_1 * y[5] * y[6];
                        // d^2\phi/dt^2 = 0.
                        dydt[7] = 0.;
                        // d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt
                        dydt[4] = (r_1 * (1. + power2(g11 * y[5])) - 2. * g11 * dydt[5]) * r_1 * y[5] / y[4];
                    };
                    return StatusCode::SUCCESS;
                } else if (motion == MotionMode::CONICAL_FIXED_ENERGY) {
                    system = [](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                                boost::numeric::ublas::bounded_vector<T, 8> &dydt, T) {
                        const T r = y[1];
                        if (r <= 2.)
                            throw std::domain_error("Schwarzschild orbital radius r <= 2");
                        const T r_1 = 1. / r, r_2 = power2(r_1);
                        const T g00 = -1. + 2. * r_1, g11 = -1. / g00;
                        dydt[0] = y[4]; // d\tau/dt
                        dydt[1] = y[5]; // dr/dt
                        dydt[2] = 0.;   // d\theta/dt = 0.
                        dydt[3] = y[7]; // d\phi/dt
                        // dE/dt = 0 ==> d^2\tau/dt^2
                        dydt[4] = 2. * y[4] * y[5] * g11 * r_2;
                        // d^2\theta/dt^2 = 0.
                        dydt[6] = 0.;
                        // d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
                        dydt[7] = 2. * y[5] * y[7] * r_1 * (g11 * r_1 - 1.);
                        // Norm = -1 ==> d^2r/dt^2
                        dydt[5] = -g00 * (r_2 * (power2(g11 * y[5]) - 2. * power2(y[4]) * g11 + 1) -
                                          power2(y[7] * std::sin(y[2])) * (2. * g11 - r));
                    };
                    return StatusCode::SUCCESS;
                } else if (motion == MotionMode::CYLINDRICAL_FIXED_ENERGY) {
                    system = [](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                                boost::numeric::ublas::bounded_vector<T, 8> &dydt, T) {
                        const T r = y[1];
                        if (r <= 2.)
                            throw std::domain_error("Schwarzschild orbital radius r <= 2");
                        const T r_1 = 1. / r, r_2 = power2(r_1);
                        const T r2 = power2(r);
                        const T sin_theta = sinTheta(y[2]), cos_theta = cosTheta(y[2]);
                        const T sin2_theta = power2(sin_theta), cos2_theta = power2(cos_theta);
                        const T g00 = -1. + 2. * r_1, g11 = -1. / g00, g33 = r2 * sin2_theta;
                        dydt[0] = y[4]; // d\tau/dt
                        dydt[1] = y[5]; // dr/dt
                        dydt[2] = y[6]; // d\theta/dt.
                        dydt[3] = y[7]; // d\phi/dt
                        // dE/dt = 0 ==> d^2\tau/dt^2
                        dydt[4] = -2. * y[4] * y[5] / (g00 * r2);
                        // dL/dt = 0 ==> d^2\phi/dt^2
                        dydt[7] = dydt[4] * y[7] / y[4];
                        // d^2\theta/dt^2 =
                        dydt[6] =
                            (-y[4] * dydt[4] + y[5] * r_2 * (1. + power2(g11 * y[5])) -
                             r * y[5] * power2(y[6]) * (1. + g11 * (1. + cos2_theta) / sin2_theta) - g33 * y[7] * dydt[7]) /
                            (-r * y[5] * sin_theta / cos_theta - g11 * y[5] * r * cos_theta / sin_theta);
                        // d^2r/dt^2 = r * (d\theta/dt)^2
                        dydt[5] = r * ((1. + cos2_theta) * power2(y[6]) - sin_theta * cos_theta * dydt[6]) / sin2_theta;
                    };
                    return StatusCode::SUCCESS;
                }
            } else if (dynamics == DynamicalSystem::HAMILTONIAN && motion == MotionMode::GEODESIC) {
                system = [](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                            boost::numeric::ublas::bounded_vector<T, 8> &dydt, T) {
                    const T r_1 = 1. / y[1], r_2 = power2(r_1), g11_1 = 1. - 2. * r_1, E = 1. - y[4], L2 = power2(y[7]);
                    const T sin_1_theta = 1. / sinTheta(y[2]), sin_2_theta = power2(sin_1_theta);
                    //[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
                    dydt[0] = g11_1 / E;                          // d\tau/dt
                    dydt[1] = g11_1 * y[5] * dydt[0];             // dr/dt
                    dydt[2] = y[6] * r_2 * dydt[0];               // d\theta/dt
                    dydt[3] = y[7] * sin_2_theta * r_2 * dydt[0]; // d\phi/dt
                    dydt[4] = 0.;
                    dydt[5] =
                        (-(power2(y[5]) + power2(E) / power2(g11_1)) + (power2(y[6]) + L2 * sin_2_theta) * r_1) * r_2 * dydt[0];
                    dydt[6] = sin_2_theta * L2 * cosTheta(y[2]) * sin_1_theta * r_2 * dydt[0];
                    dydt[7] = 0.;
                };
                return StatusCode::SUCCESS;
            }
        } else if (time == TimeSystem::TAU && dynamics == DynamicalSystem::LAGRANGIAN && motion == MotionMode::GEODESIC) {
            system = [](const boost::numeric::ublas::bounded_vector<T, 8> &y, boost::numeric::ublas::bounded_vector<T, 8> &dydt,
                        T) {
                dydt[0] = y[4]; // dt/d\tau
                dydt[1] = y[5]; // dr/d\tau
                dydt[2] = y[6]; // d\theta/d\tau
                dydt[3] = y[7]; // d\phi/d\tau
                const T r = y[1], sin_theta = sinTheta(y[2]), cos_theta = cosTheta(y[2]);
                const T rm2 = r - 2., r_1 = 1. / r;
                const T rm2r_1 = 1. / (rm2 * r);
                // d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
                dydt[4] = -2. * y[5] * rm2r_1 * y[4];
                // d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                dydt[5] =
                    -rm2 * power3(r_1) * power2(y[4]) + rm2r_1 * power2(y[5]) + rm2 * (power2(y[6]) + power2(sin_theta * y[7]));
                // d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                dydt[6] = -2. * r_1 * y[5] * y[6] + sin_theta * cos_theta * power2(y[7]);
                // d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                dydt[7] = -2. * (r_1 * y[5] + cos_theta / sin_theta * y[6]) * y[7];
            };
            return StatusCode::SUCCESS;
        }
        throw std::logic_error("Schwarzschild::getIntegrationSystem() invaild");
    }
};
} // namespace sbody
