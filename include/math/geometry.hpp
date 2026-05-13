/**
 * @file geometry.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-03
 *
 * @copyright Copyright (c) 2026
 */

#pragma once

#include <numeric>

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fmt/core.h>

#include "math/precision.hpp"
#include "math/scalar.hpp"

namespace sbody {
enum class Axis { X, Y, Z };

/// Degree, in unit of rad
constexpr real_t deg_in_rad = 1.745329251994329576923690768488612713442871888541728e-2l;
/// Minute of arc, in unit of rad
constexpr real_t min_in_rad = 2.908882086657215961539484614147687855738119814236213e-4l;
/// Second of arc, in unit of rad
constexpr real_t sec_in_rad = 4.8481368110953599358991410235794797595635330237270218e-6l;
/// Millisecond of arc, in unit of rad
constexpr real_t mas_in_rad = 4.8481368110953599358991410235794797595635330237270218e-9l;
/// Microsecond of arc, in unit of rad
constexpr real_t uas_in_rad = 4.8481368110953599358991410235794797595635330237270218e-12l;

/**
 * @brief
 *
 * @tparam T
 * @param theta
 * @return T
 */
template <typename T>
T sinTheta(T theta) {
    return std::abs(std::sin(theta));
}

/**
 * @brief
 *
 * @tparam T
 * @param theta
 * @return T
 */
template <typename T>
T cosTheta(T theta) {
    if (theta > boost::math::constants::half_pi<T>())
        return std::cos(theta);
    if (theta < 0.)
        return -std::cos(theta);
    return std::cos(theta);
}

/**
 * @brief Dot product of vector `x` and the cross product of vector `y` and vector `z`.
 * \f$\vec{x}\cdot(\vec{y}\times\vec{z})\f$.
 *
 * @param x 3 dimensional vector
 * @param y 3 dimensional vector
 * @param z 3 dimensional vector
 * @return result
 */
template <typename T>
T dotCross(const boost::numeric::ublas::bounded_vector<T, 3> &x, const boost::numeric::ublas::bounded_vector<T, 3> &y,
           const boost::numeric::ublas::bounded_vector<T, 3> &z) {
    return x[0] * (y[1] * z[2] - y[2] * z[1]) + x[1] * (y[2] * z[0] - y[0] * z[2]) + x[2] * (y[0] * z[1] - y[1] * z[0]);
}

/**
 * @brief Calculate the area of a triangle
 *
 * @param a side length
 * @param b side length
 * @param c side length
 * @return result
 */
template <typename T>
T triangleArea(const boost::numeric::ublas::bounded_vector<T, 3> &x, const boost::numeric::ublas::bounded_vector<T, 3> &y,
               const boost::numeric::ublas::bounded_vector<T, 3> &z) {
    const T a = std::sqrt(power2(y[0] - x[0]) + power2(y[1] - x[1]) + power2(y[2] - x[2])),
            b = std::sqrt(power2(z[0] - y[0]) + power2(z[1] - y[1]) + power2(z[2] - y[2])),
            c = std::sqrt(power2(x[0] - z[0]) + power2(x[1] - z[1]) + power2(x[2] - z[2])), s = 0.5 * (a + b + c);
    return std::sqrt(s * (s - a) * (s - b) * (s - c));
}

// /**
//  * @brief Calculate the area of a triangle
//  *
//  * @param a side length
//  * @param b side length
//  * @param c side length
//  * @return result
//  */
// template <typename T>
// T TriangleArea(T a, T b, T c) {
// 	const T s = 0.5 * (a + b + c);
// 	return sqrt(s * (s - a) * (s - b) * (s - c));
// }

/**
 * @brief Calculate whether point \f$(x, y)\f$ locates inside the triangle \f$\Delta abc\f$.
 *
 * @param xa x coordiante of the vertex a
 * @param ya y coordiante of the vertex a
 * @param xb x coordiante of the vertex b
 * @param yb y coordiante of the vertex b
 * @param xc x coordiante of the vertex c
 * @param yc y coordiante of the vertex c
 * @param x x coordiante of the point
 * @param y y coordiante of the point
 * @return result
 */
template <typename T>
bool pointInTriangle(T xa, T ya, T xb, T yb, T xc, T yc, T x = 0.0, T y = 0.0) {
    const T ap_cross_ab = (x - xa) * (yb - ya) - (y - ya) * (xb - xa);
    const T bp_cross_bc = (x - xb) * (yc - yb) - (y - yb) * (xc - xb);
    const T cp_cross_ca = (x - xc) * (ya - yc) - (y - yc) * (xa - xc);
    return (ap_cross_ab >= 0. && bp_cross_bc >= 0. && cp_cross_ca >= 0.) ||
           (ap_cross_ab <= 0. && bp_cross_bc <= 0. && cp_cross_ca <= 0.);
}

template <typename T>
bool pointInTriangle(const boost::numeric::ublas::bounded_vector<T, 2> &a, const boost::numeric::ublas::bounded_vector<T, 2> &b,
                     const boost::numeric::ublas::bounded_vector<T, 2> &c) {
    return PointInTriangle(a[0], a[1], b[0], b[1], c[0], c[1]);
}

/**
 * @brief Calculate whether point \f$p\f$ locates inside the triangle \f$\Delta abc\f$.
 *
 * @param a vertex a
 * @param b vertex b
 * @param c vertex c
 * @param p point p
 * @return result
 */
template <typename T>
bool pointInTriangle(const boost::numeric::ublas::bounded_vector<T, 2> &a, const boost::numeric::ublas::bounded_vector<T, 2> &b,
                     const boost::numeric::ublas::bounded_vector<T, 2> &c,
                     const boost::numeric::ublas::bounded_vector<T, 2> &p) {
    return PointInTriangle(a[0], a[1], b[0], b[1], c[0], c[1], p[0], p[1]);
}

/**
 * @brief Rotate vector `x` around the `axis` by `angle`.
 *
 * @param x 8 dimensional vector
 * @param rotate_velocity if true, rotate velocity part, else rotate position part
 * @param axis the rotation axis
 * @param angle in rad
 * @return status
 */
template <typename T>
void rotateAroundAxis(boost::numeric::ublas::bounded_vector<T, 8> &x, bool rotate_velocity, Axis axis, T angle) {
    const T sin_angle = std::sin(angle), cos_angle = std::cos(angle);
    const std::size_t offset = rotate_velocity ? 5 : 1;
    switch (axis) {
    case Axis::X:
        angle = x(offset + 1) * cos_angle - x(offset + 2) * sin_angle;
        x(offset + 2) = x(offset + 1) * sin_angle + x(offset + 2) * cos_angle;
        x(offset + 1) = angle;
        return;
    case Axis::Y:
        angle = x(offset + 2) * cos_angle - x(offset) * sin_angle;
        x(offset) = x(offset + 2) * sin_angle + x(offset) * cos_angle;
        x(offset + 2) = angle;
        return;
    case Axis::Z:
        angle = x(offset) * cos_angle - x(offset + 1) * sin_angle;
        x(offset + 1) = x(offset) * sin_angle + x(offset + 1) * cos_angle;
        x(offset) = angle;
        return;
    }
}

/**
 * @brief
 *
 * @param cartesian 4 dimensional array
 * @param spherical 4 dimensional array
 * @return status
 */
template <typename T>
void cartesianToSpherical(const boost::numeric::ublas::bounded_vector<T, 4> &cartesian,
                          boost::numeric::ublas::bounded_vector<T, 4> &spherical) {
    if (spherical(1) = std::hypot(cartesian(1), cartesian(2), cartesian(3)); spherical(1) == 0.)
        throw std::domain_error(
            fmt::format("{}#{} {}: Cartesian coordinates (0, 0, 0) is not allowed.", __FILE__, __LINE__, __func__));
    spherical(0) = cartesian(0);
    spherical(2) = std::acos(cartesian(3) / spherical(1));
    spherical(3) = std::atan2(cartesian(2), cartesian(1));
}

/**
 * @brief
 *
 * @param cartesian 8 dimensional array
 * @param spherical 8 dimensional array
 * @return status
 */
template <typename T>
void cartesianToSpherical(const boost::numeric::ublas::bounded_vector<T, 8> &cartesian,
                          boost::numeric::ublas::bounded_vector<T, 8> &spherical) {
    if (spherical[1] = std::hypot(cartesian(1), cartesian(2), cartesian(3)); spherical[1] == 0.)
        throw std::domain_error(
            fmt::format("{}#{} {}: Cartesian coordinates (0, 0, 0) is not allowed.", __FILE__, __LINE__, __func__));
    const T r_1 = 1. / spherical[1];
    spherical[0] = cartesian[0];
    spherical[2] = std::acos(cartesian[3] * r_1);
    spherical[4] = cartesian[4];
    spherical[5] = std::inner_product(cartesian.begin() + 1, cartesian.begin() + 4, cartesian.begin() + 5, T{0.}) * r_1;
    if (const T r_xy = std::hypot(cartesian(1), cartesian(2)); r_xy == 0.) {
        spherical[3] = std::atan2(cartesian[6], cartesian[5]);
        spherical[6] = std::copysign(r_1, cartesian[3]) * std::hypot(cartesian(5), cartesian(6));
        spherical[7] = 0;
    } else {
        spherical[3] = std::atan2(cartesian[2], cartesian[1]);
        spherical[6] = (-cartesian[7] + cartesian[3] * r_1 * spherical[5]) / r_xy;
        spherical[7] = (cartesian[6] * cartesian[1] - cartesian[5] * cartesian[2]) / power2(r_xy);
    }
}

/**
 * @brief
 *
 * @param x 4 or 8 dimensional vector
 * @return status
 */
template <typename T, std::size_t N>
void cartesianToSpherical(boost::numeric::ublas::bounded_vector<T, N> &x) {
    static_assert(N == 4 || N == 8, "Invalid dimension");
    const boost::numeric::ublas::bounded_vector<T, N> cartesian(x);
    cartesianToSpherical(cartesian, x);
}

/**
 * @brief
 *
 * @param spherical 4 or 8 dimensional vector
 * @param cartesian 4 or 8 dimensional vector
 * @param dimension 4 or 8.
 * @return status
 */
template <typename T>
void sphericalToCartesian(const boost::numeric::ublas::bounded_vector<T, 4> &spherical,
                          boost::numeric::ublas::bounded_vector<T, 4> &cartesian) {
    const T sin_theta = sinTheta(spherical[2]), cos_theta = cosTheta(spherical[2]);
    cartesian[0] = spherical[0];
    cartesian[1] = spherical[1] * sin_theta * std::cos(spherical[3]);
    cartesian[2] = spherical[1] * sin_theta * std::sin(spherical[3]);
    cartesian[3] = spherical[1] * cos_theta;
}

template <typename T>
void sphericalToCartesian(const boost::numeric::ublas::bounded_vector<T, 8> &spherical,
                          boost::numeric::ublas::bounded_vector<T, 8> &cartesian) {
    const T sin_theta = sinTheta(spherical[2]), cos_theta = cosTheta(spherical[2]), sin_phi = std::sin(spherical[3]),
            cos_phi = std::cos(spherical[3]);
    cartesian[0] = spherical[0];
    cartesian[1] = spherical[1] * sin_theta * cos_phi;
    cartesian[2] = spherical[1] * sin_theta * sin_phi;
    cartesian[3] = spherical[1] * cos_theta;
    cartesian[4] = spherical[4];
    cartesian[5] = spherical[5] * sin_theta * cos_phi +
                   spherical[1] * (cos_theta * cos_phi * spherical[6] - sin_theta * sin_phi * spherical[7]);
    cartesian[6] = spherical[5] * sin_theta * sin_phi +
                   spherical[1] * (cos_theta * sin_phi * spherical[6] + sin_theta * cos_phi * spherical[7]);
    cartesian[7] = spherical[5] * cos_theta - spherical[1] * sin_theta * spherical[6];
}

/**
 * @brief
 *
 * @param x 4 or 8 dimensional vector
 * @return status
 */
template <typename T, std::size_t N>
void sphericalToCartesian(boost::numeric::ublas::bounded_vector<T, N> &x) {
    static_assert(N == 4 || N == 8, "Invalid dimension");
    const boost::numeric::ublas::bounded_vector<T, N> cartesian(x);
    sphericalToCartesian(cartesian, x);
}

template <typename T>
T sphericalAngle(T cos_theta_x, T cos_theta_y, T delta_theta_xy, T delta_phi_xy) {
    const T a = power2(std::sin(0.5 * delta_theta_xy)) +
                cos_theta_x * cos_theta_y * power2(std::sin(0.5 * phiDifference(delta_phi_xy)));
    return 2. * std::atan2(squareRoot(a), squareRoot(1. - a));
}

template <typename T>
void rotationAngleOfPerpendicularCoordiante(T sin_theta_observer, T cos_theta_observer, T sin_theta0, T cos_theta0, T sin_phi0,
                                            T cos_phi0, T &iota, T &sin_iota, T &cos_iota, T &phi1) {
    if (cos_theta_observer == cos_theta0 && cos_phi0 == 1.0) {
        phi1 = 0.0;
        iota = 0.0;
        sin_iota = 0.0;
        cos_iota = 1.0;
    }
    if (cos_theta_observer == -cos_theta0 && cos_phi0 == -1.0) {
        phi1 = M_PI;
        iota = 0.0;
        sin_iota = 0.0;
        cos_iota = 1.0;
    }
    const T x = sin_theta0 * cos_phi0;
    const T y = sin_theta0 * sin_phi0;
    // const T z = cos_theta0;
    const T x1 = sin_theta_observer * x + cos_theta_observer * cos_theta0;
    // const T y1 = y;
    const T z1 = -cos_theta_observer * x + sin_theta_observer * cos_theta0;
    const T yz1 = std::hypot(y, z1); // yz1 >= 0.0
    sin_iota = z1 / yz1;
    cos_iota = y / yz1;
    iota = std::atan2(sin_iota, cos_iota);
    phi1 = std::atan2(yz1, x1); // 0.0 <= phi1 <= pi
}

template <typename T>
void rotateSphericalAxis(T sin_theta_observer, T cos_theta_observer, T sin_theta0, T cos_theta0, T sin_phi0, T cos_phi0,
                         T &theta1, T &sin_theta1, T &cos_theta1, T &phi1, T &sin_phi1, T &cos_phi1) {
    const T x = sin_theta0 * cos_phi0;
    const T y = sin_theta0 * sin_phi0;
    const T x1 = x * cos_theta_observer - cos_theta0 * sin_theta_observer;
    sin_theta1 = std::hypot(x1, y);
    cos_theta1 = x * sin_theta_observer + cos_theta0 * cos_theta_observer;
    theta1 = std::atan2(sin_theta1, cos_theta1);
    sin_phi1 = y / sin_theta1;
    cos_phi1 = x1 / sin_theta1;
    phi1 = std::atan2(y, x1);
}

template <typename T>
void rotateSphericalAxisWithIota(T sin_theta_observer, T cos_theta_observer, T sin_iota, T cos_iota, T sin_theta0, T cos_theta0,
                                 T sin_phi0, T cos_phi0, T &theta1, T &sin_theta1, T &cos_theta1, T &phi1, T &sin_phi1,
                                 T &cos_phi1) {
    const T x = sin_theta0 * cos_phi0;
    const T y = sin_theta0 * sin_phi0;
    // const T z = cos_theta0;
    const T x1 = x * sin_theta_observer + cos_theta0 * cos_theta_observer;
    const T z1 = cos_theta0 * sin_theta_observer - x * cos_theta_observer;
    const T y2 = y * cos_iota + z1 * sin_iota;
    sin_theta1 = std::hypot(x1, y2);
    cos_theta1 = z1 * cos_iota - y * sin_iota;
    theta1 = std::atan2(sin_theta1, cos_theta1);
    sin_phi1 = y2 / sin_theta1;
    cos_phi1 = x1 / sin_theta1;
    phi1 = std::atan2(y2, x1);
    if (phi1 < 0.0)
        phi1 += boost::math::constants::two_pi<T>();
}

template <typename T>
void rotateSphericalAxisToLineOfSight(T sin_theta_observer, T cos_theta_observer,
                                      const boost::numeric::ublas::bounded_vector<T, 2> &source,
                                      const boost::numeric::ublas::bounded_vector<T, 2> &target, bool parallel) {
    T sin_theta_obj = sinTheta(source[0]);
    T x = sin_theta_obj * std::cos(source[1]);
    T y = sin_theta_obj * std::sin(source[1]);
    T z = cosTheta(source[0]);
    T x_target, z_target;
    if (parallel) {
        x_target = cos_theta_observer * x - sin_theta_observer * z;
        z_target = sin_theta_observer * x + cos_theta_observer * z;
    } else {
        x_target = sin_theta_observer * x + cos_theta_observer * z;
        z_target = -cos_theta_observer * x + sin_theta_observer * z;
    }
    target[0] = std::atan2(std::hypot(x_target, y), z_target);
    target[1] = std::atan2(y, x_target);
}

/**
 * @brief The modified spherical coordiante maps the polar angle from \f$\theta\in[0,\pi]\f$ to \f[
 * \theta'\equiv\begin{cases}\theta & (\theta\leq\pi/2)\\
 * \theta-\pi & (\theta>\pi/2)\end{cases}.
 * \f]
 * @param y 8 dimensional vector
 * @param last_theta theta in the last integration step
 */
template <typename T>
void mapTheta(boost::numeric::ublas::bounded_vector<T, 8> &y, T &last_theta) {
    if (oppositeSign(last_theta, y[2])) {
        y[2] = -y[2];
        y[3] += boost::math::constants::pi<T>();
        y[6] = -y[6];
    } else if (y[2] <= -boost::math::constants::half_pi<T>())
        y[2] += boost::math::constants::pi<T>();
    else if (y[2] > boost::math::constants::half_pi<T>())
        y[2] -= boost::math::constants::pi<T>();
    last_theta = y[2];
}

template <typename Stepper, typename T>
void denseStepperMapTheta(Stepper &stepper, T &last_theta) {
    T current_theta = stepper.current_state()[2];
    if (oppositeSign(last_theta, current_theta)) {
        boost::numeric::ublas::bounded_vector<T, 8> refreshed_state = stepper.current_state();
        refreshed_state[2] = -refreshed_state[2];
        refreshed_state[3] += boost::math::constants::pi<T>();
        refreshed_state[6] = -refreshed_state[6];
        stepper.initialize(refreshed_state, stepper.current_time(), stepper.current_time_step());
        last_theta = refreshed_state[2];
    } else if (current_theta <= -boost::math::constants::half_pi<T>()) {
        boost::numeric::ublas::bounded_vector<T, 8> refreshed_state = stepper.current_state();
        refreshed_state[2] += boost::math::constants::pi<T>();
        stepper.initialize(refreshed_state, stepper.current_time(), stepper.current_time_step());
        last_theta = refreshed_state[2];
    } else if (current_theta > boost::math::constants::half_pi<T>()) {
        boost::numeric::ublas::bounded_vector<T, 8> refreshed_state = stepper.current_state();
        refreshed_state[2] -= boost::math::constants::pi<T>();
        stepper.initialize(refreshed_state, stepper.current_time(), stepper.current_time_step());
        last_theta = refreshed_state[2];
    } else
        last_theta = current_theta;
}

/**
 * @brief Return `phi` in \f$[0, 2\pi)\f$.
 *
 * @param phi \f$\phi\f$
 * @return result
 */
template <typename T>
T modBy2Pi(T phi) {
    return phi - floor(phi / boost::math::constants::two_pi<T>()) * boost::math::constants::two_pi<T>();
}

/**
 * @brief Similar to `ModBy2Pi`, but return `phi` in \f$[-\pi, \pi)\f$.
 *
 * @param phi \f$\phi\f$
 * @return result
 */
template <typename T>
T phiDifference(T phi) {
    return phi - floor(phi / boost::math::constants::two_pi<T>() + 0.5) * boost::math::constants::two_pi<T>();
}
} // namespace sbody
