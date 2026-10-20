/**
 * @file geometry_test.cpp
 * @author Tai Zhou
 * @brief Numerical tests for the MET-10 ~ MET-12 geometry helpers.
 * @version 0.1
 * @date 2024-02-21
 *
 * @copyright Copyright (c) 2024
 *
 */

#define BOOST_TEST_MODULE geometry_test

#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/test/included/unit_test.hpp>

#include "math/geometry.hpp"

using namespace std;
using namespace sbody;

namespace {

constexpr double RELATIVE_TOLERANCE = 1e-15;
constexpr double ABSOLUTE_TOLERANCE = 1e-15;

} // namespace

BOOST_AUTO_TEST_CASE(TestRotationAngleOfPerpendicularCoordiante) {
    const double sin_theta_observer = M_SQRT1_2, cos_theta_observer = M_SQRT1_2;
    double iota, sin_iota, cos_iota, phi;
    rotationAngleOfPerpendicularCoordiante(sin_theta_observer, cos_theta_observer, 1., 0., M_SQRT1_2, -M_SQRT1_2, iota,
                                           sin_iota, cos_iota, phi);
    BOOST_CHECK_CLOSE_FRACTION(iota, std::atan2(0.5, M_SQRT1_2), RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(phi, boost::math::constants::two_thirds_pi<double>(), RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestRotationAngleOfPerpendicularCoordianteDegenerate) {
    double iota, sin_iota, cos_iota, phi;
    // Target and observer coincide: phi1 = iota = 0, and the values must not be overwritten by 0/0.
    rotationAngleOfPerpendicularCoordiante(M_SQRT1_2, M_SQRT1_2, M_SQRT1_2, M_SQRT1_2, 0., 1., iota, sin_iota, cos_iota,
                                           phi);
    BOOST_CHECK_SMALL(iota, ABSOLUTE_TOLERANCE);
    BOOST_CHECK_SMALL(sin_iota, ABSOLUTE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(cos_iota, 1., RELATIVE_TOLERANCE);
    BOOST_CHECK_SMALL(phi, ABSOLUTE_TOLERANCE);

    // Antipodal target: phi1 = pi, iota = 0.
    rotationAngleOfPerpendicularCoordiante(M_SQRT1_2, M_SQRT1_2, M_SQRT1_2, -M_SQRT1_2, 0., -1., iota, sin_iota,
                                           cos_iota, phi);
    BOOST_CHECK_SMALL(iota, ABSOLUTE_TOLERANCE);
    BOOST_CHECK_SMALL(sin_iota, ABSOLUTE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(cos_iota, 1., RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(phi, M_PI, RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestSinTheta) {
    BOOST_CHECK_SMALL(sinTheta(0.), ABSOLUTE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(sinTheta(M_PI / 4.), M_SQRT1_2, RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(sinTheta(M_PI / 2.), 1., RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(sinTheta(3. * M_PI / 4.), M_SQRT1_2, RELATIVE_TOLERANCE);
    BOOST_CHECK_SMALL(sinTheta(M_PI), ABSOLUTE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(sinTheta(-M_PI / 4.), M_SQRT1_2, RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestCosTheta) {
    BOOST_CHECK_CLOSE_FRACTION(cosTheta(0.), 1., RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(cosTheta(M_PI / 4.), M_SQRT1_2, RELATIVE_TOLERANCE);
    BOOST_CHECK_SMALL(cosTheta(M_PI / 2.), ABSOLUTE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(cosTheta(3. * M_PI / 4.), -M_SQRT1_2, RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(cosTheta(M_PI), -1., RELATIVE_TOLERANCE);
    // Negative theta is mapped to pi + theta; cos(-pi/4) = cos(3 pi/4) = -sqrt(1/2).
    BOOST_CHECK_CLOSE_FRACTION(cosTheta(-M_PI / 4.), -M_SQRT1_2, RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestDotCross) {
    namespace ublas = boost::numeric::ublas;
    ublas::bounded_vector<double, 3> x, y, z;
    x(0) = 1.;
    x(1) = 0.;
    x(2) = 0.;
    y(0) = 0.;
    y(1) = 1.;
    y(2) = 0.;
    z(0) = 0.;
    z(1) = 0.;
    z(2) = 1.;
    BOOST_CHECK_CLOSE_FRACTION(dotCross(x, y, z), 1., RELATIVE_TOLERANCE);

    x(0) = 1.;
    x(1) = 2.;
    x(2) = 3.;
    y(0) = 4.;
    y(1) = 5.;
    y(2) = 6.;
    z(0) = 7.;
    z(1) = 8.;
    z(2) = 10.;
    BOOST_CHECK_CLOSE_FRACTION(dotCross(x, y, z), -3., RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestTriangleArea) {
    namespace ublas = boost::numeric::ublas;
    ublas::bounded_vector<double, 3> x, y, z;
    x(0) = 0.;
    x(1) = 0.;
    x(2) = 0.;
    y(0) = 1.;
    y(1) = 0.;
    y(2) = 0.;
    z(0) = 0.;
    z(1) = 1.;
    z(2) = 0.;
    BOOST_CHECK_CLOSE_FRACTION(triangleArea(x, y, z), 0.5, RELATIVE_TOLERANCE);

    x(0) = 0.;
    x(1) = 0.;
    x(2) = 0.;
    y(0) = 3.;
    y(1) = 0.;
    y(2) = 0.;
    z(0) = 0.;
    z(1) = 4.;
    z(2) = 0.;
    BOOST_CHECK_CLOSE_FRACTION(triangleArea(x, y, z), 6., RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestPointInTriangle) {
    BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., 0., 0.), true);
    BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., 0.5, 0.5), true);
    BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., 0.25, 0.25), true);
    BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., 1., 1.), false);
    BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., -0.1, -0.1), false);
    BOOST_CHECK_EQUAL(pointInTriangle(0., 0., 1., 0., 0., 1., 1.1, 0.), false);
}

BOOST_AUTO_TEST_CASE(TestSphericalAngle) {
    // Same direction.
    BOOST_CHECK_SMALL(sphericalAngle(1., 0., 1., 0., 0.), ABSOLUTE_TOLERANCE);
    // Equator: the angle is delta_phi.
    BOOST_CHECK_CLOSE_FRACTION(sphericalAngle(1., 0., 1., 0., 0.5), 0.5, RELATIVE_TOLERANCE);
    // North pole to south pole.
    BOOST_CHECK_CLOSE_FRACTION(sphericalAngle(0., 1., 0., -1., 0.), M_PI, RELATIVE_TOLERANCE);
    // Generic point against the law of cosines.
    const double theta_x = 0.3, theta_y = 1.2, delta_phi = 0.5;
    const double expected =
        std::acos(std::sin(theta_x) * std::sin(theta_y) * std::cos(delta_phi) + std::cos(theta_x) * std::cos(theta_y));
    BOOST_CHECK_CLOSE_FRACTION(sphericalAngle(std::sin(theta_x), std::cos(theta_x), std::sin(theta_y), std::cos(theta_y),
                                              delta_phi),
                               expected, RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestCartesianSphericalRoundTrip) {
    namespace ublas = boost::numeric::ublas;
    ublas::bounded_vector<double, 8> cartesian, spherical, restored;
    cartesian[0] = 0.;
    cartesian[1] = 1.2;
    cartesian[2] = -0.7;
    cartesian[3] = 0.9;
    cartesian[4] = 0.;
    cartesian[5] = 0.3;
    cartesian[6] = -0.4;
    cartesian[7] = 0.5;
    cartesianToSpherical(cartesian, spherical);
    sphericalToCartesian(spherical, restored);
    for (size_t i = 0; i < 8; ++i)
        BOOST_CHECK_CLOSE_FRACTION(restored[i], cartesian[i], RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestModBy2Pi) {
    BOOST_CHECK_SMALL(modBy2Pi(0.), ABSOLUTE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(M_PI), M_PI, RELATIVE_TOLERANCE);
    BOOST_CHECK_SMALL(modBy2Pi(2. * M_PI), ABSOLUTE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(3. * M_PI), M_PI, RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(-M_PI), M_PI, RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(5. * M_PI), M_PI, RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(modBy2Pi(0.5 * M_PI), 0.5 * M_PI, RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestPhiDifference) {
    BOOST_CHECK_SMALL(phiDifference(0.), ABSOLUTE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(phiDifference(M_PI), -M_PI, RELATIVE_TOLERANCE);
    BOOST_CHECK_SMALL(phiDifference(2. * M_PI), ABSOLUTE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(phiDifference(3. * M_PI), -M_PI, RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(phiDifference(-M_PI), -M_PI, RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(phiDifference(0.5 * M_PI), 0.5 * M_PI, RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(phiDifference(1.5 * M_PI), -0.5 * M_PI, RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(phiDifference(-0.5 * M_PI), -0.5 * M_PI, RELATIVE_TOLERANCE);
}
