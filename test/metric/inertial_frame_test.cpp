/**
 * @file inertial_frame_test.cpp
 * @author Tai Zhou
 * @brief Numerical comparison for MET-05; reference values come from sympy via
 *        formulas/met/metric.py
 * @version 0.1
 * @date 2026-10-14
 *
 * @copyright Copyright (c) 2026
 *
 */

#define BOOST_TEST_MODULE inertial_frame_test

#include <cmath>

#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/test/included/unit_test.hpp>

#include "metric/kerr.hpp"
#include "metric/schwarzschild.hpp"

using namespace sbody;

namespace {

namespace ublas = boost::numeric::ublas;

constexpr double RELATIVE_TOLERANCE = 1e-13;
constexpr double ABSOLUTE_TOLERANCE = 1e-14;

// [MET-05] reference values, printed by `uv run python -m formulas.met.metric`: the local inertial
// frame of the circular Schwarzschild orbit at r = 6, theta = pi/2, dtau/dt = 1/sqrt(2) and
// dphi/dt = 1/(6 sqrt 6). The rows are eps^mu_a, a = 0..3.
constexpr double FRAME_RADIUS = 6.;
constexpr double FRAME_THETA = boost::math::constants::half_pi<double>();
constexpr double FRAME_OMEGA = 0.06804138174397717;
constexpr double FRAME_REFERENCE[4][4] = {
    {1.4142135623730951, 0., 0., 0.09622504486493763},
    {0., 0.816496580927726, 0., 0.},
    {0., 0., 0.16666666666666666, 0.},
    {0.7071067811865476, 0., 0., 0.19245008972987526},
};

// The circular orbit at `radius` in the equatorial plane, in the TimeSystem::T convention:
// [tau, r, theta, phi, dtau/dt, dr/dt, dtheta/dt, dphi/dt]
ublas::bounded_vector<double, 8> makeCircularOrbitState(double radius) {
    ublas::bounded_vector<double, 8> position;
    position.clear();
    position[1] = radius;
    position[2] = boost::math::constants::half_pi<double>();
    position[4] = std::sqrt(1. - 3. / radius); // dtau/dt = 1/u^t, u^t = 1/sqrt(1 - 3M/r)
    position[7] = 1. / (radius * std::sqrt(radius)); // dphi/dt of a circular orbit
    return position;
}

// The same state in the TimeSystem::TAU convention: slot 0 is the observer time and slot 4 owns the
// four-velocity, so every velocity is scaled by u^t = 1/(dtau/dt).
ublas::bounded_vector<double, 8> toTimeSystemTau(const ublas::bounded_vector<double, 8> &state) {
    ublas::bounded_vector<double, 8> position;
    position.clear();
    position[0] = 0.;
    position[1] = state[1];
    position[2] = state[2];
    position[3] = state[3];
    const double u_t = 1. / state[4];
    position[4] = u_t;
    position[5] = u_t * state[5];
    position[6] = u_t * state[6];
    position[7] = u_t * state[7];
    return position;
}

// g(eps_i, eps_j) must be the flat metric: the frame is orthonormal and the first vector is the
// (timelike, future directed) four-velocity direction
template <typename T>
void checkFrameIsOrthonormal(const Metric<T> &metric, const ublas::bounded_vector<T, 8> &position) {
    ublas::bounded_matrix<T, 4, 4> frame;
    BOOST_REQUIRE_EQUAL(metric.localInertialFrame(position, TimeSystem::T, frame), StatusCode::SUCCESS);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const ublas::matrix_row<const ublas::bounded_matrix<T, 4, 4>> row_i(frame, i), row_j(frame, j);
            const T product = metric.dotProduct(position, row_i, row_j, true);
            if (i == j)
                BOOST_CHECK_CLOSE_FRACTION(product, i == 0 ? T{-1.} : T{1.}, RELATIVE_TOLERANCE);
            else
                BOOST_CHECK_SMALL(product, ABSOLUTE_TOLERANCE);
        }
    }
    BOOST_CHECK_GT(frame(0, 0), T{0.});
}

} // namespace

BOOST_AUTO_TEST_CASE(TestSchwarzschildFrameReference) {
    const Schwarzschild<double> schwarzschild;
    ublas::bounded_matrix<double, 4, 4> frame;
    BOOST_REQUIRE_EQUAL(schwarzschild.localInertialFrame(makeCircularOrbitState(FRAME_RADIUS), TimeSystem::T, frame),
                        StatusCode::SUCCESS);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            BOOST_CHECK_SMALL(frame(i, j) - FRAME_REFERENCE[i][j], ABSOLUTE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(frame(0, 0), std::sqrt(2.), RELATIVE_TOLERANCE); // the row is w/sqrt(-g(w, w))
}

BOOST_AUTO_TEST_CASE(TestFrameOrthonormality) {
    const Schwarzschild<double> schwarzschild;
    checkFrameIsOrthonormal(schwarzschild, makeCircularOrbitState(6.));
    checkFrameIsOrthonormal(schwarzschild, makeCircularOrbitState(20.));

    // the spin only enters through metricTensor / dotProduct, so the same algorithm must work at a > 0
    const Kerr<double> kerr(0.9);
    ublas::bounded_vector<double, 8> position;
    position.clear();
    position[1] = 20.;
    position[2] = 1.;
    position[4] = 1.;
    position[7] = 1. / 20.;
    checkFrameIsOrthonormal(kerr, position);
}

BOOST_AUTO_TEST_CASE(TestTimeSystemConvention) {
    const Schwarzschild<double> schwarzschild;
    const ublas::bounded_vector<double, 8> state = makeCircularOrbitState(FRAME_RADIUS);
    ublas::bounded_matrix<double, 4, 4> frame_t, frame_tau;
    BOOST_REQUIRE_EQUAL(schwarzschild.localInertialFrame(state, TimeSystem::T, frame_t), StatusCode::SUCCESS);
    BOOST_REQUIRE_EQUAL(schwarzschild.localInertialFrame(toTimeSystemTau(state), TimeSystem::TAU, frame_tau),
                        StatusCode::SUCCESS);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            BOOST_CHECK_SMALL(frame_t(i, j) - frame_tau(i, j), ABSOLUTE_TOLERANCE);
}
