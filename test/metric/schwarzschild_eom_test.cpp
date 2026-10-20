/**
 * @file schwarzschild_eom_test.cpp
 * @author Tai Zhou
 * @brief Structural tests for ORB-15: the six TimeSystem::T Lagrangian equations of motion
 *        preserve their defining invariants under numerical integration
 * @version 0.1
 * @date 2026-10-19
 *
 * @copyright Copyright (c) 2026
 *
 */

#define BOOST_TEST_MODULE schwarzschild_eom_test

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

#include <boost/numeric/odeint.hpp>
#include <boost/test/included/unit_test.hpp>

#include "metric/schwarzschild.hpp"

using namespace sbody;

namespace {

namespace ublas = boost::numeric::ublas;
namespace odeint = boost::numeric::odeint;

using State = ublas::bounded_vector<double, 8>;

constexpr double MAX_DRIFT = 1e-10;

struct Invariants {
    double normalization;
    double energy;
    double angular_momentum;
    double cylindrical_radius;
};

Invariants invariants(const State &state) {
    const double radius = state[1], angle = state[2], dtau_dt = state[4];
    const double radial_velocity = state[5], angular_velocity = state[6], azimuthal_velocity = state[7];
    const double f = 1. - 2. / radius;
    const double sin_angle = std::sin(angle);
    return {
        (-f + radial_velocity * radial_velocity / f + (radius * angular_velocity) * (radius * angular_velocity) +
         (radius * sin_angle * azimuthal_velocity) * (radius * sin_angle * azimuthal_velocity)) /
            (dtau_dt * dtau_dt),
        f / dtau_dt,
        (radius * sin_angle) * (radius * sin_angle) * azimuthal_velocity / dtau_dt,
        radius * sin_angle,
    };
}

State makeState(double radius, double angle, double radial_velocity, double angular_velocity,
                double azimuthal_velocity) {
    State state;
    state.clear();
    state[1] = radius;
    state[2] = angle;
    state[5] = radial_velocity;
    state[6] = angular_velocity;
    state[7] = azimuthal_velocity;
    const double f = 1. - 2. / radius;
    const double dtau_dt_squared =
        f - radial_velocity * radial_velocity / f - (radius * angular_velocity) * (radius * angular_velocity) -
        (radius * std::sin(angle) * azimuthal_velocity) * (radius * std::sin(angle) * azimuthal_velocity);
    state[4] = std::sqrt(dtau_dt_squared);
    return state;
}

using System = std::function<void(const State &, State &, double)>;

System makeSystem(MotionMode motion) {
    System system;
    const Schwarzschild<double> metric;
    if (metric.getIntegrationSystem(system, TimeSystem::T, DynamicalSystem::LAGRANGIAN, motion) !=
        StatusCode::SUCCESS)
        throw std::runtime_error("Schwarzschild::getIntegrationSystem failed");
    return system;
}

struct Drift {
    Invariants initial;
    bool check_normalization;
    bool check_energy;
    bool check_angular_momentum;
    bool check_cylindrical_radius;
    double normalization = 0.;
    double energy = 0.;
    double angular_momentum = 0.;
    double cylindrical_radius = 0.;

    void operator()(const State &state, double) {
        const Invariants current = invariants(state);
        if (check_normalization)
            normalization = std::max(normalization, std::abs(current.normalization - initial.normalization));
        if (check_energy)
            energy = std::max(energy, std::abs(current.energy - initial.energy));
        if (check_angular_momentum)
            angular_momentum = std::max(angular_momentum, std::abs(current.angular_momentum - initial.angular_momentum));
        if (check_cylindrical_radius)
            cylindrical_radius =
                std::max(cylindrical_radius, std::abs(current.cylindrical_radius - initial.cylindrical_radius));
    }
};

void integrateAndCheck(MotionMode motion, const State &initial, bool check_normalization, bool check_energy,
                       bool check_angular_momentum, bool check_cylindrical_radius) {
    State state = initial;
    Drift drift{invariants(initial), check_normalization, check_energy, check_angular_momentum,
                check_cylindrical_radius};
    if (check_normalization)
        BOOST_REQUIRE_SMALL(drift.initial.normalization + 1., 1e-12);
    odeint::integrate_const(odeint::runge_kutta4<State>(), makeSystem(motion), state, 0., 0.2, 1e-3, drift);
    if (check_normalization)
        BOOST_REQUIRE_LT(drift.normalization, MAX_DRIFT);
    if (check_energy)
        BOOST_REQUIRE_LT(drift.energy, MAX_DRIFT);
    if (check_angular_momentum)
        BOOST_REQUIRE_LT(drift.angular_momentum, MAX_DRIFT);
    if (check_cylindrical_radius)
        BOOST_REQUIRE_LT(drift.cylindrical_radius, MAX_DRIFT);
}

} // namespace

BOOST_AUTO_TEST_CASE(TestGeodesic) {
    integrateAndCheck(MotionMode::GEODESIC, makeState(6., 1., 0.01, 0.002, 0.003), true, true, true, false);
}

BOOST_AUTO_TEST_CASE(TestCircular) {
    integrateAndCheck(MotionMode::CIRCULAR, makeState(6., std::acos(-1.) / 2., 0., 0., std::pow(6., -1.5)), true, true,
                      true, false);
}

// The original conical and cylindrical modes are judged mainly by preserving g(u, u) = -1; their
// energy/angular-momentum behavior is not part of this structural assertion.
BOOST_AUTO_TEST_CASE(TestConicalNormalization) {
    integrateAndCheck(MotionMode::CONICAL, makeState(6., 1., 0.01, 0., 0.003), true, false, false, false);
}

BOOST_AUTO_TEST_CASE(TestCylindricalNormalization) {
    constexpr double RADIUS = 6., ANGLE = 1., RADIAL_VELOCITY = 0.01, AZIMUTHAL_VELOCITY = 0.003;
    const double angular_velocity = -RADIAL_VELOCITY * std::tan(ANGLE) / RADIUS;
    integrateAndCheck(MotionMode::CYLINDRICAL,
                      makeState(RADIUS, ANGLE, RADIAL_VELOCITY, angular_velocity, AZIMUTHAL_VELOCITY), true, false,
                      false, false);
}

BOOST_AUTO_TEST_CASE(TestConicalFixedEnergy) {
    integrateAndCheck(MotionMode::CONICAL_FIXED_ENERGY, makeState(6., 1., 0.01, 0., 0.003), true, true, true,
                      false);
}

BOOST_AUTO_TEST_CASE(TestCylindricalFixedEnergy) {
    constexpr double RADIUS = 6., ANGLE = 1., RADIAL_VELOCITY = 0.01, AZIMUTHAL_VELOCITY = 0.003;
    const double angular_velocity = -RADIAL_VELOCITY * std::tan(ANGLE) / RADIUS;
    integrateAndCheck(MotionMode::CYLINDRICAL_FIXED_ENERGY,
                      makeState(RADIUS, ANGLE, RADIAL_VELOCITY, angular_velocity, AZIMUTHAL_VELOCITY), true, true,
                      true, true);
}

BOOST_AUTO_TEST_CASE(TestDomainGuards) {
    System system = makeSystem(MotionMode::CONICAL);
    State state = makeState(6., 1., 0.01, 0., 0.003);
    State derivative;
    state[1] = 2.;
    BOOST_CHECK_THROW(system(state, derivative, 0.), std::domain_error);

    system = makeSystem(MotionMode::GEODESIC);
    state = makeState(6., 0., 0.01, 0., 0.01);
    BOOST_CHECK_THROW(system(state, derivative, 0.), std::domain_error);

    System unused;
    BOOST_CHECK_THROW(Schwarzschild<double>().getIntegrationSystem(
                          unused, TimeSystem::T, DynamicalSystem::LAGRANGIAN, MotionMode::MOTION_MODE_NUM),
                      std::logic_error);
}
