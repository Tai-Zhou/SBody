/**
 * @file kerr_eom_test.cpp
 * @author Tai Zhou
 * @brief Structural tests for ORB-16: the Kerr Lagrangian and Hamiltonian equations of motion
 *        preserve their defining invariants under numerical integration
 * @version 0.1
 * @date 2026-10-19
 *
 * @copyright Copyright (c) 2026
 *
 */

#define BOOST_TEST_MODULE kerr_eom_test

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

#include <boost/numeric/odeint.hpp>
#include <boost/test/included/unit_test.hpp>

#include "metric/kerr.hpp"

using namespace sbody;

namespace {

namespace ublas = boost::numeric::ublas;
namespace odeint = boost::numeric::odeint;

using State = ublas::bounded_vector<double, 8>;
using System = std::function<void(const State &, State &, double)>;

constexpr double SPIN = 0.5;
constexpr double MAX_DRIFT = 1e-9;

struct Invariants {
    double normalization;
    double energy;
    double angular_momentum;
    double carter;
};

Invariants invariants(const Kerr<double> &metric, const State &state) {
    ublas::bounded_matrix<double, 4, 4> tensor;
    if (metric.metricTensor(state, tensor) != StatusCode::SUCCESS)
        throw std::runtime_error("Kerr::metricTensor failed");
    const double four_velocity[4] = {1. / state[4], state[5] / state[4], state[6] / state[4],
                                     state[7] / state[4]};
    double normalization = 0.;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            normalization += tensor(i, j) * four_velocity[i] * four_velocity[j];
    return {
        normalization,
        metric.energy(state, TimeSystem::T, DynamicalSystem::LAGRANGIAN),
        metric.angularMomentum(state, TimeSystem::T, DynamicalSystem::LAGRANGIAN),
        metric.carterConstant(state, 1., TimeSystem::T, DynamicalSystem::LAGRANGIAN),
    };
}

State makeState(const Kerr<double> &metric, double radius, double angle, double radial_velocity,
                double angular_velocity, double azimuthal_velocity) {
    State state;
    state.clear();
    state[1] = radius;
    state[2] = angle;
    state[5] = radial_velocity;
    state[6] = angular_velocity;
    state[7] = azimuthal_velocity;
    if (metric.normalizeTimelikeGeodesic(state) != StatusCode::SUCCESS)
        throw std::runtime_error("Kerr::normalizeTimelikeGeodesic failed");
    return state;
}

System makeSystem(const Kerr<double> &metric, TimeSystem time, DynamicalSystem dynamics, MotionMode motion) {
    System system;
    if (metric.getIntegrationSystem(system, time, dynamics, motion) != StatusCode::SUCCESS)
        throw std::runtime_error("Kerr::getIntegrationSystem failed");
    return system;
}

struct KerrDrift {
    const Kerr<double> &metric;
    Invariants initial;
    bool check_normalization;
    bool check_energy;
    bool check_angular_momentum;
    bool check_carter;
    double normalization = 0.;
    double energy = 0.;
    double angular_momentum = 0.;
    double carter = 0.;

    void operator()(const State &state, double) {
        const Invariants current = invariants(metric, state);
        if (check_normalization)
            normalization = std::max(normalization, std::abs(current.normalization - initial.normalization));
        if (check_energy)
            energy = std::max(energy, std::abs(current.energy - initial.energy));
        if (check_angular_momentum)
            angular_momentum = std::max(angular_momentum, std::abs(current.angular_momentum - initial.angular_momentum));
        if (check_carter)
            carter = std::max(carter, std::abs(current.carter - initial.carter));
    }
};

void integrateAndCheck(const Kerr<double> &metric, MotionMode motion, const State &initial, bool check_normalization,
                       bool check_energy, bool check_angular_momentum, bool check_carter) {
    State state = initial;
    KerrDrift drift{metric, invariants(metric, initial), check_normalization, check_energy, check_angular_momentum,
                    check_carter};
    if (check_normalization)
        BOOST_REQUIRE_SMALL(drift.initial.normalization + 1., 1e-12);
    odeint::integrate_const(odeint::runge_kutta4<State>(), makeSystem(metric, TimeSystem::T,
                                                                      DynamicalSystem::LAGRANGIAN, motion),
                            state, 0., 0.2, 1e-3, drift);
    if (check_normalization)
        BOOST_REQUIRE_LT(drift.normalization, MAX_DRIFT);
    if (check_energy)
        BOOST_REQUIRE_LT(drift.energy, MAX_DRIFT);
    if (check_angular_momentum)
        BOOST_REQUIRE_LT(drift.angular_momentum, MAX_DRIFT);
    if (check_carter)
        BOOST_REQUIRE_LT(drift.carter, MAX_DRIFT);
}

struct HamiltonianDrift {
    const Kerr<double> &metric;
    double initial_pt;
    double initial_pphi;
    double null_condition = 0.;
    double pt = 0.;
    double pphi = 0.;

    void operator()(const State &state, double) {
        State lagrangian = state;
        if (metric.hamiltonianToLagrangian(lagrangian) != StatusCode::SUCCESS)
            throw std::runtime_error("Kerr::hamiltonianToLagrangian failed");
        ublas::bounded_matrix<double, 4, 4> tensor;
        if (metric.metricTensor(lagrangian, tensor) != StatusCode::SUCCESS)
            throw std::runtime_error("Kerr::metricTensor failed");
        const double four_velocity[4] = {1. / lagrangian[4], lagrangian[5] / lagrangian[4],
                                         lagrangian[6] / lagrangian[4], lagrangian[7] / lagrangian[4]};
        double value = 0.;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                value += tensor(i, j) * four_velocity[i] * four_velocity[j];
        null_condition = std::max(null_condition, std::abs(value));
        pt = std::max(pt, std::abs(state[4] - 1. - initial_pt));
        pphi = std::max(pphi, std::abs(state[7] - initial_pphi));
    }
};

} // namespace

BOOST_AUTO_TEST_CASE(TestKerrGeodesic) {
    const Kerr<double> metric(SPIN);
    integrateAndCheck(metric, MotionMode::GEODESIC, makeState(metric, 6., 1., 0.01, 0.002, 0.003), true, true,
                      true, true);
}

BOOST_AUTO_TEST_CASE(TestKerrCircular) {
    const Kerr<double> metric(SPIN);
    integrateAndCheck(metric, MotionMode::CIRCULAR,
                      makeState(metric, 6., std::acos(-1.) / 2., 0., 0., std::pow(6., -1.5)), true, true, true,
                      false);
}

BOOST_AUTO_TEST_CASE(TestKerrConicalNormalizationAndAngularMomentum) {
    const Kerr<double> metric(SPIN);
    integrateAndCheck(metric, MotionMode::CONICAL, makeState(metric, 6., 1., 0.01, 0., 0.003), true, false, true,
                      false);
}

BOOST_AUTO_TEST_CASE(TestKerrHamiltonianNullGeodesic) {
    const Kerr<double> metric(SPIN);
    State photon;
    double photon_time;
    BOOST_REQUIRE_EQUAL(metric.initializePhoton(photon, photon_time, 0.1, 0.2, 6., 36., 1., std::sin(1.)),
                        StatusCode::SUCCESS);
    BOOST_REQUIRE_EQUAL(metric.lagrangianToHamiltonian(photon), StatusCode::SUCCESS);
    HamiltonianDrift drift{metric, photon[4] - 1., photon[7]};
    odeint::integrate_const(odeint::runge_kutta4<State>(),
                            makeSystem(metric, TimeSystem::T, DynamicalSystem::HAMILTONIAN, MotionMode::GEODESIC),
                            photon, 0., 0.01, 1e-4, drift);
    BOOST_REQUIRE_LT(drift.null_condition, MAX_DRIFT);
    BOOST_REQUIRE_LT(drift.pt, MAX_DRIFT);
    BOOST_REQUIRE_LT(drift.pphi, MAX_DRIFT);
}

BOOST_AUTO_TEST_CASE(TestKerrUnsupportedModes) {
    const Kerr<double> metric(SPIN);
    System unused;
    BOOST_CHECK_THROW(metric.getIntegrationSystem(unused, TimeSystem::T, DynamicalSystem::LAGRANGIAN,
                                                  MotionMode::CYLINDRICAL),
                      std::logic_error);
    BOOST_CHECK_THROW(metric.getIntegrationSystem(unused, TimeSystem::T, DynamicalSystem::LAGRANGIAN,
                                                  MotionMode::CONICAL_FIXED_ENERGY),
                      std::logic_error);
}
