/**
 * @file kerr_metric_test.cpp
 * @author Tai Zhou
 * @brief Numerical comparison for MET-01/MET-02/MET-08 and the ORB-06/ORB-07 positive-rate
 *        guard; reference values come from sympy via formulas/met/metric.py
 * @version 0.1
 * @date 2026-10-14
 *
 * @copyright Copyright (c) 2026
 *
 */

#define BOOST_TEST_MODULE kerr_metric_test

#include <cmath>

#include <boost/test/included/unit_test.hpp>

#include "metric/kerr.hpp"
#include "metric/schwarzschild.hpp"

using namespace sbody;

namespace {

constexpr double RELATIVE_TOLERANCE = 1e-13;

// Standard parameter points, see CPP_POINTS in formulas/common/checks.py
struct MetricReference {
    double spin;
    double radius;
    double theta;
    double g_00;
    double g_03;
    double g_11;
    double g_22;
    double g_33;
};

struct AuxiliaryReference {
    double spin;
    double radius;
    double theta;
    double u_plus_1;
    double u_plus;
    double u_r;
    double delta;
    double rho2;
};

struct CarterReference {
    double spin;
    double l;
    double q2;
    double cos_theta;
    double sin_theta;
    double alpha;
    double beta;
};

struct IscoReference {
    double spin;
    double prograde;
    double retrograde;
};

// [MET-01] reference values, printed by `uv run python -m formulas.met.metric`
const MetricReference METRIC_REFERENCE[] = {
    {0.5, 6.0, 1.5707963267948966, -0.6666666666666666, -0.16666666666666666, 1.4845360824742269, 36.0, 36.333333333333336},
    {0.5, 20.0, 1.0, -0.9000182420830148, -0.035397212546640196, 1.11054262774582, 400.0729816454316, 283.4189175766395},
    {0.9, 6.0, 1.0, -0.6688418288075241, -0.21103586843893282, 1.4605586671180333, 36.23646053119841, 26.198668526519654},
    {0.9, 20.0, 1.5707963267948966, -0.9, -0.09, 1.1086167234832738, 400.0, 400.891},
};

// [MET-02] reference values, printed by `uv run python -m formulas.met.metric`
const AuxiliaryReference AUXILIARY_REFERENCE[] = {
    {0.5, 6.0, 1.5707963267948966, 1.8660254037844386, 0.5358983848622454, -0.5773502691896257, 24.25, 36.0},
    {0.5, 20.0, 1.0, 1.8660254037844386, 0.5358983848622454, -0.5773502691896257, 360.25, 400.0729816454316},
    {0.9, 6.0, 1.0, 1.4358898943540674, 0.6964322291925095, -1.1470786693528088, 24.81, 36.23646053119841},
    {0.9, 20.0, 1.5707963267948966, 1.4358898943540674, 0.6964322291925095, -1.1470786693528088, 360.81, 400.0},
};

// [MET-08] reference values, printed by `uv run python -m formulas.met.metric`
const CarterReference CARTER_REFERENCE[] = {
    {0.0, 0.3, 0.2, 0.5, 0.8660254037844386, -0.34641016151377546, 0.41231056256176607},
    {0.5, 0.3, 0.2, 0.5, 0.8660254037844386, -0.34641016151377546, 0.4821825380496477},
    {0.9, 0.3, 0.2, 0.5, 0.8660254037844386, -0.34641016151377546, 0.6103277807866851},
};

// [ORB-17] reference values, printed by `uv run python -m formulas.orb.orbit`
const IscoReference ISCO_REFERENCE[] = {
    {0.5, 4.2330025295308256773, 7.5545847145123584417},
    {0.9, 2.3208830417618872468, 8.7173522796064893159},
};

boost::numeric::ublas::bounded_vector<double, 8> makePosition(double radius, double theta) {
    boost::numeric::ublas::bounded_vector<double, 8> position;
    position.clear();
    position[1] = radius;
    position[2] = theta;
    return position;
}

// u_plus_1 and friends are protected members of Kerr; this test-local subclass only exposes them
class KerrProbe : public Kerr<double> {
  public:
    explicit KerrProbe(double spin) : Kerr<double>(spin) {}
    double uPlus1() const { return u_plus_1; }
    double uPlus() const { return u_plus; }
    double uMinus1() const { return u_minus_1; }
    double uMinus() const { return u_minus; }
    double uR() const { return u_r; }
};

} // namespace

BOOST_AUTO_TEST_CASE(TestKerrMetricTensor) {
    for (const MetricReference &reference : METRIC_REFERENCE) {
        const Kerr<double> kerr(reference.spin);
        boost::numeric::ublas::bounded_matrix<double, 4, 4> metric;
        BOOST_REQUIRE_EQUAL(kerr.metricTensor(makePosition(reference.radius, reference.theta), metric), StatusCode::SUCCESS);
        BOOST_CHECK_CLOSE_FRACTION(metric(0, 0), reference.g_00, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(metric(0, 3), reference.g_03, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(metric(3, 0), reference.g_03, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(metric(1, 1), reference.g_11, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(metric(2, 2), reference.g_22, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(metric(3, 3), reference.g_33, RELATIVE_TOLERANCE);
    }
}

BOOST_AUTO_TEST_CASE(TestKerrMetricTensorStructure) {
    const boost::numeric::ublas::bounded_vector<double, 8> position = makePosition(6., 1.);
    Kerr<double> kerr(0.9);
    boost::numeric::ublas::bounded_matrix<double, 4, 4> metric;
    BOOST_REQUIRE_EQUAL(kerr.metricTensor(position, metric), StatusCode::SUCCESS);
    // Only six entries are assigned after clear(), the remaining ones must be exactly 0
    BOOST_CHECK_EQUAL(metric(0, 1), 0.);
    BOOST_CHECK_EQUAL(metric(0, 2), 0.);
    BOOST_CHECK_EQUAL(metric(1, 0), 0.);
    BOOST_CHECK_EQUAL(metric(1, 2), 0.);
    BOOST_CHECK_EQUAL(metric(1, 3), 0.);
    BOOST_CHECK_EQUAL(metric(2, 3), 0.);

    // a = 0 must agree entry by entry with the independent Schwarzschild implementation
    const Kerr<double> non_spinning(0.);
    const Schwarzschild<double> schwarzschild;
    boost::numeric::ublas::bounded_matrix<double, 4, 4> kerr_metric, schwarzschild_metric;
    BOOST_REQUIRE_EQUAL(non_spinning.metricTensor(position, kerr_metric), StatusCode::SUCCESS);
    BOOST_REQUIRE_EQUAL(schwarzschild.metricTensor(position, schwarzschild_metric), StatusCode::SUCCESS);
    BOOST_CHECK_CLOSE_FRACTION(kerr_metric(0, 0), schwarzschild_metric(0, 0), RELATIVE_TOLERANCE);
    BOOST_CHECK_EQUAL(kerr_metric(0, 3), 0.);
    BOOST_CHECK_CLOSE_FRACTION(kerr_metric(1, 1), schwarzschild_metric(1, 1), RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(kerr_metric(2, 2), schwarzschild_metric(2, 2), RELATIVE_TOLERANCE);
    BOOST_CHECK_CLOSE_FRACTION(kerr_metric(3, 3), schwarzschild_metric(3, 3), RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestKerrAuxiliaryQuantities) {
    for (const AuxiliaryReference &reference : AUXILIARY_REFERENCE) {
        KerrProbe kerr(reference.spin);
        BOOST_CHECK_CLOSE_FRACTION(kerr.uPlus1(), reference.u_plus_1, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(kerr.uPlus(), reference.u_plus, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(kerr.uR(), reference.u_r, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(kerr.getDelta(reference.radius), reference.delta, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(kerr.getRho2(reference.radius, reference.theta), reference.rho2, RELATIVE_TOLERANCE);
    }
}

BOOST_AUTO_TEST_CASE(TestKerrInnerHorizonPole) {
    // [MET-02] The u-integral decomposition needs the two poles u_+- = 1/r_+- (the roots of
    // 1 - 2u + a^2 u^2 = 0): in the code u_plus = 1/r_+ (smaller root) and
    // u_minus = 1/r_- = r_+/a^2 (larger root), while u_minus_1 is the reciprocal of u_-, i.e.
    // umi = 1/u_- = r_- of geokerr.f:1575.
    for (const AuxiliaryReference &reference : AUXILIARY_REFERENCE) {
        const KerrProbe kerr(reference.spin);
        const double inner_horizon_radius = 1. - std::sqrt(1. - reference.spin * reference.spin);
        const double spin2 = reference.spin * reference.spin;
        BOOST_CHECK_CLOSE_FRACTION(kerr.uMinus1(), spin2 * kerr.uPlus(), RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(kerr.uMinus1(), inner_horizon_radius, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(kerr.uMinus1(), spin2 / (1. + std::sqrt(1. - spin2)), RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(kerr.uMinus(), 1. / inner_horizon_radius, RELATIVE_TOLERANCE);
    }
}

BOOST_AUTO_TEST_CASE(TestCarterConstantToAlphaBeta) {
    for (const CarterReference &reference : CARTER_REFERENCE) {
        double alpha = 0., beta = 0.;
        int status;
        if (reference.spin == 0.) {
            const Schwarzschild<double> schwarzschild;
            status = schwarzschild.angularMomentumCarterConstantToAlphaBeta(reference.l, reference.q2, reference.cos_theta,
                                                                            reference.sin_theta, alpha, beta);
        } else {
            const Kerr<double> kerr(reference.spin);
            status = kerr.angularMomentumCarterConstantToAlphaBeta(reference.l, reference.q2, reference.cos_theta,
                                                                   reference.sin_theta, alpha, beta);
        }
        BOOST_REQUIRE_EQUAL(status, StatusCode::SUCCESS);
        BOOST_CHECK_CLOSE_FRACTION(alpha, reference.alpha, RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(beta, reference.beta, RELATIVE_TOLERANCE);
    }
}

BOOST_AUTO_TEST_CASE(TestCarterConstantToAlphaBetaDomain) {
    const Schwarzschild<double> schwarzschild;
    const Kerr<double> kerr(0.5);
    double alpha = 0., beta = 0.;
    BOOST_CHECK_EQUAL(schwarzschild.angularMomentumCarterConstantToAlphaBeta(1., 0., 1., 0., alpha, beta),
                      StatusCode::FAILURE);
    BOOST_CHECK_EQUAL(kerr.angularMomentumCarterConstantToAlphaBeta(1., 0., 1., 0., alpha, beta), StatusCode::FAILURE);
    BOOST_CHECK_EQUAL(schwarzschild.angularMomentumCarterConstantToAlphaBeta(1., 0., 1., 1., alpha, beta),
                      StatusCode::DOMAIN_ERROR);
    BOOST_CHECK_EQUAL(kerr.angularMomentumCarterConstantToAlphaBeta(1., 0., 1., 1., alpha, beta),
                      StatusCode::DOMAIN_ERROR);
    BOOST_CHECK_EQUAL(schwarzschild.angularMomentumCarterConstantToAlphaBeta(0., 0., 0., 1., alpha, beta),
                      StatusCode::SUCCESS);
    BOOST_CHECK_SMALL(beta, RELATIVE_TOLERANCE);
    BOOST_CHECK_EQUAL(kerr.angularMomentumCarterConstantToAlphaBeta(0., 0., 0., 1., alpha, beta),
                      StatusCode::SUCCESS);
    BOOST_CHECK_SMALL(beta, RELATIVE_TOLERANCE);
}

BOOST_AUTO_TEST_CASE(TestNullNormalizationRequiresPositiveFrequency) {
    // [ORB-06] / [ORB-07] y[4] = dlambda/dt is the affine rate and must be positive
    const Kerr<double> kerr(0.5);
    const Schwarzschild<double> schwarzschild;
    boost::numeric::ublas::bounded_vector<double, 8> photon;
    photon.clear();
    photon[1] = 6.;
    photon[2] = 1.;
    photon[5] = 1.;
    photon[6] = 0.02;
    photon[7] = 0.03;
    for (double frequency : {0., -1.}) {
        auto kerr_photon = photon;
        auto schwarzschild_photon = photon;
        BOOST_CHECK_EQUAL(kerr.normalizeNullGeodesic(kerr_photon, frequency), StatusCode::DOMAIN_ERROR);
        BOOST_CHECK_EQUAL(schwarzschild.normalizeNullGeodesic(schwarzschild_photon, frequency), StatusCode::DOMAIN_ERROR);
    }
    auto kerr_photon = photon;
    auto schwarzschild_photon = photon;
    BOOST_CHECK_EQUAL(kerr.normalizeNullGeodesic(kerr_photon, 1.), StatusCode::SUCCESS);
    BOOST_CHECK_EQUAL(schwarzschild.normalizeNullGeodesic(schwarzschild_photon, 1.), StatusCode::SUCCESS);
}

BOOST_AUTO_TEST_CASE(TestKerrIsco) {
    for (const IscoReference &reference : ISCO_REFERENCE) {
        const Kerr<double> kerr(reference.spin);
        BOOST_CHECK_CLOSE_FRACTION(kerr.innermostStableCircularOrbitRadius(true), reference.prograde,
                                   RELATIVE_TOLERANCE);
        BOOST_CHECK_CLOSE_FRACTION(kerr.innermostStableCircularOrbitRadius(false), reference.retrograde,
                                   RELATIVE_TOLERANCE);
    }
}
