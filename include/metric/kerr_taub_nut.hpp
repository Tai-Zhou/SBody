/**
 * @file kerr_taub_nut.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-06
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include "infra/infra_context.hpp"
#include "metric/metric.hpp"

namespace sbody {
template <typename T>
class KerrTaubNUT : public Metric<T> {
  protected:
    const T a_, a2_, a4_;
    const T sqrt_delta_a2_;
    const T u_plus_1, u_plus, u_minus_1, u_minus, u_r;
    const T l_, l2_, l4_;

  public:
    KerrTaubNUT(T spin, T NUT)
        : a_(spin), a2_(a_ * a_), a4_(a2_ * a2_), sqrt_delta_a2_(std::sqrt(1. - a2_)), u_plus_1(1. + sqrt_delta_a2_),
          u_plus(1. / u_plus_1), u_minus_1(a2_ * u_plus), u_minus(1. / u_minus_1), u_r(-0.5 / sqrt_delta_a2_), l_(NUT),
          l2_(l_ * l_), l4_(l2_ * l2_) {}
    std::string name() const override { return "Kerr-Taub-NUT"; }
    int metricTensor(const boost::numeric::ublas::bounded_vector<T, 8> &position,
                     boost::numeric::ublas::bounded_matrix<T, 4, 4> &metric) const override {
        return StatusCode::NOT_IMPLEMENTED;
    }
    T dotProduct3D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x1, T x2, T x3, T y1, T y2,
                   T y3) const override {
        const T r = position[1], sin_theta = std::sin(position[2]), cos_theta = std::cos(position[2]);
        const T r2 = power2(r), sin2_theta = power2(sin_theta);
        const T Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + power2(l_ + a_ * cos_theta),
                chi = a_ * sin2_theta - 2. * l_ * cos_theta;
        const T rho_2 = 1. / rho2;
        return (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 +
               (power2(rho2 + a_ * chi) * sin2_theta - chi * chi * Delta) * rho_2 * x3 * y3;
    }
    T dotProduct4D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x0, T x1, T x2, T x3, T y0, T y1, T y2,
                   T y3) const override {
        const T r = position[1], sin_theta = std::sin(position[2]), cos_theta = std::cos(position[2]);
        const T r2 = power2(r), sin2_theta = power2(sin_theta);
        const T Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + power2(l_ + a_ * cos_theta),
                chi = a_ * sin2_theta - 2. * l_ * cos_theta;
        const T rho_2 = 1. / rho2;
        return (a2_ * sin2_theta - Delta) * rho_2 * x0 * y0 -
               2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * rho_2 * (x0 * y3 + x3 * y0) +
               (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 +
               (power2(rho2 + a_ * chi) * sin2_theta - chi * chi * Delta) * rho_2 * x3 * y3;
    }
    T distanceSquare(const boost::numeric::ublas::bounded_vector<T, 8> &x, const boost::numeric::ublas::bounded_vector<T, 8> &y,
                     bool include_time) const override {
        const T r = x[1], sin_theta = std::sin(x[2]), cos_theta = std::cos(x[2]), d0 = x[0] - y[0],
                d3 = phiDifference(x[3] - y[3]);
        const T r2 = power2(r), sin2_theta = power2(sin_theta);
        const T Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + power2(l_ + a_ * cos_theta),
                chi = a_ * sin2_theta - 2. * l_ * cos_theta;
        const T rho_2 = 1. / rho2;
        if (include_time)
            return (a2_ * sin2_theta - Delta) * rho_2 * power2(d0) -
                   4. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * rho_2 * d0 * d3 +
                   (rho2 / Delta) * power2(x[1] - y[1]) + rho2 * power2(x[2] - y[2]) +
                   (power2(rho2 + a_ * chi) * sin2_theta - power2(chi) * Delta) * rho_2 * power2(d3);
        return (rho2 / Delta) * power2(x[1] - y[1]) + rho2 * power2(x[2] - y[2]) +
               (power2(rho2 + a_ * chi) * sin2_theta - power2(chi) * Delta) * rho_2 * power2(d3);
    }
    int lagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T dt_dtau = 1. / y[4], r = y[1], r2 = power2(r);
        const T sin2_theta = power2(std::sin(y[2])), cos_theta = cosTheta(y[2]);
        const T Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + power2(l_ + a_ * cos_theta);
        if (rho2 == 0. || Delta == 0.)
            return StatusCode::DIVISION_BY_ZERO;
        const T rho_2 = 1. / rho2;
        y[4] = (y[4] * rho2 - Delta + a2_ * sin2_theta - 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]) *
               rho_2 * dt_dtau; // 1 + p_t
        y[5] *= rho2 / Delta * dt_dtau;
        y[6] *= rho2 * dt_dtau;
        y[7] = (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) +
                (power2(r2 + l2_ + a2_) * sin2_theta - power2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * y[7]) *
               rho_2 * dt_dtau;
        return StatusCode::SUCCESS;
    }
    int hamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T pt = y[4] - 1., r = y[1], r2 = power2(r);
        const T sin2_theta = power2(std::sin(y[2])), cos_theta = cosTheta(y[2]);
        const T Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + power2(l_ + a_ * cos_theta);
        if (rho2 == 0. || Delta == 0.)
            return StatusCode::DIVISION_BY_ZERO;
        const T rho_2 = 1. / rho2;
        y[4] = Delta * rho2 * sin2_theta /
               ((Delta * power2(a_ * sin2_theta - 2. * l_ * cos_theta) - power2(r2 + l2_ + a2_) * sin2_theta) * pt -
                2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]);
        y[5] *= Delta * rho_2 * y[4];
        y[6] *= rho_2 * y[4];
        y[7] = (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * pt + (Delta - a2_ * sin2_theta) * y[7]) /
               (Delta * rho2 * sin2_theta) * y[4];
        return StatusCode::SUCCESS;
    }
    int fastTrace(const T r_observer, const T theta_observer, const T sin_theta_observer, const T cos_theta_observer,
                  const T r_object, const T theta_object, const T phi_object, T &alpha, T &beta,
                  boost::numeric::ublas::bounded_vector<T, 8> &photon, T &photon_time, T absolute_accuracy, T relative_accuracy,
                  const InfraContext &infra_context) const override {
        return StatusCode::FAILURE;
    }
    T energy(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time, DynamicalSystem dynamics) const override {
        if (dynamics == DynamicalSystem::LAGRANGIAN) {
            const T r = y[1], r2 = power2(r);
            const T sin2_theta = power2(std::sin(y[2])), cos_theta = cosTheta(y[2]);
            const T Delta = r2 - 2. * r - l2_ + a2_;
            if (time == TimeSystem::T)
                return (Delta - a2_ * sin2_theta + 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]) /
                       ((r2 + power2(l_ + a_ * cos_theta)) * y[4]);
            // time == TimeSystem::TAU
            return ((Delta - a2_ * sin2_theta) * y[4] + 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]) /
                   (r2 + power2(l_ + a_ * cos_theta));
        } // dynamics == DynamicalSystem::HAMILTONIAN
        return 1. - y[4];
    }
    T angularMomentum(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time,
                      DynamicalSystem dynamics) const override {
        if (dynamics == DynamicalSystem::LAGRANGIAN) {
            const T r = y[1], r2 = power2(r);
            const T sin2_theta = power2(std::sin(y[2])), cos_theta = cosTheta(y[2]);
            const T Delta = r2 - 2. * r - l2_ + a2_;
            if (time == TimeSystem::T)
                return (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) +
                        (power2(r2 + l2_ + a2_) * sin2_theta - power2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * y[7]) /
                       ((r2 + power2(l_ + a_ * cos_theta)) * y[4]);
            // time == TimeSystem::TAU
            return (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[4] +
                    (power2(r2 + l2_ + a2_) * sin2_theta - power2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * y[7]) /
                   (r2 + power2(l_ + a_ * cos_theta));
        } // dynamics == DynamicalSystem::HAMILTONIAN
        return y[7];
    }
    T carterConstant(const boost::numeric::ublas::bounded_vector<T, 8> &y, const T mu2, TimeSystem time,
                     DynamicalSystem dynamics) const override {
        return std::numeric_limits<T>::quiet_NaN();
    }
    int normalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T r = y[1], r2 = power2(r);
        const T sin2_theta = power2(std::sin(y[2])), cos_theta = cosTheta(y[2]);
        const T Delta = r2 - 2. * r - l2_ + a2_;
        const T rho2 = r2 + power2(l_ + a_ * cos_theta);
        // y[7] += 2. * a_ * r / (power2(a2_ + r2) - a2_ * Delta * sin2_theta);
        y[4] = std::sqrt(
            ((Delta - a2_ * sin2_theta) + 4. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7] -
             (power2(r2 + l2_ + a2_) * sin2_theta - power2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * power2(y[7])) /
                rho2 -
            rho2 * (power2(y[5]) / Delta + power2(y[6])));
        return std::isnan(y[4]) ? StatusCode::DOMAIN_ERROR : StatusCode::SUCCESS;
    }
    int normalizeNullGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y, T frequency = 1.) const override {
        const T r = y[1], r2 = power2(r);
        const T sin2_theta = power2(std::sin(y[2])), cos_theta = cosTheta(y[2]);
        const T Delta = r2 - 2. * r - l2_ + a2_;
        const T l_a_cos_theta = l_ + a_ * cos_theta;
        const T rho2 = r2 + power2(l_a_cos_theta), rho_2 = 1. / rho2;
        const T chi = a_ * sin2_theta - 2. * l_ * cos_theta, rho2_a_chi = r2 + l2_ + a2_;
        const T a = rho2 / Delta * power2(y[5]) + rho2 * power2(y[6]) +
                    rho_2 * (power2(rho2_a_chi) * sin2_theta - power2(chi) * Delta) * power2(y[7]);
        const T b = -4. * rho_2 * ((r + l2_) * chi + l_ * cos_theta * rho2_a_chi) * y[7];
        const T c = -rho_2 * (Delta - a2_ * sin2_theta);
        const T coefficient = std::copysign(0.5 / a, frequency) * (-b + std::sqrt(b * b - 4. * a * c));
        if (std::isnan(coefficient))
            return StatusCode::DOMAIN_ERROR;
        y[4] = frequency;
        y[5] *= coefficient;
        y[6] *= coefficient;
        y[7] *= coefficient;
        return StatusCode::SUCCESS;
    }
    int getIntegrationSystem(std::function<void(const boost::numeric::ublas::bounded_vector<T, 8> &,
                                                boost::numeric::ublas::bounded_vector<T, 8> &, T)> &system,
                             TimeSystem time, DynamicalSystem dynamics,
                             MotionMode motion = MotionMode::GEODESIC) const override {
        if (motion != MotionMode::GEODESIC)
            throw std::logic_error("KerrNewman::getIntegrationSystem() invaild");
        if (time == TimeSystem::T) {
            if (dynamics == DynamicalSystem::LAGRANGIAN) {
                if (motion == MotionMode::GEODESIC) {
                    system = [a = this->a_, a2 = this->a2_, l = this->l_,
                              l2 = this->l2_](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                                              boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
                        dydt[0] = y[4]; // d\tau/dt
                        dydt[1] = y[5]; // dr/dt
                        dydt[2] = y[6]; // d\theta/dt
                        dydt[3] = y[7]; // d\phi/dt
                        const T r = y[1], r2 = power2(r);
                        const T sin_theta = sinTheta(y[2]), sin_1_theta = 1. / sin_theta, sin2_theta = power2(sin_theta),
                                cos_theta = cosTheta(y[2]), cos2_theta = power2(cos_theta);
                        const T Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
                        const T l_a_cos_theta = l + a * cos_theta, l_a_cos_theta2 = power2(l_a_cos_theta);
                        const T rho2 = r2 + l_a_cos_theta2, rho_2 = 1. / rho2, rho4 = power2(rho2), rho_4 = power2(rho_2),
                                rho_6 = rho_2 * rho_4;
                        const T chi = a * sin2_theta - 2. * l * cos_theta, rho2_a_chi = r2 + l2 + a2;
                        const T rho2_r_Delta = r * (r + 2. * l * l_a_cos_theta) - l_a_cos_theta2,
                                rho2_a_cos_theta = (2. * r + l * l_a_cos_theta) * l_a_cos_theta - r2 * l;
                        const T dydt4 =
                            2. * rho_4 *
                            (Delta_1 * rho2_a_chi * rho2_r_Delta * y[5] * (1. - chi * y[7]) -
                             sin_1_theta * chi * rho2_a_cos_theta * y[6] * (1. - chi * y[7]) -
                             2. * Delta_1 * rho2 * r * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * y[5] * y[7] -
                             sin_1_theta * rho4 * l * (1. + cos2_theta) * y[6] * y[7]);
                        // d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
                        dydt[4] = dydt4 * y[4];
                        // d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                        dydt[5] = -Delta * rho2_r_Delta * rho_6 * power2(1. - chi * y[7]) +
                                  (rho2_r_Delta - r * a2 * sin2_theta) * Delta_1 * rho_2 * power2(y[5]) +
                                  2. * a * sin_theta * l_a_cos_theta * rho_2 * y[5] * y[6] + r * Delta * rho_2 * power2(y[6]) +
                                  Delta * r * sin2_theta * rho_2 * power2(y[7]) + dydt4 * y[5];
                        // d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                        dydt[6] = rho2_a_cos_theta * sin_theta * rho_6 * (a - 2. * rho2_a_chi * y[7]) -
                                  a * sin_theta * l_a_cos_theta * rho_2 * (Delta_1 * power2(y[5]) - power2(y[6])) -
                                  2. * r * rho_2 * y[5] * y[6] -
                                  (rho2 * (Delta * chi * l_a_cos_theta - power2(rho2_a_chi) * cos_theta) -
                                   2. * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * rho2_a_chi * l_a_cos_theta) *
                                      sin_theta * rho_6 * power2(y[7]) +
                                  dydt4 * y[6];
                        // d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                        dydt[7] = -2. * a * rho2_r_Delta * Delta_1 * rho_4 * y[5] * (1. - chi * y[7]) +
                                  2. * rho2_a_cos_theta * rho_4 * sin_1_theta * y[6] * (1. - chi * y[7]) -
                                  2. * (1. - a2 * sin2_theta * Delta_1) * r * rho_2 * y[5] * y[7] -
                                  2. * cos_theta * sin_1_theta * y[6] * y[7] + dydt4 * y[7];
                    };
                    return StatusCode::SUCCESS;
                } else if (motion == MotionMode::CONICAL) {
                    system = [a = this->a_, a2 = this->a2_, l = this->l_,
                              l2 = this->l2_](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                                              boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) { // [TODO]: CHECK
                        dydt[0] = y[4];                                                                 // d\tau/dt
                        dydt[1] = y[5];                                                                 // dr/dt
                        dydt[2] = 0.;                                                                   // d\theta/dt = 0.
                        dydt[3] = y[7];                                                                 // d\phi/dt
                        const T r = y[1], r2 = power2(r);
                        const T sin2_theta = power2(std::sin(y[2])), cos_theta = cosTheta(y[2]);
                        const T Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
                        const T l_a_cos_theta2 = power2(l + a * cos_theta);
                        const T rho2 = r2 + l_a_cos_theta2, rho_2 = 1. / rho2;
                        const T chi = a * sin2_theta - 2. * l * cos_theta, rho2_a_chi = r2 + l2 + a2;
                        const T g00 = -rho_2 * (Delta - a2 * sin2_theta), g11 = rho2 * Delta_1,
                                g03 = -g00 * chi - a * sin2_theta, g33 = -g03 * chi + rho2_a_chi * sin2_theta;
                        const T L_1 = y[4] / (g03 + g33 * y[7]), L_2 = power2(L_1);
                        const T A = g33 * (g33 * L_2 + 1.), A_1 = 1. / A, B = 2. * g03 * (g33 * L_2 + 1.),
                                C = g00 + g11 * power2(y[5]) + power2(g03) * L_2;
                        const T dg00_dr = -2. * rho_2 * ((r - 1.) + g00 * r), dg11_dr = 2. * Delta_1 * (r - g11 * (r - 1.)),
                                dg03_dr = -dg00_dr * chi, dg33_dr = 2. * r * sin2_theta - dg03_dr * chi;
                        const T dA_dr = dg33_dr * (2. * g33 * L_2 + 1.),
                                dB_dr = 2. * (dg03_dr * (g33 * L_2 + 1.) + g03 * dg33_dr * L_2),
                                dC_dr = dg00_dr + dg11_dr * power2(y[5]) + 2. * g03 * dg03_dr * L_2;
                        // d^2r/dt^2 = 0.
                        dydt[5] = 0.;
                        // d^2\theta/dt^2 = 0.
                        dydt[6] = 0.;
                        // d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
                        dydt[7] =
                            y[5] * A_1 *
                            (-y[7] * dA_dr + 0.5 * (-dB_dr + (B * dB_dr - 2. * (A * dC_dr + C * dA_dr)) / (2. * A * y[7] + B)));
                        // d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt = d((g03 + g33 * d\phi/dt) / L)/dr * dr/dt
                        dydt[4] = (y[5] * (dg03_dr + dg33_dr * y[7]) + g33 * dydt[7]) * L_1;
                    };
                    return StatusCode::SUCCESS;
                }
            }
        } else if (time == TimeSystem::TAU) {
            if (dynamics == DynamicalSystem::LAGRANGIAN && motion == MotionMode::GEODESIC) {
                system = [a = this->a_, a2 = this->a2_, l = this->l_,
                          l2 = this->l2_](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                                          boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
                    dydt[0] = y[4]; // dt/d\tau
                    dydt[1] = y[5]; // dr/d\tau
                    dydt[2] = y[6]; // d\theta/d\tau
                    dydt[3] = y[7]; // d\phi/d\tau
                    const T r = y[1], r2 = power2(r);
                    const T sin_theta = sinTheta(y[2]), sin_1_theta = 1. / sin_theta, sin2_theta = power2(sin_theta),
                            cos_theta = cosTheta(y[2]), cos2_theta = power2(cos_theta);
                    const T Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
                    const T l_a_cos_theta = l + a * cos_theta, l_a_cos_theta2 = power2(l_a_cos_theta);
                    const T rho2 = r2 + l_a_cos_theta2, rho_2 = 1. / rho2, rho4 = power2(rho2), rho_4 = power2(rho_2),
                            rho_6 = rho_2 * rho_4;
                    const T chi = a * sin2_theta - 2. * l * cos_theta, rho2_a_chi = r2 + l2 + a2;
                    const T rho2_r_Delta = r * (r + 2. * l * l_a_cos_theta) - l_a_cos_theta2,
                            rho2_a_cos_theta = (2. * r + l * l_a_cos_theta) * l_a_cos_theta - r2 * l;
                    // d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
                    dydt[4] = -2. * rho_4 *
                              (Delta_1 * rho2_a_chi * rho2_r_Delta * y[5] * (y[4] - chi * y[7]) -
                               sin_1_theta * chi * rho2_a_cos_theta * y[6] * (y[4] - chi * y[7]) -
                               2. * Delta_1 * rho2 * r * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * y[5] * y[7] -
                               sin_1_theta * rho4 * l * (1. + cos2_theta) * y[6] * y[7]);
                    // d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                    dydt[5] = -Delta * rho2_r_Delta * rho_6 * power2(y[4] - chi * y[7]) +
                              (rho2_r_Delta - r * a2 * sin2_theta) * Delta_1 * rho_2 * power2(y[5]) +
                              2. * a * sin_theta * l_a_cos_theta * rho_2 * y[5] * y[6] + r * Delta * rho_2 * power2(y[6]) +
                              Delta * r * sin2_theta * rho_2 * power2(y[7]);
                    // d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                    dydt[6] = rho2_a_cos_theta * sin_theta * rho_6 * y[4] * (a * y[4] - 2. * rho2_a_chi * y[7]) -
                              a * sin_theta * l_a_cos_theta * rho_2 * (Delta_1 * power2(y[5]) - power2(y[6])) -
                              2. * r * rho_2 * y[5] * y[6] -
                              (rho2 * (Delta * chi * l_a_cos_theta - power2(rho2_a_chi) * cos_theta) -
                               2. * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * rho2_a_chi * l_a_cos_theta) *
                                  sin_theta * rho_6 * power2(y[7]);
                    // d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
                    dydt[7] = -2. * a * rho2_r_Delta * Delta_1 * rho_4 * y[5] * (y[4] - chi * y[7]) +
                              2. * rho2_a_cos_theta * rho_4 * sin_1_theta * y[6] * (y[4] - chi * y[7]) -
                              2. * (1. - a2 * sin2_theta * Delta_1) * r * rho_2 * y[5] * y[7] -
                              2. * cos_theta * sin_1_theta * y[6] * y[7];
                };
                return StatusCode::SUCCESS;
            }
        }
        throw std::logic_error("KerrNewman::getIntegrationSystem() invaild");
    }
};
} // namespace sbody
