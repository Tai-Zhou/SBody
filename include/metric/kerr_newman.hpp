/**
 * @file kerr_newman.hpp
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
class KerrNewman : public Metric<T> {
  protected:
    const T a_, a2_, a4_;
    const T sqrt_delta_a2_;
    const T u_plus_1, u_plus, u_minus_1, u_minus, u_r;
    const T r_Q_, r_Q2_, r_Q4_;

  public:
    KerrNewman(T spin, T charge)
        : a_(spin), a2_(a_ * a_), a4_(a2_ * a2_), sqrt_delta_a2_(std::sqrt(1. - a2_)), u_plus_1(1. + sqrt_delta_a2_),
          u_plus(1. / u_plus_1), u_minus_1(a2_ * u_plus), u_minus(1. / u_minus_1), u_r(-0.5 / sqrt_delta_a2_),
          r_Q_(0.5 * std::sqrt(M_1_PI) * charge), r_Q2_(r_Q_ * r_Q_), r_Q4_(r_Q2_ * r_Q2_) {}
    std::string name() const override { return "Kerr-Newman"; }
    int metricTensor(const boost::numeric::ublas::bounded_vector<T, 8> &position,
                     boost::numeric::ublas::bounded_matrix<T, 4, 4> &metric) const override {
        const T r = position[1], r2 = r * r, Delta = r2 - 2. * r + a2_ + r_Q2_, rho2 = r2 + a2_ * power2(std::cos(position[2])),
                sin2_theta = power2(std::sin(position[2]));
        if (rho2 == 0. || Delta == 0.)
            return StatusCode::DIVISION_BY_ZERO;
        const T r_rho_2 = (2. * r - r_Q2_) / rho2;
        metric.clear();
        metric(0, 0) = -(1. - r_rho_2);
        metric(0, 3) = -r_rho_2 * a_ * sin2_theta;
        metric(1, 1) = rho2 / Delta;
        metric(2, 2) = rho2;
        metric(3, 0) = -r_rho_2 * a_ * sin2_theta;
        metric(3, 3) = (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * sin2_theta;
        return StatusCode::SUCCESS;
    }
    T dotProduct3D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x1, T x2, T x3, T y1, T y2,
                   T y3) const override {
        const T r = position[1], r2 = power2(r), Delta = r2 - 2. * r + a2_ + r_Q2_,
                rho2 = r2 + a2_ * power2(std::cos(position[2])), r_rho_2 = (2. * r - r_Q2_) / rho2,
                sin2_theta = power2(std::sin(position[2]));
        return (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 +
               ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * sin2_theta * sin2_theta) * x3 * y3;
    }
    T dotProduct4D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x0, T x1, T x2, T x3, T y0, T y1, T y2,
                   T y3) const override {
        const T r = position[1], r2 = power2(r), Delta = r2 - 2. * r + a2_ + r_Q2_,
                rho2 = r2 + a2_ * power2(std::cos(position[2])), r_rho_2 = (2. * r - r_Q2_) / rho2,
                sin2_theta = power2(std::sin(position[2]));
        return (r_rho_2 - 1.) * x0 * y0 - r_rho_2 * a_ * sin2_theta * (x0 * y3 + x3 * y0) + (rho2 / Delta) * x1 * y1 +
               rho2 * x2 * y2 + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * sin2_theta * sin2_theta) * x3 * y3;
    }
    T distanceSquare(const boost::numeric::ublas::bounded_vector<T, 8> &x, const boost::numeric::ublas::bounded_vector<T, 8> &y,
                     bool include_time) const override {
        const T r = x[1], r2 = power2(r), d0 = x[0] - y[0], d3 = phiDifference(x[3] - y[3]);
        const T Delta = r2 - 2. * r + a2_ + r_Q2_, rho2 = r2 + a2_ * power2(std::cos(x[2])), r_rho_2 = (2. * r - r_Q2_) / rho2,
                sin2_theta = power2(std::sin(x[2]));
        if (include_time)
            return (r_rho_2 - 1.) * power2(d0) - 2. * r_rho_2 * a_ * sin2_theta * d0 * d3 +
                   (rho2 / Delta) * power2(x[1] - y[1]) + rho2 * power2(x[2] - y[2]) +
                   ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * power2(sin2_theta)) * power2(d3);
        return (rho2 / Delta) * power2(x[1] - y[1]) + rho2 * power2(x[2] - y[2]) +
               ((r2 + a2_) * sin2_theta + 2. * r * a2_ * power2(sin2_theta) / rho2) * power2(d3);
    }
    int lagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T dt_dtau = 1. / y[4], r2 = power2(y[1]), sin2_theta = power2(std::sin(y[2]));
        const T Delta = r2 - 2. * y[1] + a2_ + r_Q2_, rho2 = r2 + a2_ * power2(std::cos(y[2]));
        if (rho2 == 0. || Delta == 0.)
            return StatusCode::DIVISION_BY_ZERO;
        const T r_rho_2 = (2. * y[1] - r_Q2_) / rho2;
        y[4] = (y[4] - 1. + r_rho_2 * (1. - a_ * sin2_theta * y[7])) * dt_dtau; // 1 + p_t
        y[5] *= rho2 / Delta * dt_dtau;
        y[6] *= rho2 * dt_dtau;
        y[7] = (-r_rho_2 * a_ + (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * y[7]) * sin2_theta * dt_dtau;
        return StatusCode::SUCCESS;
    }
    int hamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T pt = y[4] - 1., r2 = power2(y[1]), sin2_theta = power2(std::sin(y[2]));
        const T Delta = r2 - 2. * y[1] + a2_ + r_Q2_, rho2 = r2 + a2_ * power2(std::cos(y[2]));
        if (rho2 == 0. || Delta == 0.)
            return StatusCode::DIVISION_BY_ZERO;
        const T rho_2 = 1. / rho2, r_rho_2 = (2. * y[1] - r_Q2_) / rho2;
        y[4] = -Delta / ((Delta + r_rho_2 * (a2_ + r2)) * pt + r_rho_2 * a_ * y[7]);
        y[5] *= Delta * y[4] * rho_2;
        y[6] *= y[4] * rho_2;
        if (y[7] == 0.)
            y[7] = (-r_rho_2 * a_ * pt + (1. - r_rho_2) / sin2_theta * y[7]) / Delta * y[4];
        else if (sin2_theta == 0.)
            return StatusCode::DIVISION_BY_ZERO;
        else
            y[7] = (-r_rho_2 * a_ * pt + (1. - r_rho_2) / sin2_theta * y[7]) / Delta * y[4];
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
            if (time == TimeSystem::T)
                return (1. - (2. * y[1] - r_Q2_) / (power2(y[1]) + a2_ * power2(std::cos(y[2]))) *
                                 (1. - a_ * power2(std::sin(y[2])) * y[7])) /
                       y[4];
            // time == TimeSystem::TAU
            return y[4] - (2. * y[1] - r_Q2_) / (power2(y[1]) + a2_ * power2(std::cos(y[2]))) *
                              (y[4] - a_ * power2(std::sin(y[2])) * y[7]);
        } // dynamics == DynamicalSystem::HAMILTONIAN
        return 1. - y[4];
    }
    T angularMomentum(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time,
                      DynamicalSystem dynamics) const override {
        if (dynamics == DynamicalSystem::LAGRANGIAN) {
            const T r2 = power2(y[1]), sin2_theta = power2(std::sin(y[2]));
            const T r_rho_2 = (2. * y[1] - r_Q2_) / (r2 + a2_ * power2(std::cos(y[2])));
            if (time == TimeSystem::T)
                return (-r_rho_2 * a_ + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta / y[4];
            // time == TimeSystem::TAU
            return (-r_rho_2 * a_ * y[4] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta;
        } // dynamics == DynamicalSystem::HAMILTONIAN
        return y[7];
    }
    T carterConstant(const boost::numeric::ublas::bounded_vector<T, 8> &y, const T mu2, TimeSystem time,
                     DynamicalSystem dynamics) const override {
        return std::numeric_limits<T>::quiet_NaN();
    }
    int normalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T r = y[1], r2 = power2(r);
        const T sin2_theta = power2(std::sin(y[2]));
        const T rho2 = r2 + a2_ * power2(std::cos(y[2]));
        const T r_rho_2 = (2. * r - r_Q2_) / rho2;
        y[4] = std::sqrt(1. - r_rho_2 + 2. * r_rho_2 * a_ * sin2_theta * y[7] -
                         (rho2 / (r2 - 2. * r + a2_ + r_Q2_) * power2(y[5]) + rho2 * power2(y[6]) +
                          ((a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * sin2_theta) * power2(y[7])));
        return std::isnan(y[4]) ? StatusCode::DOMAIN_ERROR : StatusCode::SUCCESS;
    }
    int normalizeNullGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y, T frequency = 1.) const override {
        const T r = y[1], r2 = power2(r);
        const T sin2_theta = power2(std::sin(y[2]));
        const T rho2 = r2 + a2_ * power2(std::cos(y[2]));
        const T r_rho_2 = (2. * r - r_Q2_) / rho2;
        const T a = rho2 / (r2 - 2. * r + a2_ + r_Q2_) * power2(y[5]) + rho2 * power2(y[6]) +
                    ((a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * sin2_theta) * power2(y[7]);
        const T b = -2. * r_rho_2 * a_ * sin2_theta * y[7];
        const T c = r_rho_2 - 1.;
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
        throw std::logic_error("KerrNewman::getIntegrationSystem() invaild");
    }
};
} // namespace sbody
