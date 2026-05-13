/**
 * @file reissner_nordstrom.hpp
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
class ReissnerNordstrom : public Metric<T> {
  protected:
    const T r_Q_, r_Q2_, r_Q4_;

  public:
    ReissnerNordstrom(T charge) : r_Q_(0.5 * std::sqrt(M_1_PI) * charge), r_Q2_(r_Q_ * r_Q_), r_Q4_(r_Q2_ * r_Q2_) {}
    std::string name() const override { return "Reissner-Nordstrom"; }
    int metricTensor(const boost::numeric::ublas::bounded_vector<T, 8> &position,
                     boost::numeric::ublas::bounded_matrix<T, 4, 4> &metric) const override {
        const T r_1 = 1. / position[1];
        metric.clear();
        metric(0, 0) = -1. + (2. - r_Q2_ * r_1) * r_1;
        metric(1, 1) = -1. / metric(0, 0);
        metric(2, 2) = position[1] * position[1];
        metric(3, 3) = power2(position[1] * std::sin(position[2]));
        return position[1] == 2. ? StatusCode::DIVISION_BY_ZERO : StatusCode::SUCCESS;
    }
    T dotProduct3D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x1, T x2, T x3, T y1, T y2,
                   T y3) const override {
        const T r_1 = 1. / position[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
        return x1 * y1 / g11_1 + position[1] * position[1] * x2 * y2 + power2(position[1] * std::sin(position[2])) * x3 * y3;
    }
    T dotProduct4D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x0, T x1, T x2, T x3, T y0, T y1, T y2,
                   T y3) const override {
        const T r_1 = 1. / position[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
        return -g11_1 * x0 * y0 + x1 * y1 / g11_1 + position[1] * position[1] * x2 * y2 +
               power2(position[1] * std::sin(position[2])) * x3 * y3;
    }
    T distanceSquare(const boost::numeric::ublas::bounded_vector<T, 8> &x, const boost::numeric::ublas::bounded_vector<T, 8> &y,
                     bool include_time) const override {
        const T r_1 = 1. / x[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
        if (include_time)
            return -g11_1 * power2(x[0] - y[0]) + power2(x[1] - y[1]) / g11_1 +
                   power2(x[1]) * (power2(x[2] - y[2]) + power2(std::sin(x[2]) * phiDifference(x[3] - y[3])));
        return power2(x[1] - y[1]) / g11_1 +
               power2(x[1]) * (power2(x[2] - y[2]) + power2(std::sin(x[2]) * phiDifference(x[3] - y[3])));
    }
    int lagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T dt_dtau = 1. / y[4], r_1 = 1. / y[1], r2 = power2(y[1]);
        const T rs_rQ2 = (2. - r_Q2_ * r_1) * r_1;
        y[4] = (y[4] - 1. + rs_rQ2) * dt_dtau; // 1 + p_t
        y[5] *= dt_dtau / (1. - rs_rQ2);
        y[6] *= r2 * dt_dtau;
        y[7] *= r2 * power2(std::sin(y[2])) * dt_dtau;
        return StatusCode::SUCCESS;
    }
    int hamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T r_1 = 1. / y[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
        y[4] = -g11_1 / (y[4] - 1.);
        y[5] *= g11_1 * y[4];
        y[6] *= power2(r_1) * y[4];
        y[7] *= power2(r_1 / std::sin(y[2])) * y[4];
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
            const T r_1 = 1. / y[1];
            if (time == TimeSystem::T)
                return (1. - (2. - r_Q2_ * r_1) * r_1) / y[4];
            // time == TimeSystem::TAU
            return (1. - (2. - r_Q2_ * r_1) * r_1) * y[4];
        } // dynamics == DynamicalSystem::HAMILTONIAN
        return 1. - y[4];
    }
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
    T carterConstant(const boost::numeric::ublas::bounded_vector<T, 8> &y, const T mu2, TimeSystem time,
                     DynamicalSystem dynamics) const override {
        if (dynamics == DynamicalSystem::LAGRANGIAN) {
            if (time == TimeSystem::T)
                return power4(y[1]) * (power2(y[6]) + power2(y[7] * std::cos(y[2]) * std::sin(y[2]))) / power2(y[4]);
            // time == TimeSystem::TAU
            return power4(y[1]) * (power2(y[6]) + power2(y[7] * std::cos(y[2]) * std::sin(y[2])));
        } // dynamics == DynamicalSystem::HAMILTONIAN
        return power2(y[6]) + power2(y[7] / tan(y[2]));
    }
    int normalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
        const T r_1 = 1. / y[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
        if (g11_1 <= 0)
            return 1;
        y[4] = std::sqrt(g11_1 - (power2(y[5]) / g11_1 + power2(y[1] * y[6]) + power2(y[1] * std::sin(y[2]) * y[7])));
        return std::isnan(y[4]) ? StatusCode::DOMAIN_ERROR : StatusCode::SUCCESS;
    }
    int normalizeNullGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y, T frequency = 1.) const override {
        const T r_1 = 1. / y[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
        if (g11_1 <= 0)
            return 1;
        const T coefficient = std::copysign(g11_1, frequency) /
                              std::sqrt(power2(y[5]) + g11_1 * (power2(y[1] * y[6]) + power2(y[1] * std::sin(y[2]) * y[7])));
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
        if (time != TimeSystem::T || dynamics != DynamicalSystem::LAGRANGIAN || motion != MotionMode::GEODESIC)
            throw std::logic_error("ReissnerNordstrom::getIntegrationSystem() invaild");
        system = [r_Q2 = this->r_Q2_](const boost::numeric::ublas::bounded_vector<T, 8> &y,
                                      boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
            dydt[0] = y[4]; // d\tau/dt
            dydt[1] = y[5]; // dr/dt
            dydt[2] = y[6]; // d\theta/dt
            dydt[3] = y[7]; // d\phi/dt
            const T r = y[1], r_1 = 1. / r, r_2 = power2(r_1);
            const T sin_theta = sinTheta(y[2]), cos_theta = cosTheta(y[2]);
            const T g11_1 = 1. - (2. - r_Q2 * r_1) * r_1, g11 = 1. / g11_1;
            const T m_r_Q_r2 = (1. - r_Q2 * r_1) * r_2, m_r_Q_r2_g11 = m_r_Q_r2 * g11;
            dydt[4] = 2. * m_r_Q_r2_g11 * y[5] * y[4];
            dydt[5] =
                -g11_1 * m_r_Q_r2 + 3. * m_r_Q_r2_g11 * power2(y[5]) + r * g11_1 * (power2(y[6]) + power2(sin_theta * y[7]));
            dydt[6] = -2. * (g11_1 * r_1 - m_r_Q_r2) * g11 * y[5] * y[6] + sin_theta * cos_theta * power2(y[7]);
            if (y[7] == 0.)
                dydt[7] = 0.;
            else if (sin_theta == 0.)
                throw std::domain_error("Non-zero angular velocity with sin(theta) == 0.");
            else
                dydt[7] = -2. * ((r_1 - m_r_Q_r2_g11) * y[5] + cos_theta / sin_theta * y[6]) * y[7];
        };
        return StatusCode::SUCCESS;
    }
};
} // namespace sbody
