/**
 * @file hayward.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-06
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include "metric/metric.hpp"

namespace sbody {
template <typename T>
class Hayward : public Metric<T> {
  protected:
    const T a_, a2_, a4_;
    const T sqrt_delta_a2_;
    const T u_plus_1, u_plus, u_minus_1, u_minus, u_r;
    const T alpha_, beta_, g_;

  public:
    Hayward(T alpha, T beta, T charge);
    std::string name() const override;
    int metricTensor(const boost::numeric::ublas::bounded_vector<T, 8> &position,
                     boost::numeric::ublas::bounded_matrix<T, 4, 4> &metric) const override;
    T dotProduct3D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x1, T x2, T x3, T y1, T y2,
                   T y3) const override;
    T dotProduct4D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x0, T x1, T x2, T x3, T y0, T y1, T y2,
                   T y3) const override;
    T distanceSquare(const boost::numeric::ublas::bounded_vector<T, 8> &x, const boost::numeric::ublas::bounded_vector<T, 8> &y,
                     bool include_time) const override;
    int lagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override;
    int hamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override;
    int fastTrace(const T r_observer, const T theta_observer, const T sin_theta_observer, const T cos_theta_observer,
                  const T r_object, const T theta_object, const T phi_object, T &alpha, T &beta,
                  boost::numeric::ublas::bounded_vector<T, 8> &photon, T &photon_time) const override;
    T energy(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time, DynamicalSystem dynamics) const override;
    T angularMomentum(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time,
                      DynamicalSystem dynamics) const override;
    T carterConstant(const boost::numeric::ublas::bounded_vector<T, 8> &y, const T mu2, TimeSystem time,
                     DynamicalSystem dynamics) const override;
    int normalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y) const override;
    int normalizeNullGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y, T frequency = 1.) const override;
    int getIntegrationSystem(std::function<void(const boost::numeric::ublas::bounded_vector<T, 8> &,
                                                boost::numeric::ublas::bounded_vector<T, 8> &, T)> &system,
                             TimeSystem time, DynamicalSystem dynamics,
                             MotionMode motion = MotionMode::GEODESIC) const override;
};
} // namespace sbody
