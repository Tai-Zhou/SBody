/**
 * @file post_newtonian.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-06
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef SBODY_POST_NEWTONIAN_HPP
#define SBODY_POST_NEWTONIAN_HPP

#include "metric/metric.hpp"

namespace sbody {
	/// Post-Newtonian
	template <typename T>
	class Newton : public Metric<T> {
	  protected:
		const int PN_;

	  public:
		Newton(int PN) : PN_(PN) {};
		std::string Name() const override {
			return "Newton";
		}
		int MetricTensor(const boost::numeric::ublas::bounded_vector<T, 8> &position, boost::numeric::ublas::bounded_matrix<T, 4, 4> &metric) const override { // FIXME: Check
			metric.clear();
			metric(0, 0) = -1.;
			metric(1, 1) = 1.;
			metric(2, 2) = position(1) * position(1);
			metric(3, 3) = power2(position[1] * std::sin(position[2]));
			return StatusCode::SUCCESS;
		}
		T DotProduct3D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x1, T x2, T x3, T y1, T y2, T y3) const override {
			return x1 * y1 + power2(position[1]) * x2 * y2 + power2(position[1] * std::sin(position[2])) * x3 * y3;
		}
		T DotProduct4D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x0, T x1, T x2, T x3, T y0, T y1, T y2, T y3) const override {
			return -x0 * y0 + x1 * y1 + power2(position[1]) * x2 * y2 + power2(position[1] * std::sin(position[2])) * x3 * y3;
		}
		T DistanceSquare(const boost::numeric::ublas::bounded_vector<T, 8> &x, const boost::numeric::ublas::bounded_vector<T, 8> &y, bool include_time) const override {
			if (include_time)
				return -power2(x[0] - y[0]) + power2(x[1] - y[1]) + power2(x[1] * (x[2] - y[2])) + power2(x[1] * std::sin(x[2]) * (x[3] - y[3]));
			return power2(x[1] - y[1]) + power2(x[1] * (x[2] - y[2])) + power2(x[1] * std::sin(x[2]) * (x[3] - y[3]));
		}
		int LagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
			return StatusCode::FAILURE;
		}
		int HamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
			return StatusCode::FAILURE;
		}
		int FastTrace(const T r_observer, const T theta_observer, const T sin_theta_observer, const T cos_theta_observer, const T r_object, const T theta_object, const T phi_object, T &alpha, T &beta, boost::numeric::ublas::bounded_vector<T, 8> &photon, T &photon_time) const override {
			const T sin_theta_object = sinTheta(theta_object);
			photon[1] = r_object * sin_theta_object * std::cos(phi_object);
			photon[2] = r_object * sin_theta_object * std::sin(phi_object);
			photon[3] = r_object * cosTheta(theta_object);
			const T dx = r_observer * sin_theta_observer - photon[1], dz = r_observer * cos_theta_observer - photon[3];
			photon[0] = photon_time = std::hypot(dx, photon[2], dz);
			if (photon[0] == 0.)
				return StatusCode::DIVISION_BY_ZERO;
			alpha = photon[2];
			beta = photon[3] * sin_theta_observer - photon[1] * cos_theta_observer;
			photon[4] = 1.;
			const T distance_1 = 1. / photon[0];
			photon[5] = dx * distance_1;
			photon[6] = photon[2] * distance_1;
			photon[7] = dz * distance_1;
			cartesianToSpherical(photon);
			return StatusCode::SUCCESS;
		}
		T Energy(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time, DynamicalSystem dynamics) const override {
			const T r_1 = 1. / y[1], r_2 = power2(r_1), r_3 = r_1 * r_2, r_4 = power2(r_2);
			const T v2 = power2(y[5]) + power2(y[1]) * (power2(y[6]) + power2(std::sin(y[2]) * y[7])), v4 = power2(v2), v6 = v2 * v4, v8 = power2(v4);
			const T rdot = y[5], rdot2 = power2(rdot);
			T E = 0.5 * v2 - r_1;
			if (PN_ & 1)
				E += 0.5 * r_2 + 0.375 * v4 + 1.5 * r_1 * v2;
			if (PN_ & 2)
				E += -0.5 * r_3 + 0.3125 * v6 + 1.75 * r_2 * v2 + 0.5 * r_2 * rdot2 + 2.625 * r_1 * v4;
			if (PN_ & 8)
				E += 0.375 * r_4 + 1.25 * r_3 * v2 + 1.5 * r_3 * rdot2 + 0.2734375 * v8 + 8.4375 * r_2 * v4 + 0.75 * r_2 * v2 * rdot2 + 3.4375 * r_1 * v6;
			return E;
		}
		T AngularMomentum(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time, DynamicalSystem dynamics) const override {
			const T r_1 = 1. / y[1], r_2 = power2(r_1), r_3 = power3(r_1);
			const T v_tan2 = power2(y[1]) * (power2(y[6]) + power2(std::sin(y[2]) * y[7]));
			const T v2 = power2(y[5]) + v_tan2, v4 = power2(v2), v6 = v2 * v4;
			const T rdot = y[5], rdot2 = power2(rdot);
			T eff = 1.;
			if (PN_ & 1)
				eff += 3. * r_1 + 0.5 * v2;
			if (PN_ & 2)
				eff += 3.5 * r_2 + 0.375 * v4 + 3.5 * r_1 * v2;
			if (PN_ & 8)
				eff += 2.5 * r_3 + 0.3125 * v6 + 11.25 * r_2 * v2 + 0.5 * r_2 * rdot2 + 4.125 * r_1 * v4;
			return y[1] * std::sqrt(v_tan2) * eff;
		}
		T CarterConstant(const boost::numeric::ublas::bounded_vector<T, 8> &y, const T mu2, TimeSystem time, DynamicalSystem dynamics) const override {
			return power2(AngularMomentum(y, time, dynamics));
		}
		T Redshift(const boost::numeric::ublas::bounded_vector<T, 8> &y, const boost::numeric::ublas::bounded_vector<T, 8> &photon_position, TimeSystem object_time, TimeSystem photon_time) const override {
			namespace ublas = boost::numeric::ublas;
			const T delta_epsilon = 1. - 2. / y[1];
			const ublas::vector_range<const ublas::bounded_vector<T, 8>> velocity(y, ublas::range(4, 8));
			return (1. - this->DotProduct(y, velocity, ublas::vector_range<const ublas::bounded_vector<T, 8>>(photon_position, ublas::range(4, 8)), false) / std::sqrt(delta_epsilon)) / std::sqrt(delta_epsilon - this->DotProduct(y, velocity, velocity, false));
		}
		int NormalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
			y[4] = 1;
			if (y[5] * y[5] + y[1] * y[1] * (y[6] * y[6] + power2(std::sin(y[2]) * y[7])) >= 1.)
				return StatusCode::FAILURE;
			return StatusCode::SUCCESS;
		}
		int NormalizeNullGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y, T frequency = 1.) const override {
			y[4] = 1;
			const T v_1 = 1. / std::sqrt(power2(y[5]) + power2(y[1]) * (power2(y[6]) + power2(std::sin(y[2]) * y[7])));
			for (int i = 5; i < 8; ++i)
				y[i] *= v_1;
			return StatusCode::SUCCESS;
		}
		int GetIntegrationSystem(std::function<void(const boost::numeric::ublas::bounded_vector<T, 8> &, boost::numeric::ublas::bounded_vector<T, 8> &, T)> &system, TimeSystem time, DynamicalSystem dynamics, MotionMode motion = MotionMode::GEODESIC) const override {
			if (time != TimeSystem::T || dynamics != DynamicalSystem::LAGRANGIAN || motion != MotionMode::GEODESIC) {
				throw std::logic_error("Newton::GetIntegrationSystem() invaild");
			}
			system = [PN = this->PN_](const boost::numeric::ublas::bounded_vector<T, 8> &y, boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
				const T r_1 = 1. / y[1], r_2 = power2(r_1), sin_theta = sinTheta(y[2]), cos_theta = cosTheta(y[2]);
				const T v_tan2 = power2(y[1]) * (power2(y[6]) + power2(sin_theta * y[7]));
				const T rdot = y[5], rdot2 = power2(y[5]);
				const T v2 = rdot2 + v_tan2;
				dydt[0] = y[4];
				dydt[1] = y[5];
				dydt[2] = y[6];
				dydt[3] = y[7];
				dydt[4] = 0.;
				T A = 1., B = 0.;
				if (PN & 1) {
					A += v2 - 4. * r_1;
					B += -4. * rdot;
				}
				if (PN & 2) {
					A += r_1 * (-2. * rdot2 + 9. * r_1);
					B += 2. * r_1 * rdot;
				}
				if (PN & 8) {
					A += -16. * power3(r_1) + rdot2 * r_2;
					B += -4. * rdot * r_2;
				}
				const T B_r2 = B * r_2;
				dydt[5] = r_1 * v_tan2 - r_2 * A - B_r2 * y[5];
				dydt[6] = sin_theta * cos_theta * power2(y[7]) - (2. * y[5] * r_1 + B_r2) * y[6];
				dydt[7] = -(2. * cos_theta / sin_theta * y[6] + 2. * y[5] * r_1 + B_r2) * y[7];
			};
			return StatusCode::SUCCESS;
		}
	};
	template <typename T>
	class PN1 : public Newton<T> {
	  protected:
		const T PN1_;

	  public:
		PN1(T fSP) : Newton<T>(1), PN1_(fSP) {};
		std::string Name() const override {
			return "Newton-PN1";
		}
		int GetIntegrationSystem(std::function<void(const boost::numeric::ublas::bounded_vector<T, 8> &, boost::numeric::ublas::bounded_vector<T, 8> &, T)> &system, TimeSystem time, DynamicalSystem dynamics, MotionMode motion = MotionMode::GEODESIC) const override {
			if (time != TimeSystem::T || dynamics != DynamicalSystem::LAGRANGIAN || motion != MotionMode::GEODESIC) {
				throw std::logic_error("PN1::GetIntegrationSystem() invaild");
			}
			system = [PN1 = this->PN1_](const boost::numeric::ublas::bounded_vector<T, 8> &y, boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
				const T r_1 = 1. / y[1], r_2 = r_1 * r_1, sin_theta = sinTheta(y[2]), cos_theta = cosTheta(y[2]);
				const T v_tan2 = y[1] * y[1] * (y[6] * y[6] + power2(sin_theta * y[7]));
				const T v2 = y[5] * y[5] + v_tan2;
				dydt[0] = y[4];
				dydt[1] = y[5];
				dydt[2] = y[6];
				dydt[3] = y[7];
				dydt[4] = 0.;
				const T A = 1. + (v2 - 4. * r_1) * PN1, B_r2 = -4. * y[5] * PN1 * r_2;
				dydt[5] = r_1 * v_tan2 - r_2 * A - B_r2 * y[5];
				dydt[6] = sin_theta * cos_theta * y[7] * y[7] - (2. * y[5] * r_1 + B_r2) * y[6];
				dydt[7] = -(2. * cos_theta / sin_theta * y[6] + 2. * y[5] * r_1 + B_r2) * y[7];
			};
			return StatusCode::SUCCESS;
		}
	};
} // namespace sbody

#endif
