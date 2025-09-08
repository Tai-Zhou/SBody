/**
 * @file Metric.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_METRIC_HPP
#define SBODY_METRIC_HPP

#include <cerrno>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <tuple>

// #include <boost/math/tools/quartic_roots.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "IO.hpp"
#include "Utility.hpp"

namespace SBody {
	/**
	 * @brief The time system used in the record of the object.
	 *
	 * - T: the time step of the integrator is the observer's time \f$t\f$, and the first element of the object vector is its intrinsic time \f$\tau\f$.
	 * - TAU: the time step of the integrator is the intrinsic time of the object, and the first element is the observer's time.
	 *
	 */
	enum TimeSystem { T,
					  TAU };

	/**
	 * @brief The dynamical system to describe the motion of the object.
	 *
	 * - LAGRANGIAN: the velocity part of the object vector is \f$\mathrm d{x^\mu}/\mathrm dt\f$ or \f$\mathrm d{x^\mu}/\mathrm d\tau\f$.
	 * - HAMILTONIAN: the velocity is recorded as \f$p_\mu\f$.
	 *
	 */
	enum DynamicalSystem { LAGRANGIAN,
						   HAMILTONIAN };

	/**
	 * @brief The mode of the motion of the object.
	 *
	 * - GEODESIC: the motion is described by the geodesic equation.
	 * - CIRCULAR: the orbit has fixed radius.
	 * - HELICAL: the motion is described by the helical orbit equation.
	 *
	 */
	enum MotionMode { GEODESIC,
					  CIRCULAR,
					  HELICAL,
					  CYLINDRICAL };

	/**
	 * @brief The base class of all kinds of metrics.
	 *
	 */
	template <typename Type>
	class Metric {
	  protected:
		virtual Type DotProduct3D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x1, Type x2, Type x3, Type y1, Type y2, Type y3) = 0;
		virtual Type DotProduct4D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x0, Type x1, Type x2, Type x3, Type y0, Type y1, Type y2, Type y3) = 0;

	  public:
		/**
		 * @brief Return the name of the metric.
		 *
		 * @return name of the metric.
		 */
		virtual std::string Name() const = 0;

		/**
		 * @brief Calculate the metric tensor at `position`, stored in `metric`.
		 *
		 * @param position 4 dimensional vector
		 * @param metric matrix with size 4×4
		 * @return status
		 */
		virtual int MetricTensor(const boost::numeric::ublas::bounded_vector<Type, 8> &position, boost::numeric::ublas::bounded_matrix<Type, 4, 4> &metric) = 0;

		/**
		 * @brief Dot product of vector `x` and `y` at `position`. \f$g_{\mu\nu}x^\mu y^\nu\f$
		 *
		 * @param position 4 dimensional vector, position to calcuate the dot product of `x` and `y`.
		 * @param x 4 dimensional vector
		 * @param y 4 dimensional vector
		 * @param include_time include time related terms if true.
		 * @return result
		 */
		template <typename VecExp1, typename VecExp2>
		Type DotProduct(const boost::numeric::ublas::bounded_vector<Type, 8> &position, const VecExp1 &x, const VecExp2 &y, bool include_time) {
			if (include_time)
				return DotProduct4D(position, x(0), x(1), x(2), x(3), y(0), y(1), y(2), y(3));
			return DotProduct3D(position, x(1), x(2), x(3), y(1), y(2), y(3));
		}

		/**
		 * @brief Calculate the square of the distance between `x` and `y` at `x`. \f$g_{\mu\nu}(x^\mu-y^\mu)(x^\nu-y^\nu)\f$
		 *
		 * @param x 8 dimensional vector
		 * @param y 8 dimensional vector
		 * @return result
		 */
		virtual Type DistanceSquare(const boost::numeric::ublas::bounded_vector<Type, 8> &x, const boost::numeric::ublas::bounded_vector<Type, 8> &y, bool include_time) = 0;

		/**
		 * @brief Calculate the local inertial frame coordinate of the object at `position`, stored in `coordinate`.
		 *
		 * @param position 8 dimensional vector
		 * @param time time systme of the object
		 * @param coordinate matrix with size 4×4
		 * @return status
		 */
		int LocalInertialFrame(const boost::numeric::ublas::bounded_vector<Type, 8> &position, TimeSystem time, boost::numeric::ublas::bounded_matrix<Type, 4, 4> &coordinate) {
			namespace ublas = boost::numeric::ublas;
			coordinate = ublas::identity_matrix<Type>(4);
			ublas::bounded_matrix<Type, 4, 4> metric, product = ublas::identity_matrix<Type>(4), product_LU;
			MetricTensor(position, metric);
			ublas::permutation_matrix<std::size_t> permutation(4);
			for (int i = 0; i < 4; ++i) {
				auto coordinate_row = ublas::matrix_row<ublas::bounded_matrix<Type, 4, 4>>(coordinate, i);
				if (i == 0) {
					if (time == T)
						std::copy(position.begin() + 5, position.end(), coordinate_row.begin() + 1);
					else
						std::copy(position.begin() + 4, position.end(), coordinate_row.begin());
				} else {
					product_LU = product;
					if (int status = ublas::lu_factorize(product_LU, permutation); status != Status::SUCCESS)
						return status;
					ublas::lu_substitute(product_LU, permutation, coordinate_row);
				}
				coordinate_row /= std::sqrt(std::abs(DotProduct(position, coordinate_row, coordinate_row, true)));
				ublas::matrix_row<ublas::bounded_matrix<Type, 4, 4>>(product, i) = ublas::prod(metric, coordinate_row);
			}
			return isnan(coordinate(3, 0)) ? Status::NUMERIC_ERROR : Status::SUCCESS;
		}

		/**
		 * @brief Convert the coordinate system from Lagrangian to Hamiltonian.
		 *
		 * @param y 8 dimensional vector
		 * @return status
		 */
		virtual int LagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<Type, 8> &y) = 0;

		/**
		 * @brief Convert the coordinate system from Hamiltonian to Lagrangian.
		 *
		 * @param y 8 dimensional vector
		 * @return status
		 */
		virtual int HamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<Type, 8> &y) = 0;

		int InitializePhoton(boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time, Type alpha, Type beta, Type r, Type r2, Type theta, Type sin_theta) {
			photon[0] = 0.;
			photon[1] = r;
			photon[4] = 1.;
			photon[5] = 1.;
			photon_time = r;
			if (sin_theta < boost::math::tools::root_epsilon<Type>()) {
				const Type k = std::hypot(alpha, beta);
				if (theta < boost::math::constants::half_pi<Type>()) {
					photon[2] = 1e-15;
					photon[3] = atan2(alpha, -beta);
					photon[6] = -k / r2;
				} else {
					photon[2] = boost::math::constants::pi<Type>() - 1e-15;
					photon[3] = atan2(alpha, beta);
					photon[6] = k / r2;
				}
				photon[7] = 0.;
			} else {
				photon[2] = theta;
				photon[3] = 0.,
				photon[6] = beta / r2;
				photon[7] = -alpha / (r2 * sin_theta);
			}
			return NormalizeNullGeodesic(photon, 1.);
		}

		int AngularMomentumCarterConstantToAlphaBeta(Type l, Type q2, Type cos_theta, Type sin_theta, Type &alpha, Type &abs_beta) {
			if (sin_theta == 0.)
				return Status::FAILURE;
			alpha = -l / sin_theta;
			const Type beta2 = q2 - Power2(cos_theta * alpha);
			if (beta2 < 0.)
				return Status::DOMAIN_ERROR;
			abs_beta = std::sqrt(beta2);
			return Status::SUCCESS;
		}

		/**
		 * @brief Trace the photon from the observer to the target, using elliptic integrals.
		 *
		 * @param r_observer radius of the observer, \f$r_\text{obs}\f$.
		 * @param theta_observer theta of the observer, \f$\theta_\text{obs}\f$.
		 * @param sin_theta_observer std::sin(theta_observer), \f$\std::sin\theta_\text{obs}\f$.
		 * @param cos_theta_observer std::cos(theta_observer), \f$\std::cos\theta_\text{obs}\f$.
		 * @param r_object radius of the target, \f$r_\text{tar}\f$.
		 * @param theta_object theta of the target, \f$\theta_\text{tar}\f$.
		 * @param phi_object phi of the target, \f$\phi_\text{tar}\f$.
		 * @param alpha x position of the target in the observer's view.
		 * @param beta y position of the target in the observer's view.
		 * @param photon 8 dimensional vector, position and the velocity of the photon traced to the target.
		 * @param photon_time look back time of the photon.
		 * @return status
		 */
		virtual int FastTrace(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type r_object, const Type theta_object, const Type phi_object, Type &alpha, Type &beta, boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time) = 0;

		virtual int FastShadow(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type alpha, const Type beta, const Type r_min) {
			return Status::FAILURE;
		}

		/**
		 * @brief Calculate the energy of the object.
		 *
		 * @param y 8 dimensional vector
		 * @param time time system of the object
		 * @param dynamics dynamical system of the object
		 * @return result
		 */
		virtual Type Energy(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) = 0;
		/**
		 * @brief Calculate the angular momentum of the object.
		 *
		 * @param y 8 dimensional vector
		 * @param time time system of the object
		 * @param dynamics dynamical system of the object
		 * @return result
		 */
		virtual Type AngularMomentum(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) = 0;

		/**
		 * @brief Calculate the Carter constant of the object.
		 *
		 * @param y 8 dimensional vector
		 * @param time time system of the object
		 * @param dynamics dynamical system of the object
		 * @return result
		 */
		virtual Type CarterConstant(const boost::numeric::ublas::bounded_vector<Type, 8> &y, const Type mu2, TimeSystem time, DynamicalSystem dynamics) = 0;

		/**
		 * @brief Calculate the redshift of the object, \f$1+z\f$
		 *
		 * @param y 8 dimensional vector of the object
		 * @param photon photon traced to the object
		 * @param object_time time system of the object
		 * @param photon_time time system of the photon
		 * @return result
		 */
		virtual Type Redshift(const boost::numeric::ublas::bounded_vector<Type, 8> &y, const boost::numeric::ublas::bounded_vector<Type, 8> &photon_position, TimeSystem object_time, TimeSystem photon_time) {
			namespace ublas = boost::numeric::ublas;
			ublas::bounded_vector<Type, 4> u, v;
			u(0) = 1.;
			std::copy(y.begin() + 5, y.end(), u.begin() + 1);
			v(0) = 1.;
			std::copy(photon_position.begin() + 5, photon_position.end(), v.begin() + 1);
			if (object_time == T) {
				if (photon_time == T)
					return -DotProduct(y, u, v, true) / (y[4] * photon_position[4]);
				return -DotProduct(y, u, ublas::vector_range<const ublas::bounded_vector<Type, 8>>(photon_position, ublas::range(4, 8)), true) / y[4];
			}
			if (photon_time == T)
				return -DotProduct(y, ublas::vector_range<const ublas::bounded_vector<Type, 8>>(y, ublas::range(4, 8)), v, true) / photon_position[4];
			return -DotProduct(y, ublas::vector_range<const ublas::bounded_vector<Type, 8>>(y, ublas::range(4, 8)), ublas::vector_range<const ublas::bounded_vector<Type, 8>>(photon_position, ublas::range(4, 8)), true);
		}

		/**
		 * @brief Normalize the timelike geodesic.
		 *
		 * @param y 8 dimensional vector
		 * @return status
		 */
		virtual int NormalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y) = 0;

		/**
		 * @brief Normalize the null geodesic.
		 *
		 * @param y 8 dimensional vector
		 * @return status
		 */
		virtual int NormalizeNullGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y, Type frequency = 1.) = 0;

		/**
		 * @brief Get the integration system to calculate the motion of the object.
		 *
		 * @param time time system of the object
		 * @param dynamics dynamical system of the object
		 * @param motion motion mode of the object
		 * @return integrator system
		 */
		virtual std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> GetIntegrationSystem(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) = 0;
	};

	/// Post-Newtonian
	template <typename Type>
	class Newton : public Metric<Type> {
	  protected:
		const int PN_;
		Type DotProduct3D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x1, Type x2, Type x3, Type y1, Type y2, Type y3) override {
			return x1 * y1 + Power2(position[1]) * x2 * y2 + Power2(position[1] * std::sin(position[2])) * x3 * y3;
		}
		Type DotProduct4D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x0, Type x1, Type x2, Type x3, Type y0, Type y1, Type y2, Type y3) override {
			return -x0 * y0 + x1 * y1 + Power2(position[1]) * x2 * y2 + Power2(position[1] * std::sin(position[2])) * x3 * y3;
		}

	  public:
		Newton(int PN) : PN_(PN) {};
		std::string Name() const override {
			return "Newton";
		}
		int MetricTensor(const boost::numeric::ublas::bounded_vector<Type, 8> &position, boost::numeric::ublas::bounded_matrix<Type, 4, 4> &metric) override { // FIXME: Check
			metric.clear();
			metric(0, 0) = -1.;
			metric(1, 1) = 1.;
			metric(2, 2) = position(1) * position(1);
			metric(3, 3) = Power2(position[1] * std::sin(position[2]));
			return Status::SUCCESS;
		}
		Type DistanceSquare(const boost::numeric::ublas::bounded_vector<Type, 8> &x, const boost::numeric::ublas::bounded_vector<Type, 8> &y, bool include_time) override {
			if (include_time)
				return -Power2(x[0] - y[0]) + Power2(x[1] - y[1]) + Power2(x[1] * (x[2] - y[2])) + Power2(x[1] * std::sin(x[2]) * (x[3] - y[3]));
			return Power2(x[1] - y[1]) + Power2(x[1] * (x[2] - y[2])) + Power2(x[1] * std::sin(x[2]) * (x[3] - y[3]));
		}
		int LagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			return Status::FAILURE;
		}
		int HamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			return Status::FAILURE;
		}
		int FastTrace(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type r_object, const Type theta_object, const Type phi_object, Type &alpha, Type &beta, boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time) override {
			const Type sin_theta_object = SinTheta(theta_object);
			photon[1] = r_object * sin_theta_object * std::cos(phi_object);
			photon[2] = r_object * sin_theta_object * std::sin(phi_object);
			photon[3] = r_object * CosTheta(theta_object);
			const Type dx = r_observer * sin_theta_observer - photon[1], dz = r_observer * cos_theta_observer - photon[3];
			photon[0] = photon_time = std::hypot(dx, photon[2], dz);
			if (photon[0] == 0.)
				return Status::DIVISION_BY_ZERO;
			alpha = photon[2];
			beta = photon[3] * sin_theta_observer - photon[1] * cos_theta_observer;
			photon[4] = 1.;
			const Type distance_1 = 1. / photon[0];
			photon[5] = dx * distance_1;
			photon[6] = photon[2] * distance_1;
			photon[7] = dz * distance_1;
			return CartesianToSpherical(photon);
		}
		Type Energy(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			const Type r_1 = 1. / y[1], r_2 = Power2(r_1), r_3 = r_1 * r_2, r_4 = Power2(r_2);
			const Type v2 = Power2(y[5]) + Power2(y[1]) * (Power2(y[6]) + Power2(std::sin(y[2]) * y[7])), v4 = Power2(v2), v6 = v2 * v4, v8 = Power2(v4);
			const Type rdot = y[5], rdot2 = Power2(rdot);
			Type E = 0.5 * v2 - r_1;
			if (PN_ & 1)
				E += 0.5 * r_2 + 0.375 * v4 + 1.5 * r_1 * v2;
			if (PN_ & 2)
				E += -0.5 * r_3 + 0.3125 * v6 + 1.75 * r_2 * v2 + 0.5 * r_2 * rdot2 + 2.625 * r_1 * v4;
			if (PN_ & 8)
				E += 0.375 * r_4 + 1.25 * r_3 * v2 + 1.5 * r_3 * rdot2 + 0.2734375 * v8 + 8.4375 * r_2 * v4 + 0.75 * r_2 * v2 * rdot2 + 3.4375 * r_1 * v6;
			return E;
		}
		Type AngularMomentum(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			const Type r_1 = 1. / y[1], r_2 = Power2(r_1), r_3 = Power3(r_1);
			const Type v_tan2 = Power2(y[1]) * (Power2(y[6]) + Power2(std::sin(y[2]) * y[7]));
			const Type v2 = Power2(y[5]) + v_tan2, v4 = Power2(v2), v6 = v2 * v4;
			const Type rdot = y[5], rdot2 = Power2(rdot);
			Type eff = 1.;
			if (PN_ & 1)
				eff += 3. * r_1 + 0.5 * v2;
			if (PN_ & 2)
				eff += 3.5 * r_2 + 0.375 * v4 + 3.5 * r_1 * v2;
			if (PN_ & 8)
				eff += 2.5 * r_3 + 0.3125 * v6 + 11.25 * r_2 * v2 + 0.5 * r_2 * rdot2 + 4.125 * r_1 * v4;
			return y[1] * std::sqrt(v_tan2) * eff;
		}
		Type CarterConstant(const boost::numeric::ublas::bounded_vector<Type, 8> &y, const Type mu2, TimeSystem time, DynamicalSystem dynamics) override {
			return Power2(AngularMomentum(y, time, dynamics));
		}
		Type Redshift(const boost::numeric::ublas::bounded_vector<Type, 8> &y, const boost::numeric::ublas::bounded_vector<Type, 8> &photon_position, TimeSystem object_time, TimeSystem photon_time) override {
			namespace ublas = boost::numeric::ublas;
			const Type delta_epsilon = 1. - 2. / y[1];
			const ublas::vector_range<const ublas::bounded_vector<Type, 8>> velocity(y, ublas::range(4, 8));
			return (1. - this->DotProduct(y, velocity, ublas::vector_range<const ublas::bounded_vector<Type, 8>>(photon_position, ublas::range(4, 8)), false) / std::sqrt(delta_epsilon)) / std::sqrt(delta_epsilon - this->DotProduct(y, velocity, velocity, false));
		}
		int NormalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			y[4] = 1;
			if (y[5] * y[5] + y[1] * y[1] * (y[6] * y[6] + Power2(std::sin(y[2]) * y[7])) >= 1.)
				return Status::FAILURE;
			return Status::SUCCESS;
		}
		int NormalizeNullGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y, Type frequency = 1.) override {
			y[4] = 1;
			const Type v_1 = 1. / std::sqrt(Power2(y[5]) + Power2(y[1]) * (Power2(y[6]) + Power2(std::sin(y[2]) * y[7])));
			for (int i = 5; i < 8; ++i)
				y[i] *= v_1;
			return Status::SUCCESS;
		}
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> GetIntegrationSystem(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override {
			if (time != T || dynamics != LAGRANGIAN || motion != GEODESIC) {
				throw std::invalid_argument(fmt::format("Newton::GetIntegrationSystem({}, {}, {}) invaild", time, dynamics, motion));
			}
			return [PN = this->PN_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
				const Type r_1 = 1. / y[1], r_2 = Power2(r_1), sin_theta = SinTheta(y[2]), cos_theta = CosTheta(y[2]);
				const Type v_tan2 = Power2(y[1]) * (Power2(y[6]) + Power2(sin_theta * y[7]));
				const Type rdot = y[5], rdot2 = Power2(y[5]);
				const Type v2 = rdot2 + v_tan2;
				dydt[0] = y[4];
				dydt[1] = y[5];
				dydt[2] = y[6];
				dydt[3] = y[7];
				dydt[4] = 0.;
				Type A = 1., B = 0.;
				if (PN & 1) {
					A += v2 - 4. * r_1;
					B += -4. * rdot;
				}
				if (PN & 2) {
					A += r_1 * (-2. * rdot2 + 9. * r_1);
					B += 2. * r_1 * rdot;
				}
				if (PN & 8) {
					A += -16. * Power3(r_1) + rdot2 * r_2;
					B += -4. * rdot * r_2;
				}
				const Type B_r2 = B * r_2;
				dydt[5] = r_1 * v_tan2 - r_2 * A - B_r2 * y[5];
				dydt[6] = sin_theta * cos_theta * Power2(y[7]) - (2. * y[5] * r_1 + B_r2) * y[6];
				dydt[7] = -(2. * cos_theta / sin_theta * y[6] + 2. * y[5] * r_1 + B_r2) * y[7];
			};
		}
	};
	template <typename Type>
	class PN1 : public Newton<Type> {
	  protected:
		const Type PN1_;

	  public:
		PN1(Type fSP) : Newton<Type>(1), PN1_(fSP) {};
		std::string Name() const override {
			return "Newton-PN1";
		}
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> GetIntegrationSystem(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override {
			if (time != T || dynamics != LAGRANGIAN || motion != GEODESIC) {
				throw std::invalid_argument(fmt::format("PN1::GetIntegrationSystem({}, {}, {}) invaild", time, dynamics, motion));
			}
			return [PN1 = this->PN1_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
				const Type r_1 = 1. / y[1], r_2 = r_1 * r_1, sin_theta = SinTheta(y[2]), cos_theta = CosTheta(y[2]);
				const Type v_tan2 = y[1] * y[1] * (y[6] * y[6] + Power2(sin_theta * y[7]));
				const Type v2 = y[5] * y[5] + v_tan2;
				dydt[0] = y[4];
				dydt[1] = y[5];
				dydt[2] = y[6];
				dydt[3] = y[7];
				dydt[4] = 0.;
				const Type A = 1. + (v2 - 4. * r_1) * PN1, B_r2 = -4. * y[5] * PN1 * r_2;
				dydt[5] = r_1 * v_tan2 - r_2 * A - B_r2 * y[5];
				dydt[6] = sin_theta * cos_theta * y[7] * y[7] - (2. * y[5] * r_1 + B_r2) * y[6];
				dydt[7] = -(2. * cos_theta / sin_theta * y[6] + 2. * y[5] * r_1 + B_r2) * y[7];
			};
		}
	};
	template <typename Type>
	class Schwarzschild : public Metric<Type> {
	  protected:
		Type DotProduct3D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x1, Type x2, Type x3, Type y1, Type y2, Type y3) override {
			return position[1] * x1 * y1 / (position[1] - 2.) + position[1] * position[1] * x2 * y2 + Power2(position[1] * std::sin(position[2])) * x3 * y3;
		}
		Type DotProduct4D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x0, Type x1, Type x2, Type x3, Type y0, Type y1, Type y2, Type y3) override {
			return -(1. - 2. / position[1]) * x0 * y0 + position[1] * x1 * y1 / (position[1] - 2.) + position[1] * position[1] * x2 * y2 + Power2(position[1] * std::sin(position[2])) * x3 * y3;
		}

	  public:
		Schwarzschild() {}
		std::string Name() const override {
			return "Schwarzschild";
		}
		int MetricTensor(const boost::numeric::ublas::bounded_vector<Type, 8> &position, boost::numeric::ublas::bounded_matrix<Type, 4, 4> &metric) override {
			metric.clear();
			metric(0, 0) = -(1. - 2. / position[1]);
			metric(1, 1) = position[1] / (position[1] - 2.);
			metric(2, 2) = position[1] * position[1];
			metric(3, 3) = Power2(position[1] * std::sin(position[2]));
			return position[1] == 2. ? Status::DIVISION_BY_ZERO : Status::SUCCESS;
		}
		Type DistanceSquare(const boost::numeric::ublas::bounded_vector<Type, 8> &x, const boost::numeric::ublas::bounded_vector<Type, 8> &y, bool include_time) override {
			if (include_time)
				return -(1. - 2. / x[1]) * Power2(x[0] - y[0]) + x[1] * Power2(x[1] - y[1]) / (x[1] - 2.) + Power2(x[1] * (x[2] - y[2])) + Power2(x[1] * std::sin(x[2]) * PhiDifference(x[3] - y[3]));
			return x[1] * Power2(x[1] - y[1]) / (x[1] - 2.) + Power2(x[1] * (x[2] - y[2])) + Power2(x[1] * std::sin(x[2]) * PhiDifference(x[3] - y[3]));
		}
		int LagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type dt_dtau = 1. / y[4], r2 = Power2(y[1]);
			y[4] = (y[4] - 1. + 2. / y[1]) * dt_dtau; // 1 + p_t
			y[5] *= y[1] / (y[1] - 2.) * dt_dtau;
			y[6] *= r2 * dt_dtau;
			y[7] *= r2 * Power2(std::sin(y[2])) * dt_dtau;
			return Status::SUCCESS;
		}
		int HamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type r_1 = 1. / y[1], g11_1 = 1. - 2. * r_1;
			y[4] = -g11_1 / (y[4] - 1.);
			y[5] *= g11_1 * y[4];
			y[6] *= r_1 * r_1 * y[4];
			y[7] *= Power2(r_1 / std::sin(y[2])) * y[4];
			return Status::SUCCESS;
		}
		int FastTrace(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type r_object, const Type theta_object, const Type phi_object, Type &alpha, Type &beta, boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time) override {
			const Type sin_theta_object = SinTheta(theta_object), cos_theta_object = CosTheta(theta_object), sin_phi_object = std::sin(phi_object), cos_phi_object = std::cos(phi_object);
			const Type cos_observer_object = sin_theta_observer * sin_theta_object * cos_phi_object + cos_theta_observer * cos_theta_object;
			const Type delta_phi = acos(cos_observer_object), sin_observer_object = std::sqrt(1. - cos_observer_object * cos_observer_object);
			const Type u0 = 1. / r_observer, u1 = 1. / r_object, u12 = u1 * u1;
			const Type g11_1 = 1. - 2. * u1;
			if (delta_phi == 0.) {
				alpha = 0.;
				beta = 0.;
				photon(0) = (r_object - r_observer + 4. * log(r_object * u0) + 4. * (u1 - u0)) * r_observer / (r_observer - 2.);
				photon(1) = r_object;
				photon(2) = theta_object;
				photon(3) = phi_object;
				photon(4) = r_object * (r_observer - 2.) / (r_observer * (r_object - 2.));
				photon(5) = (r_object - 2.) * u1;
				photon(6) = 0.;
				photon(7) = 0.;
				photon_time = r_object - r_observer + 2. * log(r_object * u0);
				return Status::SUCCESS;
			}
			Type impact_parameter_upper_limit = r_object / std::sqrt(g11_1), turning_phi;
			if (Type x02[2]; u1 * 3. < 1. && PolySolveQuadratic(2., -g11_1, -g11_1 * u1, x02) == 2) // remove (x-u1) as cubic_roots[1] = u1
				turning_phi = boost::math::constants::half_root_two<Type>() * EllipticIntegral(0, u0, u1, 0., 0., u1, -1., -x02[0], 1., x02[1], -1.);
			else { // there is no turning point on the trajectory
				impact_parameter_upper_limit = M_SQRT27 - boost::math::tools::epsilon<Type>();
				turning_phi = boost::math::constants::pi<Type>(); // make delta_phi <= turning_phi
			}
			std::pair<Type, Type> impact_root;
			std::array<Type, 3> cubic_roots;
			try {
				if (delta_phi > turning_phi)
					impact_root = boost::math::tools::bisect(
						[u0, u1, delta_phi, &cubic_roots](Type l) -> Type {
							if (l == 0.)
								return -delta_phi;
							Type half_l_2 = 0.5 / (l * l);
							if (PolySolveCubic(-0.5, 0., half_l_2, cubic_roots.begin()) == 3) {
								if (cubic_roots[1] < u1) // cubic_roots[1] < u1 < cubic_roots[1] + EPSILON
									return boost::math::constants::half_root_two<Type>() * EllipticIntegral(0, u0, u1, 0., 1., u1, -1., -cubic_roots[0], 1., cubic_roots[2], -1.) - delta_phi;
								else // u1 -> turning point -> u1 -> u0
									return boost::math::constants::root_two<Type>() * EllipticIntegral(0., u1, cubic_roots[1], 0., 1., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.) + boost::math::constants::half_root_two<Type>() * EllipticIntegral(0, u0, u1, 0., 1., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.) - delta_phi;
							}
							// impact_parameter < std::sqrt(27)
							cubic_roots[1] = std::nan("");
							return boost::math::constants::half_root_two<Type>() * EllipticIntegral2Complex(0, u0, u1, 0., 1., -half_l_2 / cubic_roots[0], cubic_roots[0] - 0.5, 1., -cubic_roots[0], 1.) - delta_phi;
						},
						Type{M_SQRT27} + boost::math::tools::root_epsilon<Type>(), impact_parameter_upper_limit, boost::math::tools::eps_tolerance<Type>());
				else
					impact_root = boost::math::tools::bisect(
						[u0, u1, delta_phi, &cubic_roots](Type l) -> Type {
							if (l == 0.)
								return -delta_phi;
							Type half_l_2 = 0.5 / (l * l);
							if (PolySolveCubic(-0.5, 0., half_l_2, cubic_roots.begin()) == 3) {
								if (cubic_roots[1] < u1) // cubic_roots[1] < u1 < cubic_roots[1] + EPSILON
									return boost::math::constants::half_root_two<Type>() * EllipticIntegral(0, u0, u1, 0., 1., u1, -1., -cubic_roots[0], 1., cubic_roots[2], -1.) - delta_phi;
								else // u1 -> u0
									return boost::math::constants::half_root_two<Type>() * EllipticIntegral(0, u0, u1, 0., 1., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.) - delta_phi;
							}
							// impact_parameter < std::sqrt(27)
							cubic_roots[1] = std::nan("");
							return boost::math::constants::half_root_two<Type>() * EllipticIntegral2Complex(0, u0, u1, 0., 1., -half_l_2 / cubic_roots[0], cubic_roots[0] - 0.5, 1., -cubic_roots[0], 1.) - delta_phi;
						},
						0., impact_parameter_upper_limit, boost::math::tools::eps_tolerance<Type>());
			} catch (const std::exception &e) {
				return Status::FAILURE;
			}
			Type impact_root_value = 0.5 * (impact_root.first + impact_root.second), impact_root_value2 = impact_root_value * impact_root_value;
			alpha = impact_root_value / sin_observer_object * sin_theta_object * sin_phi_object;
			beta = impact_root_value / sin_observer_object * (cos_theta_object * sin_theta_observer - sin_theta_object * cos_phi_object * cos_theta_observer);
			if (!isnan(cubic_roots[1])) {
				if (cubic_roots[1] < u1) { // cubic_roots[1] < u1 < cubic_roots[1] + EPSILON
					const Type ellip_int_4 = EllipticIntegral(-4, u0, u1, 0., 1., u1, -1., -cubic_roots[0], 1., cubic_roots[2], -1.);
					photon(0) = -boost::math::constants::half_root_two<Type>() * ellip_int_4 / (impact_root_value * (1. - 2. * u0));
					photon_time = -boost::math::constants::half_root_two<Type>() * (ellip_int_4 + 2. * EllipticIntegral(-2, u0, u1, 0., 1., u1, -1., -cubic_roots[0], 1., cubic_roots[2], -1.) + 4. * EllipticIntegral(-2, u0, u1, 1., -2., u1, -1., -cubic_roots[0], 1., cubic_roots[2], -1.)) / impact_root_value;
				} else if (delta_phi > turning_phi) { // u1 -> turning point -> u1 -> u0
					const Type ellip_int_4 = 2. * EllipticIntegral(-4, u1, cubic_roots[1], 0., 1., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.) + EllipticIntegral(-4, u0, u1, 0., 1., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.);
					photon(0) = -boost::math::constants::half_root_two<Type>() * ellip_int_4 / (impact_root_value * (1. - 2. * u0));
					photon_time = -boost::math::constants::half_root_two<Type>() * (ellip_int_4 + 4. * EllipticIntegral(-2, u1, cubic_roots[1], 0., 1., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.) + 2. * EllipticIntegral(-2, u0, u1, 0., 1., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.) + 8. * EllipticIntegral(-2, u1, cubic_roots[1], 1., -2., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.) + 4. * EllipticIntegral(-2, u0, u1, 1., -2., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.)) / impact_root_value;
				} else { // u1 -> u0
					const Type ellip_int_4 = EllipticIntegral(-4, u0, u1, 0., 1., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.);
					photon(0) = -boost::math::constants::half_root_two<Type>() * ellip_int_4 / (impact_root_value * (1. - 2. * u0));
					photon_time = -boost::math::constants::half_root_two<Type>() * (ellip_int_4 + 2. * EllipticIntegral(-2, u0, u1, 0., 1., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.) + 4. * EllipticIntegral(-2, u0, u1, 1., -2., cubic_roots[1], -1., -cubic_roots[0], 1., cubic_roots[2], -1.)) / impact_root_value;
				}
			} else { // impact_parameter < std::sqrt(27)
				const Type ellip_int_4 = EllipticIntegral2Complex(-4, u0, u1, 0., 1., -0.5 / (impact_root_value2 * cubic_roots[0]), cubic_roots[0] - 0.5, 1., -cubic_roots[0], 1.);
				photon(0) = -boost::math::constants::half_root_two<Type>() * ellip_int_4 / (impact_root_value * (1. - 2. * u0));
				photon_time = -boost::math::constants::half_root_two<Type>() * (ellip_int_4 + 2. * EllipticIntegral2Complex(-2, u0, u1, 0., 1., -0.5 / (impact_root_value2 * cubic_roots[0]), cubic_roots[0] - 0.5, 1., -cubic_roots[0], 1.) + 4. * EllipticIntegral2Complex(-2, u0, u1, 1., -2., -0.5 / (impact_root_value2 * cubic_roots[0]), cubic_roots[0] - 0.5, 1., -cubic_roots[0], 1.)) / impact_root_value;
			}
			photon(1) = r_object;
			photon(2) = theta_object;
			photon(3) = phi_object;
			photon(4) = g11_1 / (1. - 2. * u0); // dt/d\tau = 1. when r = r_observer
			if (delta_phi > turning_phi)
				photon(5) = -g11_1 * SquareRoot(1. - g11_1 * u12 * impact_root_value2);
			else
				photon(5) = g11_1 * SquareRoot(1. - g11_1 * u12 * impact_root_value2);
			// the direction of the angular momentum is [alpha * cos_theta_observer, beta, -alpha * sin_theta_observer],
			// the y component should have the same sign as the component perpendicular to [std::cos(phi), std::sin(phi), 0],
			// which is also the component along the [-std::sin(phi), std::cos(phi), 0] direction.
			if (alpha == 0.) {
				photon(6) = std::copysign(g11_1 * u12 * impact_root_value, beta * cos_phi_object - alpha * cos_theta_observer * sin_phi_object);
				photon(7) = 0.;
			} else {
				photon(6) = std::copysign(g11_1 * u12 * SquareRoot(impact_root_value2 - Power2(alpha * sin_theta_observer / sin_theta_object)), beta * cos_phi_object - alpha * cos_theta_observer * sin_phi_object);
				photon(7) = -g11_1 * alpha * sin_theta_observer * u12 / (sin_theta_object * sin_theta_object);
			}
			return Status::SUCCESS;
		}
		int FastShadow(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type alpha, const Type beta, const Type r_min) override {
			if (alpha * alpha + beta * beta <= 27.)
				return 0;
			return 1;
		}
		Type Energy(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				if (time == T)
					return (y[1] - 2.) / (y[1] * y[4]);
				// time == TAU
				return (1. - 2. / y[1]) * y[4];
			} // dynamics == HAMILTONIAN
			return 1. - y[4];
		}
		Type AngularMomentum(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				if (time == T)
					return Power2(y[1] * std::sin(y[2])) * y[7] / y[4];
				// time == TAU
				return Power2(y[1] * std::sin(y[2])) * y[7];
			} // dynamics == HAMILTONIAN
			return y[7];
		}
		Type CarterConstant(const boost::numeric::ublas::bounded_vector<Type, 8> &y, const Type mu2, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				if (time == T)
					return Power4(y[1]) * (Power2(y[6]) + Power2(y[7] * std::cos(y[2]) * std::sin(y[2]))) / Power2(y[4]);
				// time == TAU
				return Power4(y[1]) * (Power2(y[6]) + Power2(y[7] * std::cos(y[2]) * std::sin(y[2])));
			} // dynamics == HAMILTONIAN
			return Power2(y[6]) + Power2(y[7] / tan(y[2]));
		}
		int NormalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type g11_1 = 1. - 2. / y[1];
			if (g11_1 <= 0)
				return 1;
			y[4] = std::sqrt(g11_1 - (Power2(y[5]) / g11_1 + Power2(y[1] * y[6]) + Power2(y[1] * std::sin(y[2]) * y[7])));
			return isnan(y[4]) ? Status::DOMAIN_ERROR : Status::SUCCESS;
		}
		int NormalizeNullGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y, Type frequency = 1.) override {
			const Type g11_1 = 1. - 2. / y[1];
			if (g11_1 <= 0)
				return 1;
			const Type coefficient = std::copysign(g11_1, frequency) / std::sqrt(y[5] * y[5] + g11_1 * (Power2(y[1] * y[6]) + Power2(y[1] * std::sin(y[2]) * y[7])));
			if (isnan(coefficient))
				return Status::DOMAIN_ERROR;
			y[4] = frequency;
			y[5] *= coefficient;
			y[6] *= coefficient;
			y[7] *= coefficient;
			return Status::SUCCESS;
		}
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> GetIntegrationSystem(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override {
			if (time == T) {
				if (dynamics == LAGRANGIAN) {
					if (motion == GEODESIC)
						return [](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = y[6]; // d\theta/dt
							dydt[3] = y[7]; // d\phi/dt
							const Type r = y[1], sin_theta = SinTheta(y[2]), cos_theta = CosTheta(y[2]);
							const Type rm2 = r - 2., rm3 = r - 3.;
							const Type rm2r_1 = 1. / (rm2 * r);
							// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
							dydt[4] = 2. * y[5] * rm2r_1 * y[4];
							// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[5] = -rm2 / Power3(r) + 3. * rm2r_1 * Power2(y[5]) + rm2 * (Power2(y[6]) + Power2(sin_theta * y[7]));
							// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[6] = -2. * rm3 * rm2r_1 * y[5] * y[6] + sin_theta * cos_theta * Power2(y[7]);
							// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							if (y[7] == 0.)
								dydt[7] = 0.;
							else if (sin_theta == 0.)
								return Status::DIVISION_BY_ZERO;
							else
								dydt[7] = -2. * (rm3 * rm2r_1 * y[5] + cos_theta / sin_theta * y[6]) * y[7];
						};
					else if (motion == CIRCULAR)
						return [](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = 0.;	// dr/dt
							dydt[2] = 0.;	// d\theta/dt
							dydt[3] = y[7]; // d\phi/dt
							dydt[4] = 0.;
							dydt[5] = 0.;
							dydt[6] = 0.;
							dydt[7] = 0.;
						};
					else if (motion == HELICAL)
						return [](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
							const Type r = y[1];
							if (r <= 2.)
								throw std::domain_error("r <= 2.");
							const Type sin2_theta = Power2(std::sin(y[2]));
							const Type g00 = -1. + 2. / r, g11 = -1. / g00, g33 = r * r * sin2_theta;
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = 0.;	// d\theta/dt = 0.
							dydt[3] = y[7]; // d\phi/dt
							// d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt
							dydt[4] = sin2_theta * y[5] * (r * (-g00 - g11 * y[5] * y[5] - y[4] * y[4]) + 1. + Power2(g11 * y[5])) / (y[4] * (g33 + Power2(g33 * y[7] / y[4])));
							// d^2r/dt^2 = 0.
							dydt[5] = 0.;
							// d^2\theta/dt^2 = 0.
							dydt[6] = 0.;
							// d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
							dydt[7] = y[7] * (-2. / r * y[5] + dydt[4] / y[4]);
						};
					else if (motion == CYLINDRICAL)
						return [](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
							const Type r = y[1];
							if (r <= 2.)
								throw std::domain_error("r <= 2.");
							const Type r_1 = 1. / r;
							const Type g11 = 1. / (1. - 2. * r_1);
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = y[6]; // d\theta/dt.
							dydt[3] = y[7]; // d\phi/dt
							// d^2r/dt^2 = r * (d\theta/dt)^2
							dydt[5] = r * y[6] * y[6];
							// d^2\theta/dt^2 =
							dydt[6] = -2. * r_1 * y[5] * y[6];
							// d^2\phi/dt^2 = 0.
							dydt[7] = 0.;
							// d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt
							dydt[4] = (r_1 * (1. + Power2(g11 * y[5])) - 2. * g11 * dydt[5]) * r_1 * y[5] / y[4];
						};
				} else if (dynamics == HAMILTONIAN && motion == GEODESIC)
					return [](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
						const Type r_1 = 1. / y[1], r_2 = Power2(r_1), g11_1 = 1. - 2. * r_1, E = 1. - y[4], L2 = Power2(y[7]);
						const Type sin_1_theta = 1. / SinTheta(y[2]), sin_2_theta = Power2(sin_1_theta);
						//[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
						dydt[0] = g11_1 / E;						  // d\tau/dt
						dydt[1] = g11_1 * y[5] * dydt[0];			  // dr/dt
						dydt[2] = y[6] * r_2 * dydt[0];				  // d\theta/dt
						dydt[3] = y[7] * sin_2_theta * r_2 * dydt[0]; // d\phi/dt
						dydt[4] = 0.;
						dydt[5] = (-(Power2(y[5]) + Power2(E) / Power2(g11_1)) + (Power2(y[6]) + L2 * sin_2_theta) * r_1) * r_2 * dydt[0];
						dydt[6] = sin_2_theta * L2 * CosTheta(y[2]) * sin_1_theta * r_2 * dydt[0];
						dydt[7] = 0.;
					};
			} else if (time == TAU && dynamics == LAGRANGIAN && motion == GEODESIC)
				return [](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
					dydt[0] = y[4]; // dt/d\tau
					dydt[1] = y[5]; // dr/d\tau
					dydt[2] = y[6]; // d\theta/d\tau
					dydt[3] = y[7]; // d\phi/d\tau
					const Type r = y[1], sin_theta = SinTheta(y[2]), cos_theta = CosTheta(y[2]);
					const Type rm2 = r - 2., r_1 = 1. / r;
					const Type rm2r_1 = 1. / (rm2 * r);
					// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
					dydt[4] = -2. * y[5] * rm2r_1 * y[4];
					// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
					dydt[5] = -rm2 * Power3(r_1) * Power2(y[4]) + rm2r_1 * Power2(y[5]) + rm2 * (Power2(y[6]) + Power2(sin_theta * y[7]));
					// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
					dydt[6] = -2. * r_1 * y[5] * y[6] + sin_theta * cos_theta * Power2(y[7]);
					// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
					dydt[7] = -2. * (r_1 * y[5] + cos_theta / sin_theta * y[6]) * y[7];
				};
			throw std::invalid_argument(fmt::format("Schwarzschild::GetIntegrationSystem({}, {}, {}) invaild", time, dynamics, motion));
		}
	};
	template <typename Type>
	class ReissnerNordstrom : public Metric<Type> {
	  protected:
		Type DotProduct3D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x1, Type x2, Type x3, Type y1, Type y2, Type y3) override {
			const Type r_1 = 1. / position[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
			return x1 * y1 / g11_1 + position[1] * position[1] * x2 * y2 + Power2(position[1] * std::sin(position[2])) * x3 * y3;
		}
		Type DotProduct4D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x0, Type x1, Type x2, Type x3, Type y0, Type y1, Type y2, Type y3) override {
			const Type r_1 = 1. / position[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
			return -g11_1 * x0 * y0 + x1 * y1 / g11_1 + position[1] * position[1] * x2 * y2 + Power2(position[1] * std::sin(position[2])) * x3 * y3;
		}

	  public:
		const Type r_Q_, r_Q2_, r_Q4_;
		ReissnerNordstrom(Type charge) : r_Q_(0.5 * std::sqrt(M_1_PI) * charge), r_Q2_(r_Q_ * r_Q_), r_Q4_(r_Q2_ * r_Q2_) {}
		std::string Name() const override {
			return "Reissner-Nordstrom";
		}
		int MetricTensor(const boost::numeric::ublas::bounded_vector<Type, 8> &position, boost::numeric::ublas::bounded_matrix<Type, 4, 4> &metric) override {
			const Type r_1 = 1. / position[1];
			metric.clear();
			metric(0, 0) = -1. + (2. - r_Q2_ * r_1) * r_1;
			metric(1, 1) = -1. / metric(0, 0);
			metric(2, 2) = position[1] * position[1];
			metric(3, 3) = Power2(position[1] * std::sin(position[2]));
			return position[1] == 2. ? Status::DIVISION_BY_ZERO : Status::SUCCESS;
		}
		Type DistanceSquare(const boost::numeric::ublas::bounded_vector<Type, 8> &x, const boost::numeric::ublas::bounded_vector<Type, 8> &y, bool include_time) override {
			const Type r_1 = 1. / x[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
			if (include_time)
				return -g11_1 * Power2(x[0] - y[0]) + Power2(x[1] - y[1]) / g11_1 + Power2(x[1]) * (Power2(x[2] - y[2]) + Power2(std::sin(x[2]) * PhiDifference(x[3] - y[3])));
			return Power2(x[1] - y[1]) / g11_1 + Power2(x[1]) * (Power2(x[2] - y[2]) + Power2(std::sin(x[2]) * PhiDifference(x[3] - y[3])));
		}
		int LagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type dt_dtau = 1. / y[4], r_1 = 1. / y[1], r2 = Power2(y[1]);
			const Type rs_rQ2 = (2. - r_Q2_ * r_1) * r_1;
			y[4] = (y[4] - 1. + rs_rQ2) * dt_dtau; // 1 + p_t
			y[5] *= dt_dtau / (1. - rs_rQ2);
			y[6] *= r2 * dt_dtau;
			y[7] *= r2 * Power2(std::sin(y[2])) * dt_dtau;
			return Status::SUCCESS;
		}
		int HamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type r_1 = 1. / y[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
			y[4] = -g11_1 / (y[4] - 1.);
			y[5] *= g11_1 * y[4];
			y[6] *= Power2(r_1) * y[4];
			y[7] *= Power2(r_1 / std::sin(y[2])) * y[4];
			return Status::SUCCESS;
		}
		int FastTrace(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type r_object, const Type theta_object, const Type phi_object, Type &alpha, Type &beta, boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time) override {
			return Status::FAILURE;
		}
		Type Energy(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				const Type r_1 = 1. / y[1];
				if (time == T)
					return (1. - (2. - r_Q2_ * r_1) * r_1) / y[4];
				// time == TAU
				return (1. - (2. - r_Q2_ * r_1) * r_1) * y[4];
			} // dynamics == HAMILTONIAN
			return 1. - y[4];
		}
		Type AngularMomentum(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				if (time == T)
					return Power2(y[1] * std::sin(y[2])) * y[7] / y[4];
				// time == TAU
				return Power2(y[1] * std::sin(y[2])) * y[7];
			} // dynamics == HAMILTONIAN
			return y[7];
		}
		Type CarterConstant(const boost::numeric::ublas::bounded_vector<Type, 8> &y, const Type mu2, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				if (time == T)
					return Power4(y[1]) * (Power2(y[6]) + Power2(y[7] * std::cos(y[2]) * std::sin(y[2]))) / Power2(y[4]);
				// time == TAU
				return Power4(y[1]) * (Power2(y[6]) + Power2(y[7] * std::cos(y[2]) * std::sin(y[2])));
			} // dynamics == HAMILTONIAN
			return Power2(y[6]) + Power2(y[7] / tan(y[2]));
		}
		int NormalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type r_1 = 1. / y[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
			if (g11_1 <= 0)
				return 1;
			y[4] = std::sqrt(g11_1 - (Power2(y[5]) / g11_1 + Power2(y[1] * y[6]) + Power2(y[1] * std::sin(y[2]) * y[7])));
			return isnan(y[4]) ? Status::DOMAIN_ERROR : Status::SUCCESS;
		}
		int NormalizeNullGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y, Type frequency = 1.) override {
			const Type r_1 = 1. / y[1], g11_1 = 1. - (2. - r_Q2_ * r_1) * r_1;
			if (g11_1 <= 0)
				return 1;
			const Type coefficient = std::copysign(g11_1, frequency) / std::sqrt(Power2(y[5]) + g11_1 * (Power2(y[1] * y[6]) + Power2(y[1] * std::sin(y[2]) * y[7])));
			if (isnan(coefficient))
				return Status::DOMAIN_ERROR;
			y[4] = frequency;
			y[5] *= coefficient;
			y[6] *= coefficient;
			y[7] *= coefficient;
			return Status::SUCCESS;
		}
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> GetIntegrationSystem(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override {
			if (time != T || dynamics != LAGRANGIAN || motion != GEODESIC)
				throw std::invalid_argument(fmt::format("ReissnerNordstrom::GetIntegrationSystem({}, {}, {}) invaild", time, dynamics, motion));
			return [r_Q2 = this->r_Q2_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
				dydt[0] = y[4]; // d\tau/dt
				dydt[1] = y[5]; // dr/dt
				dydt[2] = y[6]; // d\theta/dt
				dydt[3] = y[7]; // d\phi/dt
				const Type r = y[1], r_1 = 1. / r, r_2 = Power2(r_1);
				const Type sin_theta = SinTheta(y[2]), cos_theta = CosTheta(y[2]);
				const Type g11_1 = 1. - (2. - r_Q2 * r_1) * r_1, g11 = 1. / g11_1;
				const Type m_r_Q_r2 = (1. - r_Q2 * r_1) * r_2, m_r_Q_r2_g11 = m_r_Q_r2 * g11;
				dydt[4] = 2. * m_r_Q_r2_g11 * y[5] * y[4];
				dydt[5] = -g11_1 * m_r_Q_r2 + 3. * m_r_Q_r2_g11 * Power2(y[5]) + r * g11_1 * (Power2(y[6]) + Power2(sin_theta * y[7]));
				dydt[6] = -2. * (g11_1 * r_1 - m_r_Q_r2) * g11 * y[5] * y[6] + sin_theta * cos_theta * Power2(y[7]);
				if (y[7] == 0.)
					dydt[7] = 0.;
				else if (sin_theta == 0.)
					throw std::invalid_argument("Division by zero");
				else
					dydt[7] = -2. * ((r_1 - m_r_Q_r2_g11) * y[5] + cos_theta / sin_theta * y[6]) * y[7];
			};
		}
	};

	template <typename Type>
	struct KerrFastTraceParameters;

	template <typename Type>
	class Kerr : public Metric<Type> {
	  protected:
		Type DotProduct3D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x1, Type x2, Type x3, Type y1, Type y2, Type y3) override {
			const Type r = position[1], r2 = Power2(r), Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * Power2(std::cos(position[2])), r_rho_2 = 2. * r / rho2, sin2_theta = Power2(std::sin(position[2]));
			return (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * sin2_theta * sin2_theta) * x3 * y3;
		}
		Type DotProduct4D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x0, Type x1, Type x2, Type x3, Type y0, Type y1, Type y2, Type y3) override {
			const Type r = position[1], r2 = Power2(r), Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * Power2(std::cos(position[2])), r_rho_2 = 2. * r / rho2, sin2_theta = Power2(std::sin(position[2]));
			return (r_rho_2 - 1.) * x0 * y0 - r_rho_2 * a_ * sin2_theta * (x0 * y3 + x3 * y0) + (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * sin2_theta * sin2_theta) * x3 * y3;
		}
		static int DeltaUMuPhi(const boost::numeric::ublas::bounded_vector<Type, 2> &alpha_beta, boost::numeric::ublas::bounded_vector<Type, 2> &delta_u_mu_phi, KerrFastTraceParameters<Type> &params) {
			// The photon in the observer's frame has the tetrad velocity: [1, r / R, beta / R, -alpha / R], where R = std::sqrt(r^2 + alpha^2 + beta^2).
			Type alpha = alpha_beta(0), beta = alpha_beta(1);
			if (!isfinite(alpha) || !isfinite(beta))
				return Status::NUMERIC_ERROR;
			boost::numeric::ublas::bounded_vector<Type, 8> photon;
			Type photon_time;
			if (int status = params.kerr->InitializePhoton(photon, photon_time, alpha, beta, params.r, params.r2, params.theta_obs, params.sin_theta_obs); status != Status::SUCCESS)
				return status;
			const Type a = params.kerr->a_, a2 = params.kerr->a2_;
			params.E = params.kerr->Energy(photon, T, LAGRANGIAN);
			const Type E_1 = 1. / params.E;
			params.L = params.kerr->AngularMomentum(photon, T, LAGRANGIAN);
			const Type l = params.L * E_1, l2 = Power2(l);
			params.Q = params.kerr->CarterConstant(photon, 0., T, LAGRANGIAN);
			const Type q2 = params.Q * Power2(E_1);
			Type I_u_0, I_u_1 = 0., I_u_plus_0, I_u_plus_1, I_u_minus_0, I_u_minus_1, I_u_2_0, I_u_2_1, I_u_4_0, I_u_4_1;
			if (UIntegral(a, a2, params.kerr->u_plus_1, params.kerr->u_minus_1, params.kerr->u_plus, params.kerr->u_minus, l, l2, q2, params.r, params.u_obs, params.u_obj, I_u_0, I_u_1, I_u_plus_0, I_u_plus_1, I_u_minus_0, I_u_minus_1, I_u_2_0, I_u_2_1, I_u_4_0, I_u_4_1) == Status::NO_INTERSECTION_WITH_TARGET_SURFACE) {
				// photon failed to hit the r_obj surface, the impact params need to be smaller
				delta_u_mu_phi(0) = std::min(Type{0.99999999}, I_u_0 - boost::math::tools::root_epsilon<Type>());
				return Status::SCALING_REQUIRED;
			}
			// M = q^2 + (a^2-q^2-l^2)*mu^2-a^2*mu^4
			Type M_minus_plus[2];
			Type delta_M, delta_M_plus;
			if (int root_num_M = PolySolveQuadratic(-a2, a2 - l2 - q2, q2, M_minus_plus); root_num_M == 0)
				return Status::FAILURE;
			delta_M = std::sqrt(Power2(a2 - l2 - q2) + 4. * a2 * q2) / a2;
			if (M_minus_plus[1] > 1.) {
				M_minus_plus[1] = 1.;
				delta_M_plus = 0.;
			} else if (M_minus_plus[1] > 0.9) {
				delta_M_plus = 2. * l2 / (a2 + l2 + q2 + std::sqrt(Power2(a2 + l2 + q2) - 4. * a2 * l2));
				M_minus_plus[1] = 1. - delta_M_plus;
			} else
				delta_M_plus = 1. - M_minus_plus[1];
			Type mu_plus, mu_minus;
			Type A, k, n;
			Type I_mu_0, I_mu_full_turn;
			Type I_t_mu_0, I_t_mu_full_turn;
			Type I_phi_mu_0, I_phi_mu_full_turn;
			// Here the first turning point is determined by `beta` of the observer.
			int alpha_1, alpha_2;
			if (int status = Mu0Integral(a, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, boost::math::sign(beta), params.mu_obs, mu_plus, mu_minus, A, k, n, alpha_1, I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn); status != Status::SUCCESS)
				return status;
			Type mu_f_0, mu_f_1, t_mu, phi_mu_0, phi_mu_1;
			if (int status = MuFIntegral(a, a2, l, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, mu_plus, mu_minus, I_u_0, I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn, A, k, n, alpha_1, alpha_2, mu_f_0, t_mu, phi_mu_0); status != Status::SUCCESS)
				return status;
			const Type phi_0 = phi_mu_0 + params.kerr->u_r * ((l * params.kerr->u_plus_1 + 2. * (a - l)) * I_u_plus_0 - (l * params.kerr->u_minus_1 + 2. * (a - l)) * I_u_minus_0);
			params.tau = -t_mu - I_u_4_0;
			params.t = params.tau - 2. * params.kerr->u_r * ((a * (a - l) + Power2(params.kerr->u_plus_1)) * I_u_plus_0 - (a * (a - l) + Power2(params.kerr->u_minus_1)) * I_u_minus_0 - 2. * std::sqrt(1. - a2) * I_u_2_0);
			params.u_dir = -1.;
			params.mu_dir = std::copysign(alpha_2, mu_plus);
			delta_u_mu_phi(0) = mu_f_0 - params.mu_obj;
			if (phi_0 >= 0. && params.phi_obj > boost::math::constants::half_pi<Type>())
				delta_u_mu_phi(1) = phi_0 + params.phi_obj - boost::math::constants::two_pi<Type>();
			else
				delta_u_mu_phi(1) = phi_0 + params.phi_obj;
			if (I_u_1 == 0.)
				return Status::SUCCESS;
			if (int status = MuFIntegral(a, a2, l, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, mu_plus, mu_minus, I_u_1, I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn, A, k, n, alpha_1, alpha_2, mu_f_1, t_mu, phi_mu_1); status != Status::SUCCESS)
				return status;
			const Type phi_1 = phi_mu_1 + params.kerr->u_r * ((l * params.kerr->u_plus_1 + 2. * (a - l)) * I_u_plus_1 - (l * params.kerr->u_minus_1 + 2. * (a - l)) * I_u_minus_1);
			const Type angle_diff_0 = Power2(mu_f_0 - params.mu_obj) + Power2(delta_u_mu_phi(1));
			const Type delta_phi_1 = phi_1 >= 0. && params.phi_obj > boost::math::constants::half_pi<Type>() ? phi_1 + params.phi_obj - boost::math::constants::two_pi<Type>() : phi_1 + params.phi_obj;
			const Type angle_diff_1 = Power2(mu_f_1 - params.mu_obj) + Power2(delta_phi_1);
			if (angle_diff_0 <= angle_diff_1)
				return Status::SUCCESS;
			params.tau = -t_mu - I_u_4_1;
			params.t = params.tau - 2. * params.kerr->u_r * ((a * (a - l) + Power2(params.kerr->u_plus_1)) * I_u_plus_1 - (a * (a - l) + Power2(params.kerr->u_minus_1)) * I_u_minus_1 - 2. * std::sqrt(1. - a2) * I_u_2_1);
			params.u_dir = 1.;
			params.mu_dir = std::copysign(alpha_2, mu_plus);
			delta_u_mu_phi(0) = mu_f_1 - params.mu_obj;
			delta_u_mu_phi(1) = delta_phi_1;
			return Status::SUCCESS;
		}
		static int UIntegral(Type a, Type a2, Type u_plus_1, Type u_minus_1, Type u_plus, Type u_minus, Type l, Type l2, Type q2, Type r_obs, Type u_obs, Type u_obj, Type &I_u_0, Type &I_u_1, Type &I_u_plus_0, Type &I_u_plus_1, Type &I_u_minus_0, Type &I_u_minus_1, Type &I_u_2_0, Type &I_u_2_1, Type &I_u_4_0, Type &I_u_4_1) {
			// U=1+[a^2−q^2−l^2]u^2+2[(a−l)^2+q^2]u^3−a^2q^2u^4
			const Type c = a2 - l2 - q2, d = 2. * (Power2(a - l) + q2), e = -a2 * q2;
			Type u_roots[4];
			if (std::abs(e) < absolute_accuracy) {
				// U=1+[a^2−l^2]u^2+2(a−l)^2u^3
				// This means the observer and the object are on the equatorial plane.
				if (a == l) { // U=1.
					I_u_0 = u_obj - u_obs;
					I_u_plus_0 = u_plus * log((u_obj - u_plus) / (u_obs - u_plus));
					I_u_minus_0 = u_minus * log((u_obj - u_minus) / (u_obs - u_minus));
					I_u_2_0 = log(u_obj * r_obs);
					I_u_4_0 = r_obs - 1. / u_obj;
					return Status::SUCCESS;
				}
				const Type d_1 = 1. / d, sqrt_d_1 = boost::math::constants::half_root_two<Type>() / std::abs(a - l);
				if (int root_num_U = PolySolveCubic(c * d_1, 0., d_1, u_roots); root_num_U == 1) {
					const Type f = -1. / (u_roots[0] * d), g = f / u_roots[0];
					I_u_0 = sqrt_d_1 * EllipticIntegral2Complex(0, u_obs, u_obj, 1., 0., f, g, 1., -u_roots[0], 1.);
					I_u_plus_0 = sqrt_d_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, 1., -u_plus_1, f, g, 1., -u_roots[0], 1.);
					I_u_minus_0 = sqrt_d_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, 1., -u_minus_1, f, g, 1., -u_roots[0], 1.);
					I_u_2_0 = sqrt_d_1 * EllipticIntegral2Complex(-2, u_obs, u_obj, 0., 1., f, g, 1., -u_roots[0], 1.);
					I_u_4_0 = sqrt_d_1 * EllipticIntegral2Complex(-4, u_obs, u_obj, 0., 1., f, g, 1., -u_roots[0], 1.);
					return Status::SUCCESS;
				}
				if (u_roots[1] < u_obj) {
					I_u_0 = u_roots[1] / u_obj;
					return Status::NO_INTERSECTION_WITH_TARGET_SURFACE;
				}
				if (u_roots[1] >= u_plus) { // photon falls into the BH
					I_u_0 = sqrt_d_1 * EllipticIntegral(0, u_obs, u_obj, 1., 0., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.);
					I_u_plus_0 = sqrt_d_1 * -EllipticIntegral(-2, u_obs, u_obj, 1., -u_plus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.);
					I_u_minus_0 = sqrt_d_1 * -EllipticIntegral(-2, u_obs, u_obj, 1., -u_minus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.);
					I_u_2_0 = sqrt_d_1 * EllipticIntegral(-2, u_obs, u_obj, 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.);
					I_u_4_0 = sqrt_d_1 * EllipticIntegral(-4, u_obs, u_obj, 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.);
					return Status::SUCCESS;
				}
				AMinusPlusB(sqrt_d_1 * EllipticIntegral(0, u_obs, u_roots[1], 1., 0., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.), sqrt_d_1 * EllipticIntegral(0, u_obj, u_roots[1], 1., 0., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.), I_u_0, I_u_1);
				AMinusPlusB(sqrt_d_1 * -EllipticIntegral(-2, u_obs, u_roots[1], 1., -u_plus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.), sqrt_d_1 * -EllipticIntegral(-2, u_obj, u_roots[1], 1., -u_plus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.), I_u_plus_0, I_u_plus_1);
				AMinusPlusB(sqrt_d_1 * -EllipticIntegral(-2, u_obs, u_roots[1], 1., -u_minus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.), sqrt_d_1 * -EllipticIntegral(-2, u_obj, u_roots[1], 1., -u_minus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.), I_u_minus_0, I_u_minus_1);
				AMinusPlusB(sqrt_d_1 * EllipticIntegral(-2, u_obs, u_roots[1], 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.), sqrt_d_1 * EllipticIntegral(-2, u_obj, u_roots[1], 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.), I_u_2_0, I_u_2_1);
				AMinusPlusB(sqrt_d_1 * EllipticIntegral(-4, u_obs, u_roots[1], 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.), sqrt_d_1 * EllipticIntegral(-4, u_obj, u_roots[1], 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1.), I_u_4_0, I_u_4_1);
				return Status::SUCCESS;
			}
			const Type e_1 = 1. / e, sqrt_e = std::sqrt(std::abs(e)), sqrt_e_1 = 1. / sqrt_e;
			const Type c_sqrt_e = c * sqrt_e_1, d_e = d * e_1;
			if (int root_num_U = PolySolveQuartic(d_e, c * e_1, 0., e_1, u_roots); root_num_U == 0) { // q2 < 0., e > 0.
				const Type h1 = boost::math::tools::halley_iterate(
					// Equation (21) in Dexter & Agol (2009) ApJ 696
					[c_sqrt_e, cubic_coeff = 2. * c_sqrt_e - sqrt_e * Power2(d_e)](Type h1) -> std::tuple<Type, Type, Type> {
						const Type h12 = Power2(h1), h14 = Power2(h12);
						return {(h12 - 1.) * (h14 - 1.) - h1 * (c_sqrt_e * (h14 + 1.) - cubic_coeff * h12),
								h1 * (6. * h14 - 4. * h12 - 2.) - c_sqrt_e * (5. * h14 + 1.) + 3. * cubic_coeff * h12,
								30. * h14 - (c_sqrt_e * 20. * h1 + 12.) * h12 + 6. * cubic_coeff * h1 - 2.};
					},
					// f(0) = 1, f(1) = -sqrt(e) * (d/e)^2 < 0
					0.5, 0.0, 1.0, std::numeric_limits<Type>::digits);
				const Type h2 = 1. / h1, g1 = d_e / (h2 - h1);
				I_u_0 = sqrt_e_1 * EllipticIntegral4Complex(0, u_obs, u_obj, 1., 0., sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				I_u_plus_0 = sqrt_e_1 * -EllipticIntegral4Complex(-2, u_obs, u_obj, 1., -u_plus_1, sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				I_u_minus_0 = sqrt_e_1 * -EllipticIntegral4Complex(-2, u_obs, u_obj, 1., -u_minus_1, sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				I_u_2_0 = sqrt_e_1 * EllipticIntegral4Complex(-2, u_obs, u_obj, 0., 1., sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				I_u_4_0 = sqrt_e_1 * EllipticIntegral4Complex(-4, u_obs, u_obj, 0., 1., sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				if (isnan(I_u_0) || isnan(I_u_plus_0) || isnan(I_u_minus_0) || isnan(I_u_2_0) || isnan(I_u_4_0))
					return Status::NUMERIC_ERROR;
				return Status::SUCCESS;
			} else if (root_num_U == 2) {
				const Type f = e_1 / (u_roots[0] * u_roots[1]), g = (1. / u_roots[0] + 1. / u_roots[1]) * f; // 1. / (−a^2*q^2*u_1*u_4)
				if (q2 < 0.) {
					// root[0] < root[1] < 0.
					I_u_0 = sqrt_e_1 * EllipticIntegral2Complex(0, u_obs, u_obj, 1., 0., f, g, 1., -u_roots[0], 1., -u_roots[1], 1.);
					I_u_plus_0 = sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, 1., -u_plus_1, f, g, 1., -u_roots[0], 1., -u_roots[1], 1.);
					I_u_minus_0 = sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, 1., -u_minus_1, f, g, 1., -u_roots[0], 1., -u_roots[1], 1.);
					I_u_2_0 = sqrt_e_1 * EllipticIntegral2Complex(-2, u_obs, u_obj, 0., 1., f, g, 1., -u_roots[0], 1., -u_roots[1], 1.);
					I_u_4_0 = sqrt_e_1 * EllipticIntegral2Complex(-4, u_obs, u_obj, 0., 1., f, g, 1., -u_roots[0], 1., -u_roots[1], 1.);
					if (isnan(I_u_0) || isnan(I_u_plus_0) || isnan(I_u_minus_0) || isnan(I_u_2_0) || isnan(I_u_4_0))
						return Status::NUMERIC_ERROR;
					return Status::SUCCESS;
				}
				if (u_roots[1] < u_obj) {
					I_u_0 = u_roots[1] / u_obj;
					return Status::NO_INTERSECTION_WITH_TARGET_SURFACE;
				}
				if (u_roots[1] >= u_plus) { // photon falls into the BH
					I_u_0 = sqrt_e_1 * EllipticIntegral2Complex(0, u_obs, u_obj, 1., 0., f, g, 1., u_roots[1], -1., -u_roots[0], 1.);
					I_u_plus_0 = sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, 1., -u_plus_1, f, g, 1., u_roots[1], -1., -u_roots[0], 1.);
					I_u_minus_0 = sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, 1., -u_minus_1, f, g, 1., u_roots[1], -1., -u_roots[0], 1.);
					I_u_2_0 = sqrt_e_1 * EllipticIntegral2Complex(-2, u_obs, u_obj, 0., 1., f, g, 1., u_roots[1], -1., -u_roots[0], 1.);
					I_u_4_0 = sqrt_e_1 * EllipticIntegral2Complex(-4, u_obs, u_obj, 0., 1., f, g, 1., u_roots[1], -1., -u_roots[0], 1.);
					if (isnan(I_u_0) || isnan(I_u_plus_0) || isnan(I_u_minus_0) || isnan(I_u_2_0) || isnan(I_u_4_0))
						return Status::NUMERIC_ERROR;
					return Status::SUCCESS;
				}
				AMinusPlusB(sqrt_e_1 * EllipticIntegral2Complex(0, u_obs, u_roots[1], 1., 0., f, g, 1., u_roots[1], -1., -u_roots[0], 1.), sqrt_e_1 * EllipticIntegral2Complex(0, u_obj, u_roots[1], 1., 0., f, g, 1., u_roots[1], -1., -u_roots[0], 1.), I_u_0, I_u_1);
				AMinusPlusB(sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_roots[1], 1., -u_plus_1, f, g, 1., u_roots[1], -1., -u_roots[0], 1.), sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obj, u_roots[1], 1., -u_plus_1, f, g, 1., u_roots[1], -1., -u_roots[0], 1.), I_u_plus_0, I_u_plus_1);
				AMinusPlusB(sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_roots[1], 1., -u_minus_1, f, g, 1., u_roots[1], -1., -u_roots[0], 1.), sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obj, u_roots[1], 1., -u_minus_1, f, g, 1., u_roots[1], -1., -u_roots[0], 1.), I_u_minus_0, I_u_minus_1);
				AMinusPlusB(sqrt_e_1 * EllipticIntegral2Complex(-2, u_obs, u_roots[1], 0., 1., f, g, 1., u_roots[1], -1., -u_roots[0], 1.), sqrt_e_1 * EllipticIntegral2Complex(-2, u_obj, u_roots[1], 0., 1., f, g, 1., u_roots[1], -1., -u_roots[0], 1.), I_u_2_0, I_u_2_1);
				AMinusPlusB(sqrt_e_1 * EllipticIntegral2Complex(-4, u_obs, u_roots[1], 0., 1., f, g, 1., u_roots[1], -1., -u_roots[0], 1.), sqrt_e_1 * EllipticIntegral2Complex(-4, u_obj, u_roots[1], 0., 1., f, g, 1., u_roots[1], -1., -u_roots[0], 1.), I_u_4_0, I_u_4_1);
				return Status::SUCCESS;
			}
			// root_num_U == 4
			if (u_roots[1] < u_obj) {
				I_u_0 = u_roots[1] / u_obj;
				return Status::NO_INTERSECTION_WITH_TARGET_SURFACE;
			}
			if (u_roots[1] >= 1.) { // photon falls into the BH
				I_u_0 = sqrt_e_1 * EllipticIntegral(0, u_obs, u_obj, 1., 0., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.);
				I_u_plus_0 = sqrt_e_1 * -EllipticIntegral(-2, u_obs, u_obj, 1., -u_plus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.);
				I_u_minus_0 = sqrt_e_1 * -EllipticIntegral(-2, u_obs, u_obj, 1., -u_minus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.);
				I_u_2_0 = sqrt_e_1 * EllipticIntegral(-2, u_obs, u_obj, 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.);
				I_u_4_0 = sqrt_e_1 * EllipticIntegral(-4, u_obs, u_obj, 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.);
				return Status::SUCCESS;
			}
			AMinusPlusB(sqrt_e_1 * EllipticIntegral(0, u_obs, u_roots[1], 1., 0., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.), sqrt_e_1 * EllipticIntegral(0, u_obj, u_roots[1], 1., 0., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.), I_u_0, I_u_1);
			AMinusPlusB(sqrt_e_1 * -EllipticIntegral(-2, u_obs, u_roots[1], 1., -u_plus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.), sqrt_e_1 * -EllipticIntegral(-2, u_obj, u_roots[1], 1., -u_plus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.), I_u_plus_0, I_u_plus_1);
			AMinusPlusB(sqrt_e_1 * -EllipticIntegral(-2, u_obs, u_roots[1], 1., -u_minus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.), sqrt_e_1 * -EllipticIntegral(-2, u_obj, u_roots[1], 1., -u_minus_1, u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.), I_u_minus_0, I_u_minus_1);
			AMinusPlusB(sqrt_e_1 * EllipticIntegral(-2, u_obs, u_roots[1], 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.), sqrt_e_1 * EllipticIntegral(-2, u_obj, u_roots[1], 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.), I_u_2_0, I_u_2_1);
			AMinusPlusB(sqrt_e_1 * EllipticIntegral(-4, u_obs, u_roots[1], 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.), sqrt_e_1 * EllipticIntegral(-4, u_obj, u_roots[1], 0., 1., u_roots[1], -1., -u_roots[0], 1., u_roots[2], -1., u_roots[3], -1.), I_u_4_0, I_u_4_1);
			return Status::SUCCESS;
		}
		static int Mu0Integral(Type a, Type q2, Type M_plus, Type M_minus, Type delta_M, Type delta_M_plus, int beta_sign, Type mu_obs, Type &mu_plus, Type &mu_minus, Type &A, Type &k, Type &n, int &alpha_1, Type &I_mu_0, Type &I_mu_full_turn, Type &I_t_mu_0, Type &I_t_mu_full_turn, Type &I_phi_mu_0, Type &I_phi_mu_full_turn) {
			if (std::abs(q2) < absolute_accuracy) { // M_minus == 0.
				mu_plus = std::copysign(SquareRoot(M_plus), mu_obs);
				mu_minus = 0.;
				alpha_1 = std::copysign(beta_sign, mu_obs);
				A = mu_plus * std::abs(a);
				k = 1.;
				const Type mu0_mu_plus = std::min(1., std::abs(mu_obs) / mu_plus);
				I_mu_full_turn = std::numeric_limits<Type>::infinity();
				I_t_mu_full_turn = std::numeric_limits<Type>::max();
				const Type phi = std::acos(mu0_mu_plus);
				I_t_mu_0 = SquareRoot(1. - Power2(mu0_mu_plus));
				I_mu_0 = std::log((1. + I_t_mu_0) / mu0_mu_plus) / A;
				if (M_plus == 1.)
					return Status::SUCCESS;
				n = M_plus / delta_M_plus;
				I_phi_mu_full_turn = std::numeric_limits<Type>::max();
				I_phi_mu_0 = boost::math::ellint_3(k, -n, phi); // TODO: Simplify
				return Status::SUCCESS;
			}
			if (M_minus > 0.) {
				mu_plus = std::copysign(SquareRoot(M_plus), mu_obs);
				mu_minus = std::copysign(SquareRoot(M_minus), mu_obs);
				alpha_1 = std::copysign(beta_sign, mu_obs);
				A = mu_plus * std::abs(a);
				k = SquareRoot(delta_M / M_plus);
				const Type x = std::min(1., SquareRoot((M_plus - Power2(mu_obs)) / delta_M));
				I_mu_full_turn = boost::math::ellint_1(k) / A;
				I_t_mu_full_turn = boost::math::ellint_2(k);
				const Type phi = asin(x);
				I_mu_0 = boost::math::ellint_1(k, phi) / A;
				I_t_mu_0 = boost::math::ellint_2(k, phi);
				if (M_plus == 1.)
					return Status::SUCCESS;
				n = delta_M / delta_M_plus;
				I_phi_mu_full_turn = boost::math::ellint_3(k, -n);
				I_phi_mu_0 = boost::math::ellint_3(k, -n, phi);
				return Status::SUCCESS;
			}
			// M_minus < 0.
			mu_plus = SquareRoot(M_plus);
			mu_minus = -mu_plus;
			alpha_1 = beta_sign;
			A = SquareRoot(delta_M) * std::abs(a);
			k = SquareRoot(M_plus / delta_M);
			// const Type x = SquareRoot(1. - Power2(mu_obs / mu_plus));
			I_mu_full_turn = 2. * boost::math::ellint_1(k) / A;
			I_t_mu_full_turn = 2. * boost::math::ellint_2(k);
			const Type phi = acos(std::min(1., std::abs(mu_obs) / mu_plus));
			if (mu_obs >= 0.) {
				I_mu_0 = boost::math::ellint_1(k, phi) / A;
				I_t_mu_0 = boost::math::ellint_2(k, phi);
				if (M_plus == 1.)
					return Status::SUCCESS;
				n = M_plus / delta_M_plus;
				I_phi_mu_full_turn = 2. * boost::math::ellint_3(k, -n);
				I_phi_mu_0 = boost::math::ellint_3(k, -n, phi);
				return Status::SUCCESS;
			}
			// the photon has to cross the equatorial plane.
			I_mu_0 = I_mu_full_turn - boost::math::ellint_1(k, phi) / A;
			I_t_mu_0 = I_t_mu_full_turn - boost::math::ellint_2(k, phi);
			if (M_plus == 1.)
				return Status::SUCCESS;
			n = M_plus / delta_M_plus;
			I_phi_mu_full_turn = 2. * boost::math::ellint_3(k, -n);
			I_phi_mu_0 = I_phi_mu_full_turn - boost::math::ellint_3(k, -n, phi);
			return Status::SUCCESS;
		}
		static int MuFIntegral(Type a, Type a2, Type l, Type q2, Type M_plus, Type M_minus, Type delta_M, Type delta_M_plus, Type mu_plus, Type mu_minus, Type I_u, Type I_mu_0, Type I_mu_full_turn, Type I_t_mu_0, Type I_t_mu_full_turn, Type I_phi_mu_0, Type I_phi_mu_full_turn, Type A, Type k, Type n, int alpha_1, int &alpha_2, Type &mu_f, Type &t_mu, Type &phi_mu) {
			const int N = alpha_1 > 0 ? ceil((I_u - I_mu_0) / I_mu_full_turn) : floor((I_u + I_mu_0) / I_mu_full_turn); // s_mu = alpha_1
			int alpha_3;
			// Here we follow the calculation in the code GEOKERR.
			// alpha_2 = (N & 1) ? alpha_1 : -alpha_1;
			// alpha_3 = 2 * floor((2 * N + 1 - alpha_1) / 4);
			if (N & 1) {
				alpha_2 = alpha_1;
				if (alpha_1 == 1)
					alpha_3 = N - 1;
				else
					alpha_3 = N + 1;
			} else {
				alpha_2 = -alpha_1;
				alpha_3 = N;
			}
			Type phi;
			if (std::abs(q2) < absolute_accuracy) {
				// assert alpha_3 == 0
				const Type cos_phi = 1. / cosh(std::abs(a * mu_plus) * (I_u - alpha_1 * I_mu_0));
				const Type sin_phi = tanh(std::abs(a * mu_plus) * (I_u - alpha_1 * I_mu_0));
				mu_f = mu_plus * cos_phi;
				if (cos_phi < sin_phi)
					phi = acos(cos_phi);
				else
					phi = asin(sin_phi);
				t_mu = A * (alpha_1 * I_t_mu_0 + alpha_2 * sin_phi + alpha_3 * I_t_mu_full_turn);
				if (M_plus == 1.) {
					if (alpha_1 != alpha_2)
						phi_mu = 0.;
					else
						phi_mu = std::copysign(boost::math::constants::pi<Type>(), l);
				} else
					phi_mu = l * (-I_u + (alpha_1 * I_phi_mu_0 + alpha_2 * boost::math::ellint_3(k, -n, phi) + alpha_3 * I_phi_mu_full_turn) / (A * delta_M_plus)); // TODO: simplify;
				return Status::SUCCESS;
			}
			const Type I = (I_u - alpha_1 * I_mu_0 - alpha_3 * I_mu_full_turn) / alpha_2;
			if (M_minus > 0.) {
				mu_f = mu_plus * boost::math::jacobi_dn(k, mu_plus * std::abs(a) * I);
				if (Type sin2_phi = (M_plus - Power2(mu_f)) / delta_M; sin2_phi < 0.5)
					phi = asin(SquareRoot(sin2_phi));
				else
					phi = acos(SquareRoot((Power2(mu_f) - M_minus) / delta_M));
				t_mu = A * (alpha_1 * I_t_mu_0 + alpha_2 * boost::math::ellint_2(k, phi) + alpha_3 * I_t_mu_full_turn);
				if (M_plus == 1.) {
					if (alpha_1 != alpha_2)
						phi_mu = 0.;
					else
						phi_mu = std::copysign(boost::math::constants::pi<Type>(), l);
				} else
					phi_mu = l * (-I_u + (alpha_1 * I_phi_mu_0 + alpha_2 * boost::math::ellint_3(k, -n, phi) + alpha_3 * I_phi_mu_full_turn) / (A * delta_M_plus));
				return Status::SUCCESS;
			}
			Type J;
			Type sn, cn, dn;
			// M_minus < 0.
			if (I <= 0.5 * I_mu_full_turn) {
				sn = boost::math::jacobi_elliptic(k, SquareRoot(delta_M) * std::abs(a) * I, &cn, &dn);
				mu_f = mu_plus * cn;
				if (cn < sn)
					phi = std::acos(cn);
				else
					phi = std::asin(sn);
				t_mu = a2 * M_minus * I_u + A * (alpha_1 * I_t_mu_0 + alpha_2 * boost::math::ellint_2(k, phi) + alpha_3 * I_t_mu_full_turn);
				if (M_plus == 1.) {
					if (alpha_1 != alpha_2)
						phi_mu = 0.;
					else
						phi_mu = std::copysign(boost::math::constants::pi<Type>(), l);
				} else
					phi_mu = l * (-I_u + (alpha_1 * I_phi_mu_0 + alpha_2 * boost::math::ellint_3(k, -n, phi) + alpha_3 * I_phi_mu_full_turn) / (A * delta_M_plus));
				return Status::SUCCESS;
			}
			// I > 0.5 * I_mu_full_turn, the photon has to cross the equatorial plane.
			sn = boost::math::jacobi_elliptic(k, SquareRoot(delta_M) * std::abs(a) * (I_mu_full_turn - I), &cn, &dn);
			mu_f = mu_minus * cn;
			if (cn < sn)
				phi = std::acos(cn);
			else
				phi = std::asin(sn);
			t_mu = a2 * M_minus * I_u + A * (alpha_1 * I_t_mu_0 - alpha_2 * boost::math::ellint_2(k, phi) + (alpha_2 + alpha_3) * I_t_mu_full_turn);
			if (M_plus == 1.) {
				if (alpha_1 != alpha_2)
					phi_mu = 0.;
				else
					phi_mu = std::copysign(boost::math::constants::pi<Type>(), l);
			} else
				phi_mu = l * (-I_u + (alpha_1 * I_phi_mu_0 - alpha_2 * boost::math::ellint_3(k, -n, phi) + (alpha_2 + alpha_3) * I_phi_mu_full_turn) / (A * delta_M_plus));
			return Status::SUCCESS;
		}

	  public:
		const Type a_, a2_, a4_;
		const Type sqrt_delta_a2_;
		const Type u_plus_1, u_plus, u_minus_1, u_minus, u_r;
		Kerr(Type spin) : a_(spin), a2_(a_ * a_), a4_(a2_ * a2_), sqrt_delta_a2_(std::sqrt(1. - a2_)), u_plus_1(1. + sqrt_delta_a2_), u_plus(1. / u_plus_1), u_minus_1(a2_ * u_plus), u_minus(1. / u_minus_1), u_r(-0.5 / sqrt_delta_a2_) {}
		std::string Name() const override {
			return "Kerr";
		}
		int MetricTensor(const boost::numeric::ublas::bounded_vector<Type, 8> &position, boost::numeric::ublas::bounded_matrix<Type, 4, 4> &metric) override {
			const Type r = position[1], r2 = r * r, Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * Power2(std::cos(position[2])), sin2_theta = Power2(std::sin(position[2]));
			if (rho2 == 0. || Delta == 0.)
				return Status::DIVISION_BY_ZERO;
			const Type r_rho_2 = 2. * r / rho2;
			metric.clear();
			metric(0, 0) = -(1. - r_rho_2);
			metric(0, 3) = -r_rho_2 * a_ * sin2_theta;
			metric(1, 1) = rho2 / Delta;
			metric(2, 2) = rho2;
			metric(3, 0) = -r_rho_2 * a_ * sin2_theta;
			metric(3, 3) = (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * sin2_theta;
			return Status::SUCCESS;
		}
		Type DistanceSquare(const boost::numeric::ublas::bounded_vector<Type, 8> &x, const boost::numeric::ublas::bounded_vector<Type, 8> &y, bool include_time) override {
			const Type r = x[1], r2 = r * r, d0 = x[0] - y[0], d3 = PhiDifference(x[3] - y[3]);
			const Type Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * Power2(std::cos(x[2])), r_rho_2 = 2. * r / rho2, sin2_theta = Power2(std::sin(x[2]));
			if (include_time)
				return (r_rho_2 - 1.) * Power2(d0) - 2. * r_rho_2 * a_ * sin2_theta * d0 * d3 + (rho2 / Delta) * Power2(x[1] - y[1]) + rho2 * Power2(x[2] - y[2]) + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * Power2(sin2_theta)) * Power2(d3);
			return (rho2 / Delta) * Power2(x[1] - y[1]) + rho2 * Power2(x[2] - y[2]) + ((r2 + a2_) * sin2_theta + 2. * r * a2_ * Power2(sin2_theta) / rho2) * Power2(d3);
		}
		int LagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type dt_dtau = 1. / y[4], r2 = Power2(y[1]), sin2_theta = Power2(std::sin(y[2]));
			const Type Delta = r2 - 2. * y[1] + a2_, rho2 = r2 + a2_ * Power2(std::cos(y[2]));
			if (rho2 == 0. || Delta == 0.)
				return Status::DIVISION_BY_ZERO;
			const Type r_rho_2 = 2. * y[1] / rho2;
			y[4] = (y[4] - 1. + r_rho_2 * (1. - a_ * sin2_theta * y[7])) * dt_dtau; // 1 + p_t
			y[5] *= rho2 / Delta * dt_dtau;
			y[6] *= rho2 * dt_dtau;
			y[7] = (-r_rho_2 * a_ + (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * y[7]) * sin2_theta * dt_dtau;
			return Status::SUCCESS;
		}
		int HamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type pt = y[4] - 1., r2 = Power2(y[1]), sin2_theta = Power2(std::sin(y[2]));
			const Type Delta = r2 - 2. * y[1] + a2_, rho2 = r2 + a2_ * Power2(std::cos(y[2]));
			if (rho2 == 0. || Delta == 0.)
				return Status::DIVISION_BY_ZERO;
			const Type rho_2 = 1. / rho2, r_rho_2 = 2. * y[1] * rho_2;
			y[4] = -Delta / ((Delta + r_rho_2 * (a2_ + r2)) * pt + r_rho_2 * a_ * y[7]);
			y[5] *= Delta * y[4] * rho_2;
			y[6] *= y[4] * rho_2;
			if (y[7] == 0.)
				y[7] = -r_rho_2 * a_ * pt / Delta * y[4];
			else if (sin2_theta == 0.)
				return Status::DIVISION_BY_ZERO;
			else
				y[7] = (-r_rho_2 * a_ * pt + (1. - r_rho_2) / sin2_theta * y[7]) / Delta * y[4];
			return Status::SUCCESS;
		}
		int FastTrace(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type r_object, const Type theta_object, const Type phi_object, Type &alpha, Type &beta, boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time) override {
			namespace ublas = boost::numeric::ublas;
			const Type sin_theta_object = SinTheta(theta_object), cos_theta_object = CosTheta(theta_object), sin_phi_object = std::sin(phi_object), cos_phi_object = std::cos(phi_object);
			const Type sin2_theta_object = Power2(sin_theta_object), cos2_theta_object = Power2(cos_theta_object);
			const Type cos_observer_object = sin_theta_observer * sin_theta_object * cos_phi_object + cos_theta_observer * cos_theta_object;
			const Type theta_observer_object = acos(cos_observer_object), sin_observer_object = SquareRoot(1. - Power2(cos_observer_object));
			const Type r2_object = Power2(r_object), r2_observer = Power2(r_observer), u_observer = 1. / r_observer, u_object = 1. / r_object;
			const Type Delta_object = r_object * (r_object - 2.) + a2_, rho2_object = r2_object + a2_ * cos2_theta_object;
			const Type rho_2_object = 1. / rho2_object;
			const Type r_rho_2_object = 2. * r_object * rho_2_object;
			const Type g00_object = -(1. - r_rho_2_object), g03_object = -r_rho_2_object * a_ * sin2_theta_object;
			const Type g11_object = rho2_object / Delta_object, g22_object = rho2_object, g33_object = (r2_object + a2_) * sin2_theta_object - g03_object * a_ * sin2_theta_object;
			KerrFastTraceParameters<Type> fast_trace_parameters(this, r_observer, r2_observer, u_observer, u_object, cos_theta_observer, cos_theta_object, theta_observer, sin_theta_observer, sin_theta_object, PhiDifference(phi_object) > -boost::math::constants::half_pi<Type>() ? PhiDifference(phi_object) : ModBy2Pi(phi_object));
			ublas::bounded_vector<Type, 2> alpha_beta_initial_value;
			const Type effective_radius = r_object + Power3(theta_observer_object / boost::math::constants::half_pi<Type>()) / sin_observer_object;
			const Type alpha_coefficient = sin_theta_object * sin_phi_object, beta_coefficient = cos_theta_object * sin_theta_observer - sin_theta_object * cos_phi_object * cos_theta_observer;
			alpha_beta_initial_value(0) = effective_radius * alpha_coefficient;
			alpha_beta_initial_value(1) = effective_radius * beta_coefficient;
			DNewtonRotationMultiFunctionSolver<Type, 2, KerrFastTraceParameters<Type>> alpha_beta_rotation_solver(DeltaUMuPhi, fast_trace_parameters);
			int status;
			if (status = alpha_beta_rotation_solver.Set(alpha_beta_initial_value, theta_observer, sin_theta_observer, cos_theta_observer, r_object, sin_theta_object, cos_theta_object, phi_object, sin_phi_object, cos_phi_object, false); status == Status::SUCCESS)
				if (status = alpha_beta_rotation_solver.Solve(boost::math::tools::root_epsilon<Type>()); status == Status::SUCCESS) {
					alpha = alpha_beta_rotation_solver.Root()(0);
					beta = alpha_beta_rotation_solver.Root()(1);
					photon(0) = fast_trace_parameters.tau; // not important
					photon(1) = r_object;
					photon(2) = theta_object;
					photon(3) = phi_object;
					photon(4) = (Power2(g03_object) - g00_object * g33_object) / (fast_trace_parameters.E * g33_object + fast_trace_parameters.L * g03_object);
					photon(7) = -(g00_object * fast_trace_parameters.L + g03_object * fast_trace_parameters.E) / (g03_object * fast_trace_parameters.L + g33_object * fast_trace_parameters.E);
					photon(6) = -fast_trace_parameters.mu_dir * std::sqrt(fast_trace_parameters.Q * Power2(photon[4]) - cos2_theta_object * (Power2(-r_rho_2_object * a_ + (a2_ + r2_object + r_rho_2_object * a2_ * sin2_theta_object) * photon[7]) * sin2_theta_object - a2_ * Power2(r_rho_2_object * (1. - a_ * sin2_theta_object * photon[7]) - 1.))) * rho_2_object;
					photon(5) = -fast_trace_parameters.u_dir * std::sqrt(-(g00_object + g22_object * Power2(photon[6]) + 2. * g03_object * photon[7] + g33_object * Power2(photon[7])) / g11_object); // solved by normalization
					photon_time = fast_trace_parameters.t;
					return Status::SUCCESS;
				}
			// PrintlnWarning("Kerr FastTrace() ROTATION failed with status = {}", status);
			DNewtonTranslationMultiFunctionSolver<Type, 2, KerrFastTraceParameters<Type>> alpha_beta_translation_solver(DeltaUMuPhi, fast_trace_parameters);
			if (status = alpha_beta_translation_solver.Set(alpha_beta_rotation_solver.Root(), theta_observer, sin_theta_observer, cos_theta_observer, r_object, sin_theta_object, cos_theta_object, phi_object, sin_phi_object, cos_phi_object, false); status == Status::SUCCESS)
				if (status = alpha_beta_translation_solver.Solve(boost::math::tools::root_epsilon<Type>()); status == Status::SUCCESS) {
					alpha = alpha_beta_translation_solver.Root()(0);
					beta = alpha_beta_translation_solver.Root()(1);
					photon(0) = fast_trace_parameters.tau; // not important
					photon(1) = r_object;
					photon(2) = theta_object;
					photon(3) = phi_object;
					photon(4) = (Power2(g03_object) - g00_object * g33_object) / (fast_trace_parameters.E * g33_object + fast_trace_parameters.L * g03_object);
					photon(7) = -(g00_object * fast_trace_parameters.L + g03_object * fast_trace_parameters.E) / (g03_object * fast_trace_parameters.L + g33_object * fast_trace_parameters.E);
					photon(6) = -fast_trace_parameters.mu_dir * std::sqrt(fast_trace_parameters.Q * Power2(photon[4]) - cos2_theta_object * (Power2(-r_rho_2_object * a_ + (a2_ + r2_object + r_rho_2_object * a2_ * sin2_theta_object) * photon[7]) * sin2_theta_object - a2_ * Power2(r_rho_2_object * (1. - a_ * sin2_theta_object * photon[7]) - 1.))) * rho_2_object;
					photon(5) = -fast_trace_parameters.u_dir * std::sqrt(-(g00_object + g22_object * Power2(photon[6]) + 2. * g03_object * photon[7] + g33_object * Power2(photon[7])) / g11_object); // solved by normalization
					photon_time = fast_trace_parameters.t;
					return Status::SUCCESS;
				}
			DirectionMultiFunctionSolver<Type, 2, KerrFastTraceParameters<Type>> alpha_beta_direction_solver(DeltaUMuPhi, fast_trace_parameters);
			if (status = alpha_beta_direction_solver.Set(alpha_beta_translation_solver.Root()); status != Status::SUCCESS)
				// PrintlnWarning("Kerr FastTrace() set DIRECTION failed with status = {}", status);
				if (status = alpha_beta_direction_solver.Solve(boost::math::tools::root_epsilon<Type>(), 2048); status == Status::SUCCESS) {
					// PrintlnWarning("Kerr FastTrace() DIRECTION failed with status = {}", status);
					alpha = alpha_beta_direction_solver.Root()(0);
					beta = alpha_beta_direction_solver.Root()(1);
					photon(0) = fast_trace_parameters.tau; // not important
					photon(1) = r_object;
					photon(2) = theta_object;
					photon(3) = phi_object;
					photon(4) = (Power2(g03_object) - g00_object * g33_object) / (fast_trace_parameters.E * g33_object + fast_trace_parameters.L * g03_object);
					photon(7) = -(g00_object * fast_trace_parameters.L + g03_object * fast_trace_parameters.E) / (g03_object * fast_trace_parameters.L + g33_object * fast_trace_parameters.E);
					photon(6) = -fast_trace_parameters.mu_dir * std::sqrt(fast_trace_parameters.Q * Power2(photon[4]) - cos2_theta_object * (Power2(-r_rho_2_object * a_ + (a2_ + r2_object + r_rho_2_object * a2_ * sin2_theta_object) * photon[7]) * sin2_theta_object - a2_ * Power2(r_rho_2_object * (1. - a_ * sin2_theta_object * photon[7]) - 1.))) * rho_2_object;
					photon(5) = -fast_trace_parameters.u_dir * std::sqrt(-(g00_object + g22_object * Power2(photon[6]) + 2. * g03_object * photon[7] + g33_object * Power2(photon[7])) / g11_object); // solved by normalization
					photon_time = fast_trace_parameters.t;
					return Status::SUCCESS;
				}
			return status;
		}
		template <std::size_t N>
		int CalcThetaPhi(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, Type alpha, Type beta, const boost::numeric::ublas::bounded_vector<Type, N> &u, Type theta_0[], Type theta_1[], Type phi_0[], Type phi_1[]) {
			namespace ublas = boost::numeric::ublas;
			// The photon in the observer's frame has the tetrad velocity: [1, r / R, beta / R, -alpha / R], where R = std::sqrt(r^2 + alpha^2 + beta^2).
			ublas::bounded_vector<Type, 8> photon;
			Type photon_time;
			if (int status = this->InitializePhoton(photon, photon_time, alpha, beta, r_observer, Power2(r_observer), theta_observer, sin_theta_observer); status != Status::SUCCESS)
				return status;
			const Type e_1 = 1. / Energy(photon, T, LAGRANGIAN);
			const Type l = AngularMomentum(photon, T, LAGRANGIAN) * e_1, l2 = Power2(l);
			const Type q2 = CarterConstant(photon, 0., T, LAGRANGIAN) * Power2(e_1);
			const Type u_obs = 1. / r_observer;
			const Type mu_obs = cos_theta_observer;
			// U=1+[a^2−q^2−l^2]u^2+2[(a−l)^2+q^2]u^3−a^2q^2u^4
			ublas::bounded_vector<Type, N> I_u_0, I_u_1;
			ublas::bounded_vector<Type, N> int_u_plus_0, int_u_plus_1;
			ublas::bounded_vector<Type, N> int_u_minus_0, int_u_minus_1;
			ublas::bounded_vector<Type, N> int_u_2_0, int_u_2_1;
			ublas::bounded_vector<Type, N> int_u_4_0, int_u_4_1;
			Type A_0, A_1, B_0, B_1, C_0, C_1, D_0, D_1, E_0, E_1;
			for (int i = 0; i < N; ++i) {
				A_1 = 0.;
				if (int status = UIntegral(a_, a2_, u_plus_1, u_minus_1, u_plus, u_minus, l, l2, q2, r_observer, u_obs, u(i), A_0, A_1, B_0, B_1, C_0, C_1, D_0, D_1, E_0, E_1); status != Status::SUCCESS)
					break;
				I_u_0(i) = A_0;
				int_u_plus_0(i) = B_0;
				int_u_minus_0(i) = C_0;
				int_u_2_0(i) = D_0;
				int_u_4_0(i) = E_0;
				if (A_1 != 0.) {
					I_u_1(i) = A_1;
					int_u_plus_1(i) = B_1;
					int_u_minus_1(i) = C_1;
					int_u_2_1(i) = D_1;
					int_u_4_1(i) = E_1;
				}
			}
			// M = q^2 + (a^2-q^2-l^2)*mu^2-a^2*mu^4
			Type M_minus_plus[2];
			Type delta_M, delta_M_plus;
			if (int root_num_M = PolySolveQuadratic(-a2_, a2_ - l2 - q2, q2, M_minus_plus); root_num_M == 0)
				return Status::FAILURE;
			delta_M = std::sqrt(Power2(a2_ - l2 - q2) + 4. * a2_ * q2) / a2_;
			if (M_minus_plus[1] > 1.) {
				M_minus_plus[1] = 1.;
				delta_M_plus = 0.;
			} else if (M_minus_plus[1] > 0.9)
				delta_M_plus = 2. * l2 / (a2_ + l2 + q2 + std::sqrt(Power2(a2_ + l2 + q2) - 4. * a2_ * l2));
			else
				delta_M_plus = 1. - M_minus_plus[1];
			Type mu_plus, mu_minus;
			Type A, k, n;
			Type I_mu_0, I_mu_full_turn;
			Type I_t_mu_0, I_t_mu_full_turn;
			Type I_phi_mu_0, I_phi_mu_full_turn;
			// Here the first turning point is determined by `beta` of the observer.
			int alpha_1, alpha_2;
			Mu0Integral(a_, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, boost::math::sign(beta), mu_obs, mu_plus, mu_minus, A, k, n, alpha_1, I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn);
			Type mu_f_0, mu_f_1, t_mu, phi_mu_0, phi_mu_1;
			for (int i = 0; i < I_u_0.size(); ++i) {
				if (int status = MuFIntegral(a_, a2_, l, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, mu_plus, mu_minus, I_u_0[i], I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn, A, k, n, alpha_1, alpha_2, mu_f_0, t_mu, phi_mu_0); status != Status::SUCCESS)
					return status;
				theta_0[i] = acos(mu_f_0);
				const Type phi_u_0 = u_r * ((l * u_plus_1 + 2. * (a_ - l)) * int_u_plus_0[i] - (l * u_minus_1 + 2. * (a_ - l)) * int_u_minus_0[i]);
				phi_0[i] = phi_u_0 + phi_mu_0;
				if (isnan(phi_0[i]))
					return Status::FAILURE;
			}
			for (int i = 0; i < I_u_1.size(); ++i) {
				if (int status = MuFIntegral(a_, a2_, l, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, mu_plus, mu_minus, I_u_1[i], I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn, A, k, n, alpha_1, alpha_2, mu_f_1, t_mu, phi_mu_1); status != Status::SUCCESS)
					return status;
				theta_1[i] = acos(mu_f_1);
				const Type phi_u_1 = u_r * ((l * u_plus_1 + 2. * (a_ - l)) * int_u_plus_1[i] - (l * u_minus_1 + 2. * (a_ - l)) * int_u_minus_1[i]);
				phi_1[i] = phi_u_1 + phi_mu_1;
				if (isnan(phi_1[i]))
					return Status::FAILURE;
			}
			return Status::SUCCESS;
		}
		int FastShadow(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type alpha, const Type beta, const Type r_min) override {
			boost::numeric::ublas::bounded_vector<Type, 8> photon;
			Type photon_time;
			if (int status = this->InitializePhoton(photon, photon_time, alpha, beta, r_observer, Power2(r_observer), theta_observer, sin_theta_observer); status != Status::SUCCESS)
				return status;
			const Type E = Energy(photon, T, LAGRANGIAN), E_1 = 1. / E;
			const Type L = AngularMomentum(photon, T, LAGRANGIAN), l = L * E_1, l2 = Power2(l);
			const Type Q = CarterConstant(photon, 0., T, LAGRANGIAN), q2 = Q * Power2(E_1);
			const Type c = a2_ - l2 - q2, d = 2. * (Power2(a_ - l) + q2), e = -a2_ * q2;
			Type u_roots[4];
			if (std::abs(e) < absolute_accuracy) {
				// U=1+[a^2−l^2]u^2+2(a−l)^2u^3
				// This means the observer and the object are on the equatorial plane.
				if (a_ == l) // U=1.
					return 0;
				const Type d_1 = 1. / d, sqrt_d_1 = boost::math::constants::half_root_two<Type>() / std::abs(a_ - l);
				if (int root_num_U = PolySolveCubic(c * d_1, 0., d_1, u_roots); root_num_U == 1) {
					return 0;
				}
				if (u_roots[1] >= 1.) { // photon falls into the BH
					return 0;
				}
				return 1;
			}
			const Type e_1 = 1. / e, sqrt_e = std::sqrt(std::abs(e)), sqrt_e_1 = 1. / sqrt_e;
			const Type c_sqrt_e = c * sqrt_e_1, d_e = d * e_1;
			if (int root_num_U = PolySolveQuartic(d_e, c * e_1, 0., e_1, u_roots); root_num_U == 0) { // q2 < 0., e > 0.
				return 0;
			} else if (root_num_U == 2) {
				const Type f = e_1 / (u_roots[0] * u_roots[1]), g = (1. / u_roots[0] + 1. / u_roots[1]) * f; // 1. / (−a^2*q^2*u_1*u_4)
				if (q2 < 0.) {
					// root[0] < root[1] < 0.
					return 0;
				}
				if (u_roots[1] >= 1.) { // photon falls into the BH
					return 0;
				}
				return 1;
			}
			// root_num_U == 4
			if (u_roots[1] >= 1.) { // photon falls into the BH
				return 0;
			}
			return 1;
		}
		Type Energy(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				if (time == T)
					return (1. - 2. * y[1] / (Power2(y[1]) + a2_ * Power2(std::cos(y[2]))) * (1. - a_ * Power2(std::sin(y[2])) * y[7])) / y[4];
				// time == TAU
				return y[4] - 2. * y[1] / (Power2(y[1]) + a2_ * Power2(std::cos(y[2]))) * (y[4] - a_ * Power2(std::sin(y[2])) * y[7]);
			} // dynamics == HAMILTONIAN
			return 1. - y[4];
		}
		Type AngularMomentum(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				const Type r2 = Power2(y[1]), sin2_theta = Power2(std::sin(y[2]));
				const Type r_rho_2 = 2. * y[1] / (r2 + a2_ * Power2(std::cos(y[2])));
				if (time == T)
					return (-r_rho_2 * a_ + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta / y[4];
				// time == TAU
				return (-r_rho_2 * a_ * y[4] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta;
			} // dynamics == HAMILTONIAN
			return y[7];
		}
		Type CarterConstant(const boost::numeric::ublas::bounded_vector<Type, 8> &y, const Type mu2, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				const Type r2 = Power2(y[1]), sin2_theta = Power2(std::sin(y[2])), cos2_theta = Power2(std::cos(y[2]));
				const Type rho2 = r2 + a2_ * cos2_theta;
				const Type r_rho_2 = 2. * y[1] / rho2;
				if (time == T)
					return mu2 * cos2_theta * a2_ + (Power2(rho2 * y[6]) + cos2_theta * (Power2(-r_rho_2 * a_ + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta - a2_ * Power2(r_rho_2 * (1. - a_ * sin2_theta * y[7]) - 1.))) / Power2(y[4]);
				// time == TAU
				return mu2 * cos2_theta * a2_ + Power2(rho2 * y[6]) + cos2_theta * (Power2(-r_rho_2 * a_ * y[4] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta - a2_ * Power2(r_rho_2 * (y[4] - a_ * sin2_theta * y[7]) - y[4]));
			} // dynamics == HAMILTONIAN
			return Power2(y[6]) + Power2(std::cos(y[2])) * (a2_ * (mu2 - Power2(1. - y[4])) + Power2(y[7] / std::sin(y[2])));
		}
		int NormalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type r = y[1], r2 = Power2(r), a2_r2 = a2_ + r2;
			const Type sin2_theta = Power2(std::sin(y[2])), sin4_theta = Power2(sin2_theta);
			const Type rho2 = r2 + a2_ * Power2(std::cos(y[2]));
			const Type r_rho_2 = 2. * r / rho2;
			// y[7] += 2. * a_ * r / (Power2(a2_r2) - a2_ * Delta * sin2_theta);
			y[4] = std::sqrt(1. - r_rho_2 + 2. * r_rho_2 * a_ * sin2_theta * y[7] - (rho2 / (r2 - 2. * r + a2_) * Power2(y[5]) + rho2 * Power2(y[6]) + (a2_r2 * sin2_theta + r_rho_2 * a2_ * sin4_theta) * Power2(y[7])));
			return isnan(y[4]) ? Status::DOMAIN_ERROR : Status::SUCCESS;
		}
		int NormalizeNullGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y, Type frequency = 1.) override {
			const Type r = y[1], r2 = r * r;
			const Type sin2_theta = Power2(std::sin(y[2]));
			const Type rho2 = r2 + a2_ * Power2(std::cos(y[2]));
			const Type r_rho_2 = 2. * r / rho2;
			const Type a = rho2 / (r2 - 2. * r + a2_) * y[5] * y[5] + rho2 * y[6] * y[6] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * sin2_theta * y[7] * y[7];
			const Type b = -2. * r_rho_2 * a_ * sin2_theta * y[7];
			const Type c = r_rho_2 - 1.;
			const Type coefficient = std::copysign(0.5 / a, frequency) * (-b + std::sqrt(b * b - 4. * a * c));
			if (isnan(coefficient))
				return Status::DOMAIN_ERROR;
			y[4] = frequency;
			y[5] *= coefficient;
			y[6] *= coefficient;
			y[7] *= coefficient;
			return Status::SUCCESS;
		}
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> GetIntegrationSystem(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override {
			if (time == T) {
				if (dynamics == LAGRANGIAN) {
					if (motion == GEODESIC)
						return [a = this->a_, a2 = this->a2_, a4 = this->a4_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = y[6]; // d\theta/dt
							dydt[3] = y[7]; // d\phi/dt
							const Type r = y[1], r2 = Power2(r), r4 = Power2(r2), a2_r2 = a2 + r2;
							const Type sin_theta = SinTheta(y[2]), sin2_theta = Power2(sin_theta), sin4_theta = Power2(sin2_theta), cos_theta = CosTheta(y[2]), cos2_theta = Power2(cos_theta), sin_theta_cos_theta = sin_theta * cos_theta, cot_theta = cos_theta / sin_theta;
							const Type Delta = r2 - 2. * r + a2, Delta_1 = 1. / Delta;
							const Type rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho4 = Power2(rho2), rho_4 = Power2(rho_2), rho_6 = rho_2 * rho_4, r2_a2_cos2_theta = r2 - a2 * cos2_theta;
							const Type dydt4 = 2. * rho_4 * (Delta_1 * a2_r2 * r2_a2_cos2_theta * y[5] - 2. * a2 * r * sin_theta_cos_theta * y[6] * (1. - a * sin2_theta * y[7]) - Delta_1 * a * (2. * r4 + r2 * rho2 + a2 * r2_a2_cos2_theta) * sin2_theta * y[5] * y[7]);
							// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
							dydt[4] = dydt4 * y[4];
							// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[5] = (-Delta * r2_a2_cos2_theta * Power2(1. - a * sin2_theta * y[7]) - (r * (a2 * sin2_theta - r) + a2 * cos2_theta) * Delta_1 * rho4 * Power2(y[5]) + 2. * a2 * sin_theta_cos_theta * rho4 * y[5] * y[6] + r * Delta * rho4 * Power2(y[6]) + Delta * sin2_theta * r * rho4 * Power2(y[7])) * rho_6 + dydt4 * y[5];
							// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[6] = (2. * a2 * r * sin_theta_cos_theta - a2 * sin_theta_cos_theta * rho4 * (Delta_1 * Power2(y[5]) - Power2(y[6])) - 2. * r * rho4 * y[5] * y[6] - 4. * a * r * sin_theta_cos_theta * a2_r2 * y[7] + sin_theta_cos_theta * (2. * a4 * r * sin4_theta + 4. * a2 * r * sin2_theta * rho2 + a2_r2 * rho4) * Power2(y[7])) * rho_6 + dydt4 * y[6];
							// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[7] = (-2. * a * r2_a2_cos2_theta * Delta_1 * y[5] + 4. * a * r * cot_theta * y[6] - 2. * Delta_1 * (r * rho4 - 2. * r2 * rho2 - r2_a2_cos2_theta * a2 * sin2_theta) * y[5] * y[7] - 2. * cot_theta * (rho4 + 2. * a2 * r * sin2_theta) * y[6] * y[7]) * rho_4 + dydt4 * y[7];
						};
					else if (motion == CIRCULAR)
						return [](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = 0.;	// dr/dt
							dydt[2] = 0.;	// d\theta/dt
							dydt[3] = y[7]; // d\phi/dt
							dydt[4] = 0.;
							dydt[5] = 0.;
							dydt[6] = 0.;
							dydt[7] = 0.;
						};
					else if (motion == HELICAL)
						return [a = this->a_, a2 = this->a2_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = 0.;	// d\theta/dt = 0.
							dydt[3] = y[7]; // d\phi/dt
							const Type r = y[1], r2 = Power2(r);
							const Type sin2_theta = Power2(std::sin(y[2]));
							const Type Delta_1 = 1. / (r2 - 2. * r + a2);
							const Type rho2 = r2 + a2 * Power2(std::cos(y[2])), rho_2 = 1. / rho2;
							const Type r_rho_2 = 2. * y[1] * rho_2;
							const Type g00 = -(1. - r_rho_2), g11 = rho2 * Delta_1, g03 = -r_rho_2 * a * sin2_theta, g33 = (r2 + a2 - g03 * a) * sin2_theta;
							const Type L_1 = y[4] / (g03 + g33 * y[7]), L_2 = Power2(L_1);
							const Type A = g33 * (g33 * L_2 + 1.), A_1 = 1. / A, B = 2. * g03 * (g33 * L_2 + 1.), C = g00 + g11 * Power2(y[5]) + Power2(g03) * L_2;
							const Type dg00_dr = (2. - 4. * r2 * rho_2) * rho_2, dg11_dr = 2. * (r - g11 * (r - 1.)) * Delta_1, dg03_dr = -dg00_dr * a * sin2_theta, dg33_dr = 2. * r * sin2_theta - dg03_dr * a * sin2_theta;
							const Type dA_dr = dg33_dr * (2. * g33 * L_2 + 1.), dB_dr = 2. * (dg03_dr * (g33 * L_2 + 1.) + g03 * dg33_dr * L_2), dC_dr = dg00_dr + dg11_dr * Power2(y[5]) + 2. * g03 * dg03_dr * L_2;
							// d^2r/dt^2 = 0.
							dydt[5] = 0.;
							// d^2\theta/dt^2 = 0.
							dydt[6] = 0.;
							// d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
							dydt[7] = y[5] * A_1 * (-y[7] * dA_dr + 0.5 * (-dB_dr + (B * dB_dr - 2. * (A * dC_dr + C * dA_dr)) / (2. * A * y[7] + B)));
							// d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt = d((g03 + g33 * d\phi/dt) / L)/dr * dr/dt
							dydt[4] = (y[5] * (dg03_dr + dg33_dr * y[7]) + g33 * dydt[7]) * L_1;
						};
				} else if (dynamics == HAMILTONIAN && motion == GEODESIC)
					return [a = this->a_, a2 = this->a2_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
						const Type r = y[1], r2 = y[1] * y[1], a2_r2 = a2 + r2, pr2 = y[5] * y[5], ptheta2 = y[6] * y[6];
						const Type E = 1. - y[4], E2 = E * E, delta_E2 = (2. - y[4]) * y[4], L2 = y[7] * y[7];
						const Type sin_theta = SinTheta(y[2]), sin2_theta = sin_theta * sin_theta, sin_2_theta = 1. / sin2_theta, sin_4_theta = sin_2_theta * sin_2_theta, cos_theta = CosTheta(y[2]), cos2_theta = cos_theta * cos_theta;
						const Type Delta = a2_r2 - 2. * r, Delta_1 = 1. / Delta, Delta_2 = Delta_1 * Delta_1;
						const Type rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho_4 = rho_2 * rho_2;
						const Type Q = ptheta2 + cos2_theta * (a2 * delta_E2 + L2 * sin_2_theta);
						const Type R = -a2_r2 * r2 * delta_E2 + E2 * 2. * r * a2 + 2. * r * r2 - Delta * Q - (r2 - 2. * r) * L2 - 4. * r * a * E * y[7]; // R = Power2(E * a2_r2 - a * y[7]) - Delta * (r2 + Power2(y[7] - a * E) + Q);
						//[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
						dydt[0] = rho2 * Delta / (a2_r2 * (E * a2_r2 - a * y[7]) + a * Delta * (y[7] - a * E * sin2_theta));			// d\tau/dt
						dydt[1] = Delta * rho_2 * y[5] * dydt[0];																		// dr/dt
						dydt[2] = rho_2 * y[6] * dydt[0];																				// d\theta/dt
						dydt[3] = (2. * E * a * r - y[7] * (a2 - Delta * (1. + cos2_theta * sin_2_theta))) * Delta_1 * rho_2 * dydt[0]; // d\phi/dt
						dydt[4] = 0.;
						dydt[5] = ((a2 * cos2_theta + r * a2 * sin2_theta - r2) * pr2 + ((-2. * delta_E2 * r2 - a2 * delta_E2 + 3. * r - L2 - Q) * r + a2 * E2 + L2 - 2. * a * E * y[7] + Q) * rho2 * Delta_1 - ((r - 1.) * rho2 + Delta * r) * R * Delta_2) * rho_4 * dydt[0];
						dydt[6] = (-(Delta * pr2 - R * Delta_1) * a2 * rho_2 + a2 * delta_E2 + L2 * sin_2_theta + cos2_theta * sin_4_theta * L2) * sin_theta * cos_theta * rho_2 * dydt[0];
						dydt[7] = 0.;
					};
			} else if (time == TAU) {
				if (dynamics == LAGRANGIAN)
					return [a = this->a_, a2 = this->a2_, a4 = this->a4_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
						dydt[0] = y[4]; // dt/d\tau
						dydt[1] = y[5]; // dr/d\tau
						dydt[2] = y[6]; // d\theta/d\tau
						dydt[3] = y[7]; // d\phi/d\tau
						const Type r = y[1], r2 = Power2(r), a2_r2 = a2 + r2;
						const Type sin_theta = SinTheta(y[2]), sin2_theta = Power2(sin_theta), sin4_theta = Power2(sin2_theta), cos_theta = CosTheta(y[2]), cos2_theta = Power2(cos_theta), sin_theta_cos_theta = sin_theta * cos_theta, cot_theta = cos_theta / sin_theta;
						const Type Delta = r2 - 2. * r + a2, Delta_1 = 1. / Delta;
						const Type rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho4 = Power2(rho2), rho_4 = Power2(rho_2), rho_6 = rho_2 * rho_4, r2_a2_cos2_theta = r2 - a2 * cos2_theta;
						dydt[4] = -2. * rho_4 * (Delta_1 * a2_r2 * r2_a2_cos2_theta * y[5] * (y[4] - a * sin2_theta * y[7]) - 2. * a2 * r * sin_theta_cos_theta * y[6] * (y[4] - a * sin2_theta * y[7]) - 2. * Delta_1 * r2 * rho2 * a * sin2_theta * y[5] * y[7]);
						// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
						dydt[5] = (-Delta * r2_a2_cos2_theta * Power2(y[4] - a * sin2_theta * y[7]) - (r * (a2 * sin2_theta - r) + a2 * cos2_theta) * Delta_1 * rho4 * Power2(y[5]) + 2. * a2 * sin_theta_cos_theta * rho4 * y[5] * y[6] + r * Delta * rho4 * Power2(y[6]) + Delta * sin2_theta * r * rho4 * Power2(y[7])) * rho_6;
						// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
						dydt[6] = (2. * a2 * r * sin_theta_cos_theta * Power2(y[4]) - a2 * sin_theta_cos_theta * rho4 * (Delta_1 * Power2(y[5]) - Power2(y[6])) - 2. * r * rho4 * y[5] * y[6] - 4. * a * r * sin_theta_cos_theta * a2_r2 * y[4] * y[7] + sin_theta_cos_theta * (2. * a4 * r * sin4_theta + 4. * a2 * r * sin2_theta * rho2 + a2_r2 * rho4) * Power2(y[7])) * rho_6;
						// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
						dydt[7] = (-2. * a * r2_a2_cos2_theta * Delta_1 * y[4] * y[5] + 4. * a * r * cot_theta * y[4] * y[6] - 2. * Delta_1 * (r * rho4 - 2. * r2 * rho2 - r2_a2_cos2_theta * a2 * sin2_theta) * y[5] * y[7] - 2. * cot_theta * (rho4 + 2. * a2 * r * sin2_theta) * y[6] * y[7]) * rho_4;
					};
				else if (dynamics == HAMILTONIAN)
					return [a = this->a_, a2 = this->a2_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
						const Type r = y[1], r2 = Power2(y[1]), a2_r2 = a2 + r2, pr2 = Power2(y[5]), ptheta2 = Power2(y[6]);
						const Type E = 1. - y[4], E2 = Power2(E), delta_E2 = (2. - y[4]) * y[4], L2 = Power2(y[7]);
						const Type sin_theta = SinTheta(y[2]), sin2_theta = Power2(sin_theta), sin_2_theta = 1. / sin2_theta, sin_4_theta = Power2(sin_2_theta), cos_theta = CosTheta(y[2]), cos2_theta = Power2(cos_theta);
						const Type Delta = a2_r2 - 2. * r, Delta_1 = 1. / Delta, Delta_2 = Power2(Delta_1);
						const Type rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho_4 = Power2(rho_2);
						const Type Q = ptheta2 + cos2_theta * (a2 * delta_E2 + L2 * sin_2_theta);
						const Type R = -a2_r2 * r2 * delta_E2 + E2 * 2. * r * a2 + 2. * r * r2 - Delta * Q - (r2 - 2. * r) * L2 - 4. * r * a * E * y[7]; // R = Power2(E * a2_r2 - a * y[7]) - Delta * (r2 + Power2(y[7] - a * E) + Q);
						//[t,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
						dydt[0] = (a2_r2 * (E * a2_r2 - a * y[7]) * Delta_1 + a * (y[7] - a * E * sin2_theta)) * rho_2;		  // dt/d\tau
						dydt[1] = Delta * rho_2 * y[5];																		  // dr/d\tau
						dydt[2] = rho_2 * y[6];																				  // d\theta/d\tau
						dydt[3] = (2. * E * a * r - y[7] * (a2 - Delta * (1. + cos2_theta * sin_2_theta))) * Delta_1 * rho_2; // d\phi/d\tau
						dydt[4] = 0.;
						dydt[5] = ((a2 * cos2_theta + r * a2 * sin2_theta - r2) * pr2 + ((-2. * delta_E2 * r2 - a2 * delta_E2 + 3. * r - L2 - Q) * r + a2 * E2 + L2 - 2. * a * E * y[7] + Q) * rho2 * Delta_1 - ((r - 1.) * rho2 + Delta * r) * R * Delta_2) * rho_4;
						dydt[6] = (-(Delta * pr2 - R * Delta_1) * a2 * rho_2 + a2 * delta_E2 + L2 * sin_2_theta + cos2_theta * sin_4_theta * L2) * sin_theta * cos_theta * rho_2;
						dydt[7] = 0.;
					};
			}
			throw std::invalid_argument(fmt::format("Kerr::GetIntegrationSystem({}, {}, {}) invaild", time, dynamics, motion));
		}
	};

	template <typename Type>
	struct KerrFastTraceParameters {
		Kerr<Type> *const kerr;
		const Type r, r2;
		const Type u_obs, u_obj;
		const Type mu_obs, mu_obj;
		const Type theta_obs;
		const Type sin_theta_obs, sin_theta_obj;
		const Type phi_obj;
		Type E, L, Q;
		Type t;
		Type tau;
		int u_dir, mu_dir;
		KerrFastTraceParameters(Kerr<Type> *kerr, Type r, Type r2, Type u_obs, Type u_obj, Type mu_obs, Type mu_obj, Type theta_obs, Type sin_theta_obs, Type sin_theta_obj, Type phi_obj) : kerr(kerr), r(r), r2(r2), u_obs(u_obs), u_obj(u_obj), mu_obs(mu_obs), mu_obj(mu_obj), theta_obs(theta_obs), sin_theta_obs(sin_theta_obs), sin_theta_obj(sin_theta_obj), phi_obj(phi_obj) {}
	};

	template <typename Type>
	class KerrNewman : public Metric<Type> {
	  protected:
		Type DotProduct3D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x1, Type x2, Type x3, Type y1, Type y2, Type y3) override {
			const Type r = position[1], r2 = Power2(r), Delta = r2 - 2. * r + a2_ + r_Q2_, rho2 = r2 + a2_ * Power2(std::cos(position[2])), r_rho_2 = (2. * r - r_Q2_) / rho2, sin2_theta = Power2(std::sin(position[2]));
			return (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * sin2_theta * sin2_theta) * x3 * y3;
		}
		Type DotProduct4D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x0, Type x1, Type x2, Type x3, Type y0, Type y1, Type y2, Type y3) override {
			const Type r = position[1], r2 = Power2(r), Delta = r2 - 2. * r + a2_ + r_Q2_, rho2 = r2 + a2_ * Power2(std::cos(position[2])), r_rho_2 = (2. * r - r_Q2_) / rho2, sin2_theta = Power2(std::sin(position[2]));
			return (r_rho_2 - 1.) * x0 * y0 - r_rho_2 * a_ * sin2_theta * (x0 * y3 + x3 * y0) + (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * sin2_theta * sin2_theta) * x3 * y3;
		}

	  public:
		const Type a_, a2_, a4_;
		const Type sqrt_delta_a2_;
		const Type u_plus_1, u_plus, u_minus_1, u_minus, u_r;
		const Type r_Q_, r_Q2_, r_Q4_;
		KerrNewman(Type spin, Type charge) : a_(spin), a2_(a_ * a_), a4_(a2_ * a2_), sqrt_delta_a2_(std::sqrt(1. - a2_)), u_plus_1(1. + sqrt_delta_a2_), u_plus(1. / u_plus_1), u_minus_1(a2_ * u_plus), u_minus(1. / u_minus_1), u_r(-0.5 / sqrt_delta_a2_), r_Q_(0.5 * std::sqrt(M_1_PI) * charge), r_Q2_(r_Q_ * r_Q_), r_Q4_(r_Q2_ * r_Q2_) {}
		std::string Name() const override {
			return "Kerr-Newman";
		}
		int MetricTensor(const boost::numeric::ublas::bounded_vector<Type, 8> &position, boost::numeric::ublas::bounded_matrix<Type, 4, 4> &metric) override {
			const Type r = position[1], r2 = r * r, Delta = r2 - 2. * r + a2_ + r_Q2_, rho2 = r2 + a2_ * Power2(std::cos(position[2])), sin2_theta = Power2(std::sin(position[2]));
			if (rho2 == 0. || Delta == 0.)
				return Status::DIVISION_BY_ZERO;
			const Type r_rho_2 = (2. * r - r_Q2_) / rho2;
			metric.clear();
			metric(0, 0) = -(1. - r_rho_2);
			metric(0, 3) = -r_rho_2 * a_ * sin2_theta;
			metric(1, 1) = rho2 / Delta;
			metric(2, 2) = rho2;
			metric(3, 0) = -r_rho_2 * a_ * sin2_theta;
			metric(3, 3) = (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * sin2_theta;
			return Status::SUCCESS;
		}
		Type DistanceSquare(const boost::numeric::ublas::bounded_vector<Type, 8> &x, const boost::numeric::ublas::bounded_vector<Type, 8> &y, bool include_time) override {
			const Type r = x[1], r2 = Power2(r), d0 = x[0] - y[0], d3 = PhiDifference(x[3] - y[3]);
			const Type Delta = r2 - 2. * r + a2_ + r_Q2_, rho2 = r2 + a2_ * Power2(std::cos(x[2])), r_rho_2 = (2. * r - r_Q2_) / rho2, sin2_theta = Power2(std::sin(x[2]));
			if (include_time)
				return (r_rho_2 - 1.) * Power2(d0) - 2. * r_rho_2 * a_ * sin2_theta * d0 * d3 + (rho2 / Delta) * Power2(x[1] - y[1]) + rho2 * Power2(x[2] - y[2]) + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * Power2(sin2_theta)) * Power2(d3);
			return (rho2 / Delta) * Power2(x[1] - y[1]) + rho2 * Power2(x[2] - y[2]) + ((r2 + a2_) * sin2_theta + 2. * r * a2_ * Power2(sin2_theta) / rho2) * Power2(d3);
		}
		int LagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type dt_dtau = 1. / y[4], r2 = Power2(y[1]), sin2_theta = Power2(std::sin(y[2]));
			const Type Delta = r2 - 2. * y[1] + a2_ + r_Q2_, rho2 = r2 + a2_ * Power2(std::cos(y[2]));
			if (rho2 == 0. || Delta == 0.)
				return Status::DIVISION_BY_ZERO;
			const Type r_rho_2 = (2. * y[1] - r_Q2_) / rho2;
			y[4] = (y[4] - 1. + r_rho_2 * (1. - a_ * sin2_theta * y[7])) * dt_dtau; // 1 + p_t
			y[5] *= rho2 / Delta * dt_dtau;
			y[6] *= rho2 * dt_dtau;
			y[7] = (-r_rho_2 * a_ + (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * y[7]) * sin2_theta * dt_dtau;
			return Status::SUCCESS;
		}
		int HamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type pt = y[4] - 1., r2 = Power2(y[1]), sin2_theta = Power2(std::sin(y[2]));
			const Type Delta = r2 - 2. * y[1] + a2_ + r_Q2_, rho2 = r2 + a2_ * Power2(std::cos(y[2]));
			if (rho2 == 0. || Delta == 0.)
				return Status::DIVISION_BY_ZERO;
			const Type rho_2 = 1. / rho2, r_rho_2 = (2. * y[1] - r_Q2_) / rho2;
			y[4] = -Delta / ((Delta + r_rho_2 * (a2_ + r2)) * pt + r_rho_2 * a_ * y[7]);
			y[5] *= Delta * y[4] * rho_2;
			y[6] *= y[4] * rho_2;
			if (y[7] == 0.)
				y[7] = (-r_rho_2 * a_ * pt + (1. - r_rho_2) / sin2_theta * y[7]) / Delta * y[4];
			else if (sin2_theta == 0.)
				return Status::DIVISION_BY_ZERO;
			else
				y[7] = (-r_rho_2 * a_ * pt + (1. - r_rho_2) / sin2_theta * y[7]) / Delta * y[4];
			return Status::SUCCESS;
		}
		int FastTrace(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type r_object, const Type theta_object, const Type phi_object, Type &alpha, Type &beta, boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time) override {
			return Status::FAILURE;
		}
		Type Energy(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				if (time == T)
					return (1. - (2. * y[1] - r_Q2_) / (Power2(y[1]) + a2_ * Power2(std::cos(y[2]))) * (1. - a_ * Power2(std::sin(y[2])) * y[7])) / y[4];
				// time == TAU
				return y[4] - (2. * y[1] - r_Q2_) / (Power2(y[1]) + a2_ * Power2(std::cos(y[2]))) * (y[4] - a_ * Power2(std::sin(y[2])) * y[7]);
			} // dynamics == HAMILTONIAN
			return 1. - y[4];
		}
		Type AngularMomentum(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				const Type r2 = Power2(y[1]), sin2_theta = Power2(std::sin(y[2]));
				const Type r_rho_2 = (2. * y[1] - r_Q2_) / (r2 + a2_ * Power2(std::cos(y[2])));
				if (time == T)
					return (-r_rho_2 * a_ + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta / y[4];
				// time == TAU
				return (-r_rho_2 * a_ * y[4] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta;
			} // dynamics == HAMILTONIAN
			return y[7];
		}
		Type CarterConstant(const boost::numeric::ublas::bounded_vector<Type, 8> &y, const Type mu2, TimeSystem time, DynamicalSystem dynamics) override {
			return std::numeric_limits<Type>::quiet_NaN();
		}
		int NormalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type r = y[1], r2 = Power2(r);
			const Type sin2_theta = Power2(std::sin(y[2]));
			const Type rho2 = r2 + a2_ * Power2(std::cos(y[2]));
			const Type r_rho_2 = (2. * r - r_Q2_) / rho2;
			y[4] = std::sqrt(1. - r_rho_2 + 2. * r_rho_2 * a_ * sin2_theta * y[7] - (rho2 / (r2 - 2. * r + a2_ + r_Q2_) * Power2(y[5]) + rho2 * Power2(y[6]) + ((a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * sin2_theta) * Power2(y[7])));
			return isnan(y[4]) ? Status::DOMAIN_ERROR : Status::SUCCESS;
		}
		int NormalizeNullGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y, Type frequency = 1.) override {
			const Type r = y[1], r2 = Power2(r);
			const Type sin2_theta = Power2(std::sin(y[2]));
			const Type rho2 = r2 + a2_ * Power2(std::cos(y[2]));
			const Type r_rho_2 = (2. * r - r_Q2_) / rho2;
			const Type a = rho2 / (r2 - 2. * r + a2_ + r_Q2_) * Power2(y[5]) + rho2 * Power2(y[6]) + ((a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * sin2_theta) * Power2(y[7]);
			const Type b = -2. * r_rho_2 * a_ * sin2_theta * y[7];
			const Type c = r_rho_2 - 1.;
			const Type coefficient = std::copysign(0.5 / a, frequency) * (-b + std::sqrt(b * b - 4. * a * c));
			if (isnan(coefficient))
				return Status::DOMAIN_ERROR;
			y[4] = frequency;
			y[5] *= coefficient;
			y[6] *= coefficient;
			y[7] *= coefficient;
			return Status::SUCCESS;
		}
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> GetIntegrationSystem(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override {
			throw std::invalid_argument(fmt::format("KerrNewman::GetIntegrationSystem({}, {}, {}) invaild", time, dynamics, motion));
		}
	};

	template <typename Type>
	class KerrTaubNUT : public Metric<Type> {
	  protected:
		Type DotProduct3D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x1, Type x2, Type x3, Type y1, Type y2, Type y3) override {
			const Type r = position[1], sin_theta = std::sin(position[2]), cos_theta = std::cos(position[2]);
			const Type r2 = Power2(r), sin2_theta = Power2(sin_theta);
			const Type Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + Power2(l_ + a_ * cos_theta), chi = a_ * sin2_theta - 2. * l_ * cos_theta;
			const Type rho_2 = 1. / rho2;
			return (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 + (Power2(rho2 + a_ * chi) * sin2_theta - chi * chi * Delta) * rho_2 * x3 * y3;
		}
		Type DotProduct4D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x0, Type x1, Type x2, Type x3, Type y0, Type y1, Type y2, Type y3) override {
			const Type r = position[1], sin_theta = std::sin(position[2]), cos_theta = std::cos(position[2]);
			const Type r2 = Power2(r), sin2_theta = Power2(sin_theta);
			const Type Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + Power2(l_ + a_ * cos_theta), chi = a_ * sin2_theta - 2. * l_ * cos_theta;
			const Type rho_2 = 1. / rho2;
			return (a2_ * sin2_theta - Delta) * rho_2 * x0 * y0 - 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * rho_2 * (x0 * y3 + x3 * y0) + (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 + (Power2(rho2 + a_ * chi) * sin2_theta - chi * chi * Delta) * rho_2 * x3 * y3;
		}

	  public:
		const Type a_, a2_, a4_;
		const Type sqrt_delta_a2_;
		const Type u_plus_1, u_plus, u_minus_1, u_minus, u_r;
		const Type l_, l2_, l4_;
		KerrTaubNUT(Type spin, Type NUT) : a_(spin), a2_(a_ * a_), a4_(a2_ * a2_), sqrt_delta_a2_(std::sqrt(1. - a2_)), u_plus_1(1. + sqrt_delta_a2_), u_plus(1. / u_plus_1), u_minus_1(a2_ * u_plus), u_minus(1. / u_minus_1), u_r(-0.5 / sqrt_delta_a2_), l_(NUT), l2_(l_ * l_), l4_(l2_ * l2_) {}
		std::string Name() const override {
			return "Kerr-Taub-NUT";
		}
		int MetricTensor(const boost::numeric::ublas::bounded_vector<Type, 8> &position, boost::numeric::ublas::bounded_matrix<Type, 4, 4> &metric) override {
			return Status::NOT_IMPLEMENTED;
		}
		Type DistanceSquare(const boost::numeric::ublas::bounded_vector<Type, 8> &x, const boost::numeric::ublas::bounded_vector<Type, 8> &y, bool include_time) override {
			const Type r = x[1], sin_theta = std::sin(x[2]), cos_theta = std::cos(x[2]), d0 = x[0] - y[0], d3 = PhiDifference(x[3] - y[3]);
			const Type r2 = Power2(r), sin2_theta = Power2(sin_theta);
			const Type Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + Power2(l_ + a_ * cos_theta), chi = a_ * sin2_theta - 2. * l_ * cos_theta;
			const Type rho_2 = 1. / rho2;
			if (include_time)
				return (a2_ * sin2_theta - Delta) * rho_2 * Power2(d0) - 4. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * rho_2 * d0 * d3 + (rho2 / Delta) * Power2(x[1] - y[1]) + rho2 * Power2(x[2] - y[2]) + (Power2(rho2 + a_ * chi) * sin2_theta - Power2(chi) * Delta) * rho_2 * Power2(d3);
			return (rho2 / Delta) * Power2(x[1] - y[1]) + rho2 * Power2(x[2] - y[2]) + (Power2(rho2 + a_ * chi) * sin2_theta - Power2(chi) * Delta) * rho_2 * Power2(d3);
		}
		int LagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type dt_dtau = 1. / y[4], r = y[1], r2 = Power2(r);
			const Type sin2_theta = Power2(std::sin(y[2])), cos_theta = CosTheta(y[2]);
			const Type Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + Power2(l_ + a_ * cos_theta);
			if (rho2 == 0. || Delta == 0.)
				return Status::DIVISION_BY_ZERO;
			const Type rho_2 = 1. / rho2;
			y[4] = (y[4] * rho2 - Delta + a2_ * sin2_theta - 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]) * rho_2 * dt_dtau; // 1 + p_t
			y[5] *= rho2 / Delta * dt_dtau;
			y[6] *= rho2 * dt_dtau;
			y[7] = (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) + (Power2(r2 + l2_ + a2_) * sin2_theta - Power2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * y[7]) * rho_2 * dt_dtau;
			return Status::SUCCESS;
		}
		int HamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type pt = y[4] - 1., r = y[1], r2 = Power2(r);
			const Type sin2_theta = Power2(std::sin(y[2])), cos_theta = CosTheta(y[2]);
			const Type Delta = r2 - 2. * r - l2_ + a2_, rho2 = r2 + Power2(l_ + a_ * cos_theta);
			if (rho2 == 0. || Delta == 0.)
				return Status::DIVISION_BY_ZERO;
			const Type rho_2 = 1. / rho2;
			y[4] = Delta * rho2 * sin2_theta / ((Delta * Power2(a_ * sin2_theta - 2. * l_ * cos_theta) - Power2(r2 + l2_ + a2_) * sin2_theta) * pt - 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]);
			y[5] *= Delta * rho_2 * y[4];
			y[6] *= rho_2 * y[4];
			y[7] = (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * pt + (Delta - a2_ * sin2_theta) * y[7]) / (Delta * rho2 * sin2_theta) * y[4];
			return Status::SUCCESS;
		}
		int FastTrace(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type r_object, const Type theta_object, const Type phi_object, Type &alpha, Type &beta, boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time) override {
			return Status::FAILURE;
		}
		Type Energy(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				const Type r = y[1], r2 = Power2(r);
				const Type sin2_theta = Power2(std::sin(y[2])), cos_theta = CosTheta(y[2]);
				const Type Delta = r2 - 2. * r - l2_ + a2_;
				if (time == T)
					return (Delta - a2_ * sin2_theta + 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]) / ((r2 + Power2(l_ + a_ * cos_theta)) * y[4]);
				// time == TAU
				return ((Delta - a2_ * sin2_theta) * y[4] + 2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7]) / (r2 + Power2(l_ + a_ * cos_theta));
			} // dynamics == HAMILTONIAN
			return 1. - y[4];
		}
		Type AngularMomentum(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override {
			if (dynamics == LAGRANGIAN) {
				const Type r = y[1], r2 = Power2(r);
				const Type sin2_theta = Power2(std::sin(y[2])), cos_theta = CosTheta(y[2]);
				const Type Delta = r2 - 2. * r - l2_ + a2_;
				if (time == T)
					return (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) + (Power2(r2 + l2_ + a2_) * sin2_theta - Power2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * y[7]) / ((r2 + Power2(l_ + a_ * cos_theta)) * y[4]);
				// time == TAU
				return (-2. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[4] + (Power2(r2 + l2_ + a2_) * sin2_theta - Power2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * y[7]) / (r2 + Power2(l_ + a_ * cos_theta));
			} // dynamics == HAMILTONIAN
			return y[7];
		}
		Type CarterConstant(const boost::numeric::ublas::bounded_vector<Type, 8> &y, const Type mu2, TimeSystem time, DynamicalSystem dynamics) override {
			return std::numeric_limits<Type>::quiet_NaN();
		}
		int NormalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y) override {
			const Type r = y[1], r2 = Power2(r);
			const Type sin2_theta = Power2(std::sin(y[2])), cos_theta = CosTheta(y[2]);
			const Type Delta = r2 - 2. * r - l2_ + a2_;
			const Type rho2 = r2 + Power2(l_ + a_ * cos_theta);
			// y[7] += 2. * a_ * r / (Power2(a2_ + r2) - a2_ * Delta * sin2_theta);
			y[4] = std::sqrt(((Delta - a2_ * sin2_theta) + 4. * ((r + l2_) * a_ * sin2_theta + Delta * l_ * cos_theta) * y[7] - (Power2(r2 + l2_ + a2_) * sin2_theta - Power2(a_ * sin2_theta - 2. * l_ * cos_theta) * Delta) * Power2(y[7])) / rho2 - rho2 * (Power2(y[5]) / Delta + Power2(y[6])));
			return isnan(y[4]) ? Status::DOMAIN_ERROR : Status::SUCCESS;
		}
		int NormalizeNullGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y, Type frequency = 1.) override {
			const Type r = y[1], r2 = Power2(r);
			const Type sin2_theta = Power2(std::sin(y[2])), cos_theta = CosTheta(y[2]);
			const Type Delta = r2 - 2. * r - l2_ + a2_;
			const Type l_a_cos_theta = l_ + a_ * cos_theta;
			const Type rho2 = r2 + Power2(l_a_cos_theta), rho_2 = 1. / rho2;
			const Type chi = a_ * sin2_theta - 2. * l_ * cos_theta, rho2_a_chi = r2 + l2_ + a2_;
			const Type a = rho2 / Delta * Power2(y[5]) + rho2 * Power2(y[6]) + rho_2 * (Power2(rho2_a_chi) * sin2_theta - Power2(chi) * Delta) * Power2(y[7]);
			const Type b = -4. * rho_2 * ((r + l2_) * chi + l_ * cos_theta * rho2_a_chi) * y[7];
			const Type c = -rho_2 * (Delta - a2_ * sin2_theta);
			const Type coefficient = std::copysign(0.5 / a, frequency) * (-b + std::sqrt(b * b - 4. * a * c));
			if (isnan(coefficient))
				return Status::DOMAIN_ERROR;
			y[4] = frequency;
			y[5] *= coefficient;
			y[6] *= coefficient;
			y[7] *= coefficient;
			return Status::SUCCESS;
		}
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> GetIntegrationSystem(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override {
			if (motion != GEODESIC)
				throw std::invalid_argument(fmt::format("KerrNewman::GetIntegrationSystem({}, {}, {}) invaild", time, dynamics, motion));
			if (time == T) {
				if (dynamics == LAGRANGIAN)
					if (motion == GEODESIC)
						return [a = this->a_, a2 = this->a2_, l = this->l_, l2 = this->l2_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = y[6]; // d\theta/dt
							dydt[3] = y[7]; // d\phi/dt
							const Type r = y[1], r2 = Power2(r);
							const Type sin_theta = SinTheta(y[2]), sin_1_theta = 1. / sin_theta, sin2_theta = Power2(sin_theta), cos_theta = CosTheta(y[2]), cos2_theta = Power2(cos_theta);
							const Type Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
							const Type l_a_cos_theta = l + a * cos_theta, l_a_cos_theta2 = Power2(l_a_cos_theta);
							const Type rho2 = r2 + l_a_cos_theta2, rho_2 = 1. / rho2, rho4 = Power2(rho2), rho_4 = Power2(rho_2), rho_6 = rho_2 * rho_4;
							const Type chi = a * sin2_theta - 2. * l * cos_theta, rho2_a_chi = r2 + l2 + a2;
							const Type rho2_r_Delta = r * (r + 2. * l * l_a_cos_theta) - l_a_cos_theta2, rho2_a_cos_theta = (2. * r + l * l_a_cos_theta) * l_a_cos_theta - r2 * l;
							const Type dydt4 = 2. * rho_4 * (Delta_1 * rho2_a_chi * rho2_r_Delta * y[5] * (1. - chi * y[7]) - sin_1_theta * chi * rho2_a_cos_theta * y[6] * (1. - chi * y[7]) - 2. * Delta_1 * rho2 * r * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * y[5] * y[7] - sin_1_theta * rho4 * l * (1. + cos2_theta) * y[6] * y[7]);
							// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
							dydt[4] = dydt4 * y[4];
							// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[5] = -Delta * rho2_r_Delta * rho_6 * Power2(1. - chi * y[7]) + (rho2_r_Delta - r * a2 * sin2_theta) * Delta_1 * rho_2 * Power2(y[5]) + 2. * a * sin_theta * l_a_cos_theta * rho_2 * y[5] * y[6] + r * Delta * rho_2 * Power2(y[6]) + Delta * r * sin2_theta * rho_2 * Power2(y[7]) + dydt4 * y[5];
							// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[6] = rho2_a_cos_theta * sin_theta * rho_6 * (a - 2. * rho2_a_chi * y[7]) - a * sin_theta * l_a_cos_theta * rho_2 * (Delta_1 * Power2(y[5]) - Power2(y[6])) - 2. * r * rho_2 * y[5] * y[6] - (rho2 * (Delta * chi * l_a_cos_theta - Power2(rho2_a_chi) * cos_theta) - 2. * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * rho2_a_chi * l_a_cos_theta) * sin_theta * rho_6 * Power2(y[7]) + dydt4 * y[6];
							// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[7] = -2. * a * rho2_r_Delta * Delta_1 * rho_4 * y[5] * (1. - chi * y[7]) + 2. * rho2_a_cos_theta * rho_4 * sin_1_theta * y[6] * (1. - chi * y[7]) - 2. * (1. - a2 * sin2_theta * Delta_1) * r * rho_2 * y[5] * y[7] - 2. * cos_theta * sin_1_theta * y[6] * y[7] + dydt4 * y[7];
						};
					else if (motion == HELICAL)
						return [a = this->a_, a2 = this->a2_, l = this->l_, l2 = this->l2_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) { // [TODO]: CHECK
							dydt[0] = y[4];																																												   // d\tau/dt
							dydt[1] = y[5];																																												   // dr/dt
							dydt[2] = 0.;																																												   // d\theta/dt = 0.
							dydt[3] = y[7];																																												   // d\phi/dt
							const Type r = y[1], r2 = Power2(r);
							const Type sin2_theta = Power2(std::sin(y[2])), cos_theta = CosTheta(y[2]);
							const Type Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
							const Type l_a_cos_theta2 = Power2(l + a * cos_theta);
							const Type rho2 = r2 + l_a_cos_theta2, rho_2 = 1. / rho2;
							const Type chi = a * sin2_theta - 2. * l * cos_theta, rho2_a_chi = r2 + l2 + a2;
							const Type g00 = -rho_2 * (Delta - a2 * sin2_theta), g11 = rho2 * Delta_1, g03 = -g00 * chi - a * sin2_theta, g33 = -g03 * chi + rho2_a_chi * sin2_theta;
							const Type L_1 = y[4] / (g03 + g33 * y[7]), L_2 = Power2(L_1);
							const Type A = g33 * (g33 * L_2 + 1.), A_1 = 1. / A, B = 2. * g03 * (g33 * L_2 + 1.), C = g00 + g11 * Power2(y[5]) + Power2(g03) * L_2;
							const Type dg00_dr = -2. * rho_2 * ((r - 1.) + g00 * r), dg11_dr = 2. * Delta_1 * (r - g11 * (r - 1.)), dg03_dr = -dg00_dr * chi, dg33_dr = 2. * r * sin2_theta - dg03_dr * chi;
							const Type dA_dr = dg33_dr * (2. * g33 * L_2 + 1.), dB_dr = 2. * (dg03_dr * (g33 * L_2 + 1.) + g03 * dg33_dr * L_2), dC_dr = dg00_dr + dg11_dr * Power2(y[5]) + 2. * g03 * dg03_dr * L_2;
							// d^2r/dt^2 = 0.
							dydt[5] = 0.;
							// d^2\theta/dt^2 = 0.
							dydt[6] = 0.;
							// d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
							dydt[7] = y[5] * A_1 * (-y[7] * dA_dr + 0.5 * (-dB_dr + (B * dB_dr - 2. * (A * dC_dr + C * dA_dr)) / (2. * A * y[7] + B)));
							// d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt = d((g03 + g33 * d\phi/dt) / L)/dr * dr/dt
							dydt[4] = (y[5] * (dg03_dr + dg33_dr * y[7]) + g33 * dydt[7]) * L_1;
						};
			} else if (time == TAU) {
				if (dynamics == LAGRANGIAN && motion == GEODESIC)
					return [a = this->a_, a2 = this->a2_, l = this->l_, l2 = this->l2_](const boost::numeric::ublas::bounded_vector<Type, 8> &y, boost::numeric::ublas::bounded_vector<Type, 8> &dydt, const Type t) {
						dydt[0] = y[4]; // dt/d\tau
						dydt[1] = y[5]; // dr/d\tau
						dydt[2] = y[6]; // d\theta/d\tau
						dydt[3] = y[7]; // d\phi/d\tau
						const Type r = y[1], r2 = Power2(r);
						const Type sin_theta = SinTheta(y[2]), sin_1_theta = 1. / sin_theta, sin2_theta = Power2(sin_theta), cos_theta = CosTheta(y[2]), cos2_theta = Power2(cos_theta);
						const Type Delta = r2 - 2. * r - l2 + a2, Delta_1 = 1. / Delta;
						const Type l_a_cos_theta = l + a * cos_theta, l_a_cos_theta2 = Power2(l_a_cos_theta);
						const Type rho2 = r2 + l_a_cos_theta2, rho_2 = 1. / rho2, rho4 = Power2(rho2), rho_4 = Power2(rho_2), rho_6 = rho_2 * rho_4;
						const Type chi = a * sin2_theta - 2. * l * cos_theta, rho2_a_chi = r2 + l2 + a2;
						const Type rho2_r_Delta = r * (r + 2. * l * l_a_cos_theta) - l_a_cos_theta2, rho2_a_cos_theta = (2. * r + l * l_a_cos_theta) * l_a_cos_theta - r2 * l;
						// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
						dydt[4] = -2. * rho_4 * (Delta_1 * rho2_a_chi * rho2_r_Delta * y[5] * (y[4] - chi * y[7]) - sin_1_theta * chi * rho2_a_cos_theta * y[6] * (y[4] - chi * y[7]) - 2. * Delta_1 * rho2 * r * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * y[5] * y[7] - sin_1_theta * rho4 * l * (1. + cos2_theta) * y[6] * y[7]);
						// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
						dydt[5] = -Delta * rho2_r_Delta * rho_6 * Power2(y[4] - chi * y[7]) + (rho2_r_Delta - r * a2 * sin2_theta) * Delta_1 * rho_2 * Power2(y[5]) + 2. * a * sin_theta * l_a_cos_theta * rho_2 * y[5] * y[6] + r * Delta * rho_2 * Power2(y[6]) + Delta * r * sin2_theta * rho_2 * Power2(y[7]);
						// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
						dydt[6] = rho2_a_cos_theta * sin_theta * rho_6 * y[4] * (a * y[4] - 2. * rho2_a_chi * y[7]) - a * sin_theta * l_a_cos_theta * rho_2 * (Delta_1 * Power2(y[5]) - Power2(y[6])) - 2. * r * rho_2 * y[5] * y[6] - (rho2 * (Delta * chi * l_a_cos_theta - Power2(rho2_a_chi) * cos_theta) - 2. * ((r + l2) * a * sin2_theta + Delta * l * cos_theta) * rho2_a_chi * l_a_cos_theta) * sin_theta * rho_6 * Power2(y[7]);
						// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
						dydt[7] = -2. * a * rho2_r_Delta * Delta_1 * rho_4 * y[5] * (y[4] - chi * y[7]) + 2. * rho2_a_cos_theta * rho_4 * sin_1_theta * y[6] * (y[4] - chi * y[7]) - 2. * (1. - a2 * sin2_theta * Delta_1) * r * rho_2 * y[5] * y[7] - 2. * cos_theta * sin_1_theta * y[6] * y[7];
					};
			}
			throw std::invalid_argument(fmt::format("KerrNewman::GetIntegrationSystem({}, {}, {}) invaild", time, dynamics, motion));
		}
	};

	template <typename Type>
	class Hayward : public Metric<Type> {
	  protected:
		Type DotProduct3D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x1, Type x2, Type x3, Type y1, Type y2, Type y3) override;
		Type DotProduct4D(const boost::numeric::ublas::bounded_vector<Type, 8> &position, Type x0, Type x1, Type x2, Type x3, Type y0, Type y1, Type y2, Type y3) override;

	  public:
		const Type a_, a2_, a4_;
		const Type sqrt_delta_a2_;
		const Type u_plus_1, u_plus, u_minus_1, u_minus, u_r;
		const Type alpha_, beta_, g_;
		Hayward(Type alpha, Type beta, Type charge);
		std::string Name() override;
		int MetricTensor(const boost::numeric::ublas::bounded_vector<Type, 8> &position, boost::numeric::ublas::bounded_matrix<Type, 4, 4> &metric) override;
		Type DistanceSquare(const boost::numeric::ublas::bounded_vector<Type, 8> &x, const boost::numeric::ublas::bounded_vector<Type, 8> &y, bool include_time) override;
		int LagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override;
		int HamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<Type, 8> &y) override;
		int FastTrace(const Type r_observer, const Type theta_observer, const Type sin_theta_observer, const Type cos_theta_observer, const Type r_object, const Type theta_object, const Type phi_object, Type &alpha, Type &beta, boost::numeric::ublas::bounded_vector<Type, 8> &photon, Type &photon_time) override;
		Type Energy(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override;
		Type AngularMomentum(const boost::numeric::ublas::bounded_vector<Type, 8> &y, TimeSystem time, DynamicalSystem dynamics) override;
		Type CarterConstant(const boost::numeric::ublas::bounded_vector<Type, 8> &y, const Type mu2, TimeSystem time, DynamicalSystem dynamics) override;
		int NormalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y) override;
		int NormalizeNullGeodesic(boost::numeric::ublas::bounded_vector<Type, 8> &y, Type frequency = 1.) override;
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> GetIntegrationSystem(TimeSystem time, DynamicalSystem dynamics, MotionMode motion = GEODESIC) override;
	};
} // namespace SBody

#endif
