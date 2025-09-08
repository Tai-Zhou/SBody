/**
 * @file Object.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_OBJECT_HPP
#define SBODY_OBJECT_HPP

#include <array>
#include <functional>
#include <memory>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include "Metric.hpp"
#include "Utility.hpp"

namespace SBody {
	/**
	 * @brief
	 *
	 */
	template <typename Type>
	class Object {
	  protected:
		std::shared_ptr<Metric<Type>> metric_;

	  public:
		/**
		 * @brief Construct a new Object object
		 *
		 * @param metric
		 */
		Object(std::shared_ptr<Metric<Type>> metric) : metric_(metric) {}
		/**
		 * @brief
		 *
		 * @param current
		 * @param last
		 * @return int
		 */
		virtual int Hit(const boost::numeric::ublas::bounded_vector<Type, 8> &current, const boost::numeric::ublas::bounded_vector<Type, 8> &last) = 0;

		/**
		 * @brief
		 * \f[z=\frac{E_\mathrm{obj}}{E_\mathrm{obs}}-1\f]
		 * where
		 * \f[E_\mathrm{obj}=-g_{\mu\nu}u_\mathrm{ph}^\mu u_\mathrm{obj}^\nu\f]
		 * For the distant observer, the equation is the same. With the tetrad velocity \f$u_\mathrm{obs}=(1,0,0,0)\f$ and metric
		 * \f[
		 * 	g_\mathrm{obs}=
		 * 	\begin{pmatrix}
		 * 		-1 & 0 & 0 & 0 \\
		 * 		0  & 1 & 0 & 0 \\
		 * 		0  & 0 & 1 & 0 \\
		 * 		0  & 0 & 0 & 1
		 * 	\end{pmatrix},
		 * \f]
		 * the energy of the photon can be simplified to
		 * \f[
		 * 	E_\mathrm{obs}=u_\mathrm{ph}^0.
		 * \f]
		 *
		 * @param photon 8 dimensional information of photon
		 * @param photon_time
		 * @return Type
		 */
		virtual Type Redshift(const boost::numeric::ublas::bounded_vector<Type, 8> &photon, TimeSystem photon_time) = 0;
	};

	/**
	 * @brief
	 *
	 */
	template <typename Type>
	class Particle : public Object<Type> {
	  protected:
		TimeSystem time_;
		DynamicalSystem coordinate_;
		/// 8 dimensional information
		boost::numeric::ublas::bounded_vector<Type, 8> position_;
		/// if the position fixed
		const bool fixed_;
		/// integration
		std::function<void(const boost::numeric::ublas::bounded_vector<Type, 8> &, boost::numeric::ublas::bounded_vector<Type, 8> &, const Type)> integration_system_;
		boost::numeric::odeint::controlled_runge_kutta<boost::numeric::odeint::runge_kutta_dopri5<boost::numeric::ublas::bounded_vector<Type, 8>>> integration_stepper_;

	  public:
		/**
		 * @brief Construct a new star object
		 *
		 * @param metric
		 * @param radius radius
		 * @param fixed whether the position of the star is fixed
		 */
		Particle(std::shared_ptr<Metric<Type>> metric, TimeSystem time, DynamicalSystem coordinate, bool fixed = false) : Object<Type>(metric), time_(time), coordinate_(coordinate), fixed_(fixed) {}
		int Position(Type *position) {
			std::copy(position_.begin(), position_.end(), position);
			return Status::SUCCESS;
		}
		int Position(boost::numeric::ublas::bounded_vector<Type, 8> &position) {
			position = position_;
			return Status::SUCCESS;
		}
		int InitializeKeplerian(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type true_anomaly, Type observer_inclination = 0., Type observer_rotation = 0.) {
			integration_system_ = this->metric_->GetIntegrationSystem(time_, coordinate_, GEODESIC);
			position_[0] = 0.;
			periapsis += true_anomaly;
			Type r = a * (1. - e * e) / (1. + e * std::cos(true_anomaly));
			const Type tp1 = -r * std::cos(periapsis), tp2 = -r * std::sin(periapsis) * std::cos(inclination);
			const Type xp1 = tp1 * std::cos(ascending_node) - tp2 * std::sin(ascending_node), xp2 = tp2 * std::cos(ascending_node) + tp1 * std::sin(ascending_node), xp3 = -r * std::sin(periapsis) * std::sin(inclination);
			position_[1] = (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::cos(observer_inclination) + xp3 * std::sin(observer_inclination);
			position_[2] = xp2 * std::cos(observer_rotation) - xp1 * std::sin(observer_rotation);
			position_[3] = xp3 * std::cos(observer_inclination) - (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::sin(observer_inclination);
			if (fixed_) {
				position_[5] = 0.;
				position_[6] = 0.;
				position_[7] = 0.;
			} else {
				const Type vphi = std::sqrt((1 - e * e) * a) / r, vr = std::copysign(SquareRoot(2. / r - 1. / a - vphi * vphi), PhiDifference(true_anomaly));
				const Type tp5 = vphi * std::sin(periapsis) - vr * std::cos(periapsis), tp6 = -(vphi * std::cos(periapsis) + vr * std::sin(periapsis)) * std::cos(inclination);
				const Type xp5 = tp5 * std::cos(ascending_node) - tp6 * std::sin(ascending_node), xp6 = tp5 * std::sin(ascending_node) + tp6 * std::cos(ascending_node), xp7 = -(vphi * std::cos(periapsis) + vr * std::sin(periapsis)) * std::sin(inclination);
				position_[5] = (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::cos(observer_inclination) + xp7 * std::sin(observer_inclination);
				position_[6] = xp6 * std::cos(observer_rotation) - xp5 * std::sin(observer_rotation);
				position_[7] = xp7 * std::cos(observer_inclination) - (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::sin(observer_inclination);
			}
			CartesianToSpherical(position_);
			return this->metric_->NormalizeTimelikeGeodesic(position_);
		}
		int InitializeGeodesic(Type orbital_radius, Type inclination, Type periapsis, Type ascending_node, Type v_r, Type v_phi, Type observer_inclination = 0., Type observer_rotation = 0.) {
			integration_system_ = this->metric_->GetIntegrationSystem(time_, coordinate_, GEODESIC);
			position_[0] = 0.;
			const Type tp1 = -orbital_radius * std::cos(periapsis), tp2 = -orbital_radius * std::sin(periapsis) * std::cos(inclination);
			const Type xp1 = tp1 * std::cos(ascending_node) - tp2 * std::sin(ascending_node), xp2 = tp2 * std::cos(ascending_node) + tp1 * std::sin(ascending_node), xp3 = -orbital_radius * std::sin(periapsis) * std::sin(inclination);
			position_[1] = (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::cos(observer_inclination) + xp3 * std::sin(observer_inclination);
			position_[2] = xp2 * std::cos(observer_rotation) - xp1 * std::sin(observer_rotation);
			position_[3] = xp3 * std::cos(observer_inclination) - (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::sin(observer_inclination);
			if (fixed_) {
				position_[5] = 0.;
				position_[6] = 0.;
				position_[7] = 0.;
			} else {
				const Type tp5 = v_phi * std::sin(periapsis) - v_r * std::cos(periapsis), tp6 = -(v_phi * std::cos(periapsis) + v_r * std::sin(periapsis)) * std::cos(inclination);
				const Type xp5 = tp5 * std::cos(ascending_node) - tp6 * std::sin(ascending_node), xp6 = tp5 * std::sin(ascending_node) + tp6 * std::cos(ascending_node), xp7 = -(v_phi * std::cos(periapsis) + v_r * std::sin(periapsis)) * std::sin(inclination);
				position_[5] = (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::cos(observer_inclination) + xp7 * std::sin(observer_inclination);
				position_[6] = xp6 * std::cos(observer_rotation) - xp5 * std::sin(observer_rotation);
				position_[7] = xp7 * std::cos(observer_inclination) - (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::sin(observer_inclination);
			}
			CartesianToSpherical(position_);
			return this->metric_->NormalizeTimelikeGeodesic(position_);
		}
		int InitializeSchwarzschildKeplerianPericenter(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type observer_inclination, Type observer_rotation) {
			integration_system_ = this->metric_->GetIntegrationSystem(time_, coordinate_, GEODESIC);
			position_[0] = 0.;
			Type r = a * (1. - e); // pericenter
			const Type tp1 = -r * std::cos(periapsis), tp2 = -r * std::sin(periapsis) * std::cos(inclination);
			const Type xp1 = tp1 * std::cos(ascending_node) - tp2 * std::sin(ascending_node), xp2 = tp2 * std::cos(ascending_node) + tp1 * std::sin(ascending_node), xp3 = -r * std::sin(periapsis) * std::sin(inclination);
			position_[1] = (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::cos(observer_inclination) + xp3 * std::sin(observer_inclination);
			position_[2] = xp2 * std::cos(observer_rotation) - xp1 * std::sin(observer_rotation);
			position_[3] = xp3 * std::cos(observer_inclination) - (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::sin(observer_inclination);
			if (fixed_) {
				position_[5] = 0.;
				position_[6] = 0.;
				position_[7] = 0.;
			} else {
				const Type E = (a - 2.) / std::sqrt(a * (a - 3.));
				const Type vphi = std::sqrt(E * E * r * (r - 2.) - Power2(r - 2.)) / (r * E);
				const Type tp5 = vphi * std::sin(periapsis), tp6 = -vphi * std::cos(periapsis) * std::cos(inclination);
				const Type xp5 = tp5 * std::cos(ascending_node) - tp6 * std::sin(ascending_node), xp6 = tp5 * std::sin(ascending_node) + tp6 * std::cos(ascending_node), xp7 = -vphi * std::cos(periapsis) * std::sin(inclination);
				position_[5] = (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::cos(observer_inclination) + xp7 * std::sin(observer_inclination);
				position_[6] = xp6 * std::cos(observer_rotation) - xp5 * std::sin(observer_rotation);
				position_[7] = xp7 * std::cos(observer_inclination) - (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::sin(observer_inclination);
			}
			CartesianToSpherical(position_);
			return this->metric_->NormalizeTimelikeGeodesic(position_);
		}
		int InitializeSchwarzschildKeplerianApocenter(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type observer_inclination, Type observer_rotation) {
			integration_system_ = this->metric_->GetIntegrationSystem(time_, coordinate_, GEODESIC);
			position_[0] = 0.;
			Type r = a * (1. + e); // apocenter
			const Type tp1 = r * std::cos(periapsis), tp2 = r * std::sin(periapsis) * std::cos(inclination);
			const Type xp1 = tp1 * std::cos(ascending_node) - tp2 * std::sin(ascending_node), xp2 = tp2 * std::cos(ascending_node) + tp1 * std::sin(ascending_node), xp3 = r * std::sin(periapsis) * std::sin(inclination);
			position_[1] = (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::cos(observer_inclination) + xp3 * std::sin(observer_inclination);
			position_[2] = xp2 * std::cos(observer_rotation) - xp1 * std::sin(observer_rotation);
			position_[3] = xp3 * std::cos(observer_inclination) - (xp1 * std::cos(observer_rotation) + xp2 * std::sin(observer_rotation)) * std::sin(observer_inclination);
			if (fixed_) {
				position_[5] = 0.;
				position_[6] = 0.;
				position_[7] = 0.;
			} else {
				const Type E = (a - 2.) / std::sqrt(a * (a - 3.));
				const Type vphi = std::sqrt(E * E * r * (r - 2.) - Power2(r - 2.)) / (r * E);
				const Type tp5 = -vphi * std::sin(periapsis), tp6 = vphi * std::cos(periapsis) * std::cos(inclination);
				const Type xp5 = tp5 * std::cos(ascending_node) - tp6 * std::sin(ascending_node), xp6 = tp5 * std::sin(ascending_node) + tp6 * std::cos(ascending_node), xp7 = vphi * std::cos(periapsis) * std::sin(inclination);
				position_[5] = (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::cos(observer_inclination) + xp7 * std::sin(observer_inclination);
				position_[6] = xp6 * std::cos(observer_rotation) - xp5 * std::sin(observer_rotation);
				position_[7] = xp7 * std::cos(observer_inclination) - (xp5 * std::cos(observer_rotation) + xp6 * std::sin(observer_rotation)) * std::sin(observer_inclination);
			}
			CartesianToSpherical(position_);
			return this->metric_->NormalizeTimelikeGeodesic(position_);
		}
		int InitializeKeplerianHarmonic(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type true_anomaly, Type observer_inclination, Type observer_rotation) {
			InitializeKeplerian(a, e, inclination, periapsis, ascending_node, true_anomaly, observer_inclination, observer_rotation);
			position_[1] += 1.;
			return this->metric_->NormalizeTimelikeGeodesic(position_);
		}
		int InitializeCircular(Type r, Type phi, Type v_phi_ratio) {
			integration_system_ = this->metric_->GetIntegrationSystem(time_, coordinate_, CIRCULAR);
			position_[0] = 0.;
			position_[1] = r;
			position_[2] = boost::math::constants::half_pi<Type>();
			position_[3] = phi;
			if (fixed_) {
				position_[5] = 0.;
				position_[6] = 0.;
				position_[7] = 0.;
			} else {
				position_[5] = 0.;
				position_[6] = 0.;
				position_[7] = v_phi_ratio / (r * std::sqrt(r));
			}
			return this->metric_->NormalizeTimelikeGeodesic(position_);
		}
		int InitializeHelical(Type r, Type theta, Type phi, Type v_r, Type v_phi) {
			integration_system_ = this->metric_->GetIntegrationSystem(time_, coordinate_, HELICAL);
			position_[0] = 0.;
			position_[1] = r;
			position_[2] = theta;
			position_[3] = phi;
			if (fixed_) {
				position_[5] = 0.;
				position_[6] = 0.;
				position_[7] = 0.;
			} else {
				position_[5] = v_r;
				position_[6] = 0.;
				position_[7] = v_phi;
			}
			return this->metric_->NormalizeTimelikeGeodesic(position_);
		}
		int InitializeCylindrical(Type r, Type phi, Type z, Type v_z, Type v_phi) {
			integration_system_ = this->metric_->GetIntegrationSystem(time_, coordinate_, CYLINDRICAL);
			position_[0] = 0.;
			position_[1] = std::hypot(r, z);
			position_[2] = std::atan2(r, z);
			position_[3] = phi;
			if (fixed_) {
				position_[5] = 0.;
				position_[6] = 0.;
				position_[7] = 0.;
			} else {
				position_[5] = v_z * std::cos(position_[2]);
				position_[6] = -v_z * std::sin(position_[2]) / position_[1];
				position_[7] = v_phi;
			}
			return this->metric_->NormalizeTimelikeGeodesic(position_);
		}
		int IntegratorApply(Type &t, Type t1, Type &h) {
			if (fixed_) {
				t = t1;
				return Status::SUCCESS;
			}
			if (h == 0.)
				h = std::copysign(std::max(std::sqrt(std::abs(t1 - t)), 1.), t1 - t);
			Type last_theta = position_[2];
			while ((h > 0 ? t1 - t : t - t1) > 0.) {
				if (std::abs(h) > std::abs(t1 - t))
					h = t1 - t;
				try {
					integration_stepper_.try_step(integration_system_, position_, t, h);
				} catch (const std::exception &e) {
					return Status::FAILURE;
				}
				MapTheta(position_, last_theta);
			}
			return Status::SUCCESS;
		}
		int IntegratorApplyStep(Type &t, const Type h) {
			if (fixed_) {
				t += h;
				return Status::SUCCESS;
			}
			Type start_time = t;
			try {
				integration_stepper_.try_step(integration_system_, position_, t, h);
			} catch (const std::exception &e) {
				return Status::FAILURE;
			}
			if (t == start_time)
				return Status::FAILURE;
			return Status::SUCCESS;
		}
		int Normalize() {
			return this->metric_->NormalizeTimelikeGeodesic(position_);
		}

		/**
		 * @brief Check if the star hit by the photon
		 *
		 * @param current
		 * @param last
		 * @return int
		 */
		int Hit(const boost::numeric::ublas::bounded_vector<Type, 8> &current, const boost::numeric::ublas::bounded_vector<Type, 8> &last) override {
			return false;
		}

		/**
		 * @brief return frequency
		 *
		 * @param photon
		 * @param photon_time
		 * @return Type
		 */
		Type Redshift(const boost::numeric::ublas::bounded_vector<Type, 8> &photon, TimeSystem photon_time) override {
			return this->metric_->Redshift(position_, photon, time_, photon_time);
		}

		int MetricTensor(boost::numeric::ublas::bounded_matrix<Type, 4, 4> &metric) {
			return this->metric_->MetricTensor(position_, metric);
		}
		template <typename VecExp1, typename VecExp2>
		Type DotProduct(const VecExp1 &x, const VecExp2 &y, bool include_time) {
			return this->metric_->DotProduct(position_, x, y, include_time);
		}
		int LocalInertialFrame(boost::numeric::ublas::bounded_matrix<Type, 4, 4> &coordinate) {
			return this->metric_->LocalInertialFrame(position_, time_, coordinate);
		}
		Type Energy() {
			return this->metric_->Energy(position_, time_, coordinate_);
		}
		Type AngularMomentum() {
			return this->metric_->AngularMomentum(position_, time_, coordinate_);
		}
		Type CarterConstant() {
			return this->metric_->CarterConstant(position_, 1., time_, coordinate_);
		}
	};

	template <typename Type>
	class Star : public Particle<Type> {
	  protected:
		/// radius
		const Type radius_;
		/// radius square
		const Type radius_square_;

	  public:
		Star(std::shared_ptr<Metric<Type>> metric, TimeSystem time, DynamicalSystem coordinate, Type radius = 0, bool fixed = false) : Particle<Type>(metric, time, coordinate, fixed), radius_(radius), radius_square_(radius * radius) {}
		int Hit(const boost::numeric::ublas::bounded_vector<Type, 8> &current, const boost::numeric::ublas::bounded_vector<Type, 8> &last) {
			Type a2 = this->metric_->DistanceSquare(this->position_, current, 3);
			if (a2 <= radius_square_)
				return 1;
			Type b2 = this->metric_->DistanceSquare(this->position_, last, 3), c2 = this->metric_->DistanceSquare(current, last, 3);
			if (a2 + c2 > b2 && b2 + c2 > a2 && 2 * a2 * c2 - Power2(a2 - b2 + c2) <= 4 * c2 * radius_square_)
				return 1; // if min distance between current and last < radius, return 1;
			return Status::SUCCESS;
		}
	};

	template <typename Type>
	class HotSpot : public Particle<Type> {
	  protected:
		const Type spectral_index_;
		const Type luminosity_;
		const Type luminosity_mu_;
		const Type luminosity_sigma_;

	  public:
		HotSpot(std::shared_ptr<Metric<Type>> metric, TimeSystem time, DynamicalSystem coordinate, Type spectral_index, Type luminosity, Type luminosity_mu, Type luminosity_sigma, bool fixed = false) : Particle<Type>(metric, time, coordinate, fixed), spectral_index_(spectral_index), luminosity_(luminosity), luminosity_mu_(luminosity_mu), luminosity_sigma_(luminosity_sigma) {}
		Type Luminosity(Type t) {
			Type intrinsic_time = this->time_ == T ? this->position_[0] : t;
			return luminosity_ * exp(-0.5 * Power2((intrinsic_time - luminosity_mu_) / luminosity_sigma_));
		}
		Type SpectralDensity(Type t, Type redshift) {
			Type intrinsic_time = this->time_ == T ? this->position_[0] : t;
			return luminosity_ * pow(redshift, spectral_index_) * exp(-0.5 * Power2((intrinsic_time - luminosity_mu_) / luminosity_sigma_));
		}
	};

	template <typename Type>
	class Disk : public Object<Type> {
	  protected:
		const Type inner_radius;
		const Type outer_radius;
		const Type inner_velocity;
		const Type velocity_slope_index;

	  public:
		int Hit(const boost::numeric::ublas::bounded_vector<Type, 8> &current, const boost::numeric::ublas::bounded_vector<Type, 8> &last) {
			if (!OppositeSign(current[2], last[2]))
				return Status::SUCCESS;
			if (std::abs(current[2]) < boost::math::constants::third_pi<Type>())
				return Status::SUCCESS;
			return Status::FAILURE;
		}
	};
} // namespace SBody

#endif
