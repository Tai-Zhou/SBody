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

#ifndef SBODY_OBJECT_H
#define SBODY_OBJECT_H

#include <array>
#include <memory>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

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
		static std::vector<Object *> object_list_;
		/**
		 * @brief Construct a new Object object
		 *
		 * @param metric
		 */
		Object(std::shared_ptr<Metric<Type>> metric);
		/**
		 * @brief
		 *
		 * @param current
		 * @param last
		 * @return int
		 */
		virtual int Hit(const Type current[], const Type last[]) = 0;
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
		 * @param photon 8 dimensional information of photon
		 * @return double
		 */
		virtual Type Redshift(const Type photon[], TimeSystem photon_time) = 0;
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
		std::array<Type, 8> position_;
		/// if the position fixed
		const bool fixed_;
		/// integrator
		std::unique_ptr<Integrator> integrator_;

	  public:
		/**
		 * @brief Construct a new star object
		 *
		 * @param metric
		 * @param radius radius
		 * @param fixed whether the position of the star is fixed
		 */
		Particle(std::shared_ptr<Metric<Type>> metric, TimeSystem time, DynamicalSystem coordinate, bool fixed = false);
		int Position(Type *position);
		int Position(std::array<Type, 8> &position);
		int InitializeKeplerian(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type true_anomaly, Type observer_inclination = 0., Type observer_rotation = 0.);
		int InitializeGeodesic(Type orbital_radius, Type inclination, Type periapsis, Type ascending_node, Type v_r, Type v_phi, Type observer_inclination = 0., Type observer_rotation = 0.);
		int InitializeSchwarzschildKeplerianPericenter(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type observer_inclination, Type observer_rotation);
		int InitializeSchwarzschildKeplerianApocenter(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type observer_inclination, Type observer_rotation);
		int InitializeKeplerianHarmonic(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type true_anomaly, Type observer_inclination, Type observer_rotation);
		int InitializeCircular(Type r, Type phi, Type v_phi_ratio);
		int InitializeHelical(Type r, Type theta, Type phi, Type v_r, Type v_phi);
		int IntegratorApply(Type *t, Type t1, Type *h);
		int IntegratorApplyFixedStep(Type *t, const Type h);
		int IntegratorReset();

		/**
		 * @brief Check if the star hit by the photon
		 *
		 * @param current
		 * @param last
		 * @return int
		 */
		int Hit(const Type current[], const Type last[] = nullptr) override;

		/**
		 * @brief return frequency
		 *
		 * @param photon
		 * @return double
		 */
		Type Redshift(const Type photon[], TimeSystem photon_time) override;

		int MetricTensor(gsl_matrix *metric);
		Type DotProduct(const Type x[], const Type y[], const size_t dimension);
		int LocalInertialFrame(gsl_matrix *coordinate);
		Type Energy();
		Type AngularMomentum();
		Type CarterConstant();
	};

	template <typename Type>
	class Star : public Particle<Type> {
	  protected:
		/// radius
		const Type radius_;
		/// radius square
		const Type radius_square_;

	  public:
		Star(std::shared_ptr<Metric<Type>> metric, TimeSystem time, DynamicalSystem coordinate, Type radius = 0, bool fixed = false);
		int Hit(const Type current[], const Type last[]);
	};

	template <typename Type>
	class HotSpot : public Particle<Type> {
	  protected:
		const Type spectral_index_;
		const Type luminosity_;
		const Type luminosity_mu_;
		const Type luminosity_sigma_;

	  public:
		HotSpot(std::shared_ptr<Metric<Type>> metric, TimeSystem time, DynamicalSystem coordinate, Type spectral_index, Type luminosity, Type luminosity_mu, Type luminosity_sigma, bool fixed = false);
		Type Luminosity(Type t);
		Type SpectralDensity(Type t, Type redshift);
	};

	template <typename Type>
	class Disk : public Object<Type> {
	  protected:
		const Type inner_radius;
		const Type outer_radius;
		const Type inner_velocity;
		const Type velocity_slope_index;

	  public:
		int Hit(const Type current[], const Type last[]);
	};

	template <typename Type>
	std::vector<Object<Type> *> Object<Type>::object_list_;
	template <typename Type>
	Object<Type>::Object(std::shared_ptr<Metric<Type>> metric) : metric_(metric) {
		object_list_.push_back(this);
	}
	template <typename Type>
	Particle<Type>::Particle(std::shared_ptr<Metric<Type>> metric, TimeSystem time, DynamicalSystem coordinate, bool fixed) : Object<Type>(metric), time_(time), coordinate_(coordinate), fixed_(fixed) {
		position_.fill(0.);
	}
	template <typename Type>
	int Particle<Type>::Position(Type *position) {
		std::copy(position_.begin(), position_.end(), position);
		return Status::SUCCESS;
	}
	template <typename Type>
	int Particle<Type>::Position(std::array<Type, 8> &position) {
		std::copy(position_.begin(), position_.end(), position.begin());
		return Status::SUCCESS;
	}
	template <typename Type>
	int Particle<Type>::InitializeKeplerian(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type true_anomaly, Type observer_inclination, Type observer_rotation) {
		integrator_ = this->metric_->GetIntegrator(time_, coordinate_, GEODESIC);
		position_[0] = 0.;
		periapsis += true_anomaly;
		Type r = a * (1. - e * e) / (1. + e * cos(true_anomaly));
		const Type tp1 = -r * cos(periapsis), tp2 = -r * sin(periapsis) * cos(inclination);
		const Type xp1 = tp1 * cos(ascending_node) - tp2 * sin(ascending_node), xp2 = tp2 * cos(ascending_node) + tp1 * sin(ascending_node), xp3 = -r * sin(periapsis) * sin(inclination);
		position_[1] = (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * cos(observer_inclination) + xp3 * sin(observer_inclination);
		position_[2] = xp2 * cos(observer_rotation) - xp1 * sin(observer_rotation);
		position_[3] = xp3 * cos(observer_inclination) - (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * sin(observer_inclination);
		if (fixed_) {
			position_[5] = 0.;
			position_[6] = 0.;
			position_[7] = 0.;
		} else {
			const Type vphi = sqrt((1 - e * e) * a) / r, vr = GSL_SIGN(PhiDifference(true_anomaly)) * sqrt(std::max(0., 2. / r - 1. / a - vphi * vphi));
			const Type tp5 = vphi * sin(periapsis) - vr * cos(periapsis), tp6 = -(vphi * cos(periapsis) + vr * sin(periapsis)) * cos(inclination);
			const Type xp5 = tp5 * cos(ascending_node) - tp6 * sin(ascending_node), xp6 = tp5 * sin(ascending_node) + tp6 * cos(ascending_node), xp7 = -(vphi * cos(periapsis) + vr * sin(periapsis)) * sin(inclination);
			position_[5] = (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * cos(observer_inclination) + xp7 * sin(observer_inclination);
			position_[6] = xp6 * cos(observer_rotation) - xp5 * sin(observer_rotation);
			position_[7] = xp7 * cos(observer_inclination) - (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * sin(observer_inclination);
		}
		CartesianToSpherical(position_.data());
		return this->metric_->NormalizeTimelikeGeodesic(position_.data());
	}
	template <typename Type>
	int Particle<Type>::InitializeGeodesic(Type orbital_radius, Type inclination, Type periapsis, Type ascending_node, Type v_r, Type v_phi, Type observer_inclination, Type observer_rotation) {
		integrator_ = this->metric_->GetIntegrator(time_, coordinate_, GEODESIC);
		position_[0] = 0.;
		const Type tp1 = -orbital_radius * cos(periapsis), tp2 = -orbital_radius * sin(periapsis) * cos(inclination);
		const Type xp1 = tp1 * cos(ascending_node) - tp2 * sin(ascending_node), xp2 = tp2 * cos(ascending_node) + tp1 * sin(ascending_node), xp3 = -orbital_radius * sin(periapsis) * sin(inclination);
		position_[1] = (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * cos(observer_inclination) + xp3 * sin(observer_inclination);
		position_[2] = xp2 * cos(observer_rotation) - xp1 * sin(observer_rotation);
		position_[3] = xp3 * cos(observer_inclination) - (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * sin(observer_inclination);
		if (fixed_) {
			position_[5] = 0.;
			position_[6] = 0.;
			position_[7] = 0.;
		} else {
			const Type tp5 = v_phi * sin(periapsis) - v_r * cos(periapsis), tp6 = -(v_phi * cos(periapsis) + v_r * sin(periapsis)) * cos(inclination);
			const Type xp5 = tp5 * cos(ascending_node) - tp6 * sin(ascending_node), xp6 = tp5 * sin(ascending_node) + tp6 * cos(ascending_node), xp7 = -(v_phi * cos(periapsis) + v_r * sin(periapsis)) * sin(inclination);
			position_[5] = (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * cos(observer_inclination) + xp7 * sin(observer_inclination);
			position_[6] = xp6 * cos(observer_rotation) - xp5 * sin(observer_rotation);
			position_[7] = xp7 * cos(observer_inclination) - (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * sin(observer_inclination);
		}
		CartesianToSpherical(position_.data());
		return this->metric_->NormalizeTimelikeGeodesic(position_.data());
	}
	template <typename Type>
	int Particle<Type>::InitializeSchwarzschildKeplerianPericenter(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type observer_inclination, Type observer_rotation) {
		integrator_ = this->metric_->GetIntegrator(time_, coordinate_, GEODESIC);
		position_[0] = 0.;
		Type r = a * (1. - e); // pericenter
		const Type tp1 = -r * cos(periapsis), tp2 = -r * sin(periapsis) * cos(inclination);
		const Type xp1 = tp1 * cos(ascending_node) - tp2 * sin(ascending_node), xp2 = tp2 * cos(ascending_node) + tp1 * sin(ascending_node), xp3 = -r * sin(periapsis) * sin(inclination);
		position_[1] = (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * cos(observer_inclination) + xp3 * sin(observer_inclination);
		position_[2] = xp2 * cos(observer_rotation) - xp1 * sin(observer_rotation);
		position_[3] = xp3 * cos(observer_inclination) - (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * sin(observer_inclination);
		if (fixed_) {
			position_[5] = 0.;
			position_[6] = 0.;
			position_[7] = 0.;
		} else {
			const Type E = (a - 2.) / sqrt(a * (a - 3.));
			const Type vphi = sqrt(E * E * r * (r - 2.) - Power2(r - 2.)) / (r * E);
			const Type tp5 = vphi * sin(periapsis), tp6 = -vphi * cos(periapsis) * cos(inclination);
			const Type xp5 = tp5 * cos(ascending_node) - tp6 * sin(ascending_node), xp6 = tp5 * sin(ascending_node) + tp6 * cos(ascending_node), xp7 = -vphi * cos(periapsis) * sin(inclination);
			position_[5] = (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * cos(observer_inclination) + xp7 * sin(observer_inclination);
			position_[6] = xp6 * cos(observer_rotation) - xp5 * sin(observer_rotation);
			position_[7] = xp7 * cos(observer_inclination) - (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * sin(observer_inclination);
		}
		CartesianToSpherical(position_);
		return this->metric_->NormalizeTimelikeGeodesic(position_);
	}
	template <typename Type>
	int Particle<Type>::InitializeSchwarzschildKeplerianApocenter(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type observer_inclination, Type observer_rotation) {
		integrator_ = this->metric_->GetIntegrator(time_, coordinate_, GEODESIC);
		position_[0] = 0.;
		Type r = a * (1. + e); // apocenter
		const Type tp1 = r * cos(periapsis), tp2 = r * sin(periapsis) * cos(inclination);
		const Type xp1 = tp1 * cos(ascending_node) - tp2 * sin(ascending_node), xp2 = tp2 * cos(ascending_node) + tp1 * sin(ascending_node), xp3 = r * sin(periapsis) * sin(inclination);
		position_[1] = (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * cos(observer_inclination) + xp3 * sin(observer_inclination);
		position_[2] = xp2 * cos(observer_rotation) - xp1 * sin(observer_rotation);
		position_[3] = xp3 * cos(observer_inclination) - (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * sin(observer_inclination);
		if (fixed_) {
			position_[5] = 0.;
			position_[6] = 0.;
			position_[7] = 0.;
		} else {
			const Type E = (a - 2.) / sqrt(a * (a - 3.));
			const Type vphi = sqrt(E * E * r * (r - 2.) - Power2(r - 2.)) / (r * E);
			const Type tp5 = -vphi * sin(periapsis), tp6 = vphi * cos(periapsis) * cos(inclination);
			const Type xp5 = tp5 * cos(ascending_node) - tp6 * sin(ascending_node), xp6 = tp5 * sin(ascending_node) + tp6 * cos(ascending_node), xp7 = vphi * cos(periapsis) * sin(inclination);
			position_[5] = (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * cos(observer_inclination) + xp7 * sin(observer_inclination);
			position_[6] = xp6 * cos(observer_rotation) - xp5 * sin(observer_rotation);
			position_[7] = xp7 * cos(observer_inclination) - (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * sin(observer_inclination);
		}
		CartesianToSpherical(position_.data());
		return this->metric_->NormalizeTimelikeGeodesic(position_.data());
	}
	template <typename Type>
	int Particle<Type>::InitializeKeplerianHarmonic(Type a, Type e, Type inclination, Type periapsis, Type ascending_node, Type true_anomaly, Type observer_inclination, Type observer_rotation) {
		InitializeKeplerian(a, e, inclination, periapsis, ascending_node, true_anomaly, observer_inclination, observer_rotation);
		position_[1] += 1.;
		return this->metric_->NormalizeTimelikeGeodesic(position_.data());
	}
	template <typename Type>
	int Particle<Type>::InitializeCircular(Type r, Type phi, Type v_phi_ratio) {
		integrator_ = this->metric_->GetIntegrator(time_, coordinate_, CIRCULAR);
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
			position_[7] = v_phi_ratio / (r * sqrt(r));
		}
		return this->metric_->NormalizeTimelikeGeodesic(position_.data());
	}
	template <typename Type>
	int Particle<Type>::InitializeHelical(Type r, Type theta, Type phi, Type v_r, Type v_phi) {
		integrator_ = this->metric_->GetIntegrator(time_, coordinate_, HELICAL);
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
		return this->metric_->NormalizeTimelikeGeodesic(position_.data());
	}

	template <typename Type>
	int Particle<Type>::IntegratorApply(Type *t, Type t1, Type *h) {
		if (fixed_) {
			*t = t1;
			return Status::SUCCESS;
		}
		return integrator_->Apply(t, t1, h, position_.data());
	}
	template <typename Type>
	int Particle<Type>::IntegratorApplyFixedStep(Type *t, const Type h) {
		if (fixed_) {
			*t += h;
			return Status::SUCCESS;
		}
		return integrator_->ApplyFixedStep(t, h, position_);
	}
	template <typename Type>
	int Particle<Type>::IntegratorReset() {
		return integrator_->Reset();
	}
	template <typename Type>
	int Particle<Type>::Hit(const Type current[], const Type last[]) {
		return false;
	}
	template <typename Type>
	Type Particle<Type>::Redshift(const Type photon[], TimeSystem photon_time) {
		return this->metric_->Redshift(position_.data(), photon, time_, photon_time);
	}
	template <typename Type>
	int Particle<Type>::MetricTensor(gsl_matrix *metric) {
		return this->metric_->MetricTensor(position_, metric);
	}
	template <typename Type>
	Type Particle<Type>::DotProduct(const Type x[], const Type y[], const size_t dimension) {
		return this->metric_->DotProduct(position_, x, y, dimension);
	}
	template <typename Type>
	int Particle<Type>::LocalInertialFrame(gsl_matrix *coordinate) {
		return this->metric_->LocalInertialFrame(position_, time_, coordinate);
	}
	template <typename Type>
	Type Particle<Type>::Energy() {
		return this->metric_->Energy(position_.data(), time_, coordinate_);
	}
	template <typename Type>
	Type Particle<Type>::AngularMomentum() {
		return this->metric_->AngularMomentum(position_.data(), time_, coordinate_);
	}
	template <typename Type>
	Type Particle<Type>::CarterConstant() {
		return this->metric_->CarterConstant(position_.data(), 1., time_, coordinate_);
	}
	template <typename Type>
	Star<Type>::Star(std::shared_ptr<Metric<Type>> metric, TimeSystem time, DynamicalSystem coordinate, Type radius, bool fixed) : Particle<Type>(metric, time, coordinate, fixed), radius_(radius), radius_square_(radius * radius) {}
	template <typename Type>
	int Star<Type>::Hit(const Type current[], const Type last[]) {
		Type a2 = this->metric_->DistanceSquare(this->position_, current, 3);
		if (a2 <= radius_square_)
			return 1;
		Type b2 = this->metric_->DistanceSquare(this->position_, last, 3), c2 = this->metric_->DistanceSquare(current, last, 3);
		if (a2 + c2 > b2 && b2 + c2 > a2 && 2 * a2 * c2 - Power2(a2 - b2 + c2) <= 4 * c2 * radius_square_)
			return 1; // if min distance between current and last < radius, return 1;
		return Status::SUCCESS;
	}
	template <typename Type>
	HotSpot<Type>::HotSpot(std::shared_ptr<Metric<Type>> metric, TimeSystem time, DynamicalSystem coordinate, Type spectral_index, Type luminosity, Type luminosity_mu, Type luminosity_sigma, bool fixed) : Particle<Type>(metric, time, coordinate, fixed), spectral_index_(spectral_index), luminosity_(luminosity), luminosity_mu_(luminosity_mu), luminosity_sigma_(luminosity_sigma) {}
	template <typename Type>
	Type HotSpot<Type>::Luminosity(Type t) {
		Type intrinsic_time = this->time_ == T ? this->position_[0] : t;
		return luminosity_ * exp(-0.5 * Power2((intrinsic_time - luminosity_mu_) / luminosity_sigma_));
	}
	template <typename Type>
	Type HotSpot<Type>::SpectralDensity(Type t, Type redshift) {
		Type intrinsic_time = this->time_ == T ? this->position_[0] : t;
		return luminosity_ * pow(redshift, spectral_index_) * exp(-0.5 * Power2((intrinsic_time - luminosity_mu_) / luminosity_sigma_));
	}
	template <typename Type>
	int Disk<Type>::Hit(const Type current[], const Type last[]) {
		if (!OppositeSign(current[2], last[2]))
			return Status::SUCCESS;
		if (abs(current[2]) < M_PI_4)
			return Status::SUCCESS;
		return Status::FAILURE;
	}
} // namespace SBody

#endif
