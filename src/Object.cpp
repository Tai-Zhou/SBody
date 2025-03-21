/**
 * @file Object.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Object.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include "Metric.h"
#include "Utility.h"

using namespace std;

namespace SBody {
	vector<Object *> Object::object_list_;
	Object::Object(shared_ptr<Metric> metric) : metric_(metric) {
		object_list_.push_back(this);
	}
	Particle::Particle(shared_ptr<Metric> metric, TimeSystem time, DynamicalSystem coordinate, bool fixed) : Object(metric), time_(time), coordinate_(coordinate), fixed_(fixed) {
		for (int i = 0; i < 8; ++i)
			position_[i] = 0;
	}
	int Particle::Position(double *position) {
		copy(position_, position_ + 8, position);
		return GSL_SUCCESS;
	}
	int Particle::InitializeKeplerian(double a, double e, double inclination, double periapsis, double ascending_node, double true_anomaly, double observer_inclination, double observer_rotation) {
		integrator_ = metric_->GetIntegrator(time_, coordinate_, GEODESIC);
		position_[0] = 0.;
		periapsis += true_anomaly;
		double r = a * (1. - e * e) / (1. + e * cos(true_anomaly));
		const double tp1 = -r * cos(periapsis), tp2 = -r * sin(periapsis) * cos(inclination);
		const double xp1 = tp1 * cos(ascending_node) - tp2 * sin(ascending_node), xp2 = tp2 * cos(ascending_node) + tp1 * sin(ascending_node), xp3 = -r * sin(periapsis) * sin(inclination);
		position_[1] = (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * cos(observer_inclination) + xp3 * sin(observer_inclination);
		position_[2] = xp2 * cos(observer_rotation) - xp1 * sin(observer_rotation);
		position_[3] = xp3 * cos(observer_inclination) - (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * sin(observer_inclination);
		if (fixed_) {
			position_[5] = 0.;
			position_[6] = 0.;
			position_[7] = 0.;
		} else {
			const double vphi = sqrt((1 - e * e) * a) / r, vr = GSL_SIGN(PhiDifference(true_anomaly)) * sqrt(max(0., 2. / r - 1. / a - vphi * vphi));
			const double tp5 = vphi * sin(periapsis) - vr * cos(periapsis), tp6 = -(vphi * cos(periapsis) + vr * sin(periapsis)) * cos(inclination);
			const double xp5 = tp5 * cos(ascending_node) - tp6 * sin(ascending_node), xp6 = tp5 * sin(ascending_node) + tp6 * cos(ascending_node), xp7 = -(vphi * cos(periapsis) + vr * sin(periapsis)) * sin(inclination);
			position_[5] = (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * cos(observer_inclination) + xp7 * sin(observer_inclination);
			position_[6] = xp6 * cos(observer_rotation) - xp5 * sin(observer_rotation);
			position_[7] = xp7 * cos(observer_inclination) - (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * sin(observer_inclination);
		}
		CartesianToSpherical(position_);
		return metric_->NormalizeTimelikeGeodesic(position_);
	}
	int Particle::InitializeGeodesic(double orbital_radius, double inclination, double periapsis, double ascending_node, double v_r, double v_phi, double observer_inclination, double observer_rotation) {
		integrator_ = metric_->GetIntegrator(time_, coordinate_, GEODESIC);
		position_[0] = 0.;
		const double tp1 = -orbital_radius * cos(periapsis), tp2 = -orbital_radius * sin(periapsis) * cos(inclination);
		const double xp1 = tp1 * cos(ascending_node) - tp2 * sin(ascending_node), xp2 = tp2 * cos(ascending_node) + tp1 * sin(ascending_node), xp3 = -orbital_radius * sin(periapsis) * sin(inclination);
		position_[1] = (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * cos(observer_inclination) + xp3 * sin(observer_inclination);
		position_[2] = xp2 * cos(observer_rotation) - xp1 * sin(observer_rotation);
		position_[3] = xp3 * cos(observer_inclination) - (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * sin(observer_inclination);
		if (fixed_) {
			position_[5] = 0.;
			position_[6] = 0.;
			position_[7] = 0.;
		} else {
			const double tp5 = v_phi * sin(periapsis) - v_r * cos(periapsis), tp6 = -(v_phi * cos(periapsis) + v_r * sin(periapsis)) * cos(inclination);
			const double xp5 = tp5 * cos(ascending_node) - tp6 * sin(ascending_node), xp6 = tp5 * sin(ascending_node) + tp6 * cos(ascending_node), xp7 = -(v_phi * cos(periapsis) + v_r * sin(periapsis)) * sin(inclination);
			position_[5] = (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * cos(observer_inclination) + xp7 * sin(observer_inclination);
			position_[6] = xp6 * cos(observer_rotation) - xp5 * sin(observer_rotation);
			position_[7] = xp7 * cos(observer_inclination) - (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * sin(observer_inclination);
		}
		CartesianToSpherical(position_);
		return metric_->NormalizeTimelikeGeodesic(position_);
	}
	int Particle::InitializeSchwarzschildKeplerianPericenter(double a, double e, double inclination, double periapsis, double ascending_node, double observer_inclination, double observer_rotation) {
		integrator_ = metric_->GetIntegrator(time_, coordinate_, GEODESIC);
		position_[0] = 0.;
		double r = a * (1. - e); // pericenter
		const double tp1 = -r * cos(periapsis), tp2 = -r * sin(periapsis) * cos(inclination);
		const double xp1 = tp1 * cos(ascending_node) - tp2 * sin(ascending_node), xp2 = tp2 * cos(ascending_node) + tp1 * sin(ascending_node), xp3 = -r * sin(periapsis) * sin(inclination);
		position_[1] = (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * cos(observer_inclination) + xp3 * sin(observer_inclination);
		position_[2] = xp2 * cos(observer_rotation) - xp1 * sin(observer_rotation);
		position_[3] = xp3 * cos(observer_inclination) - (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * sin(observer_inclination);
		if (fixed_) {
			position_[5] = 0.;
			position_[6] = 0.;
			position_[7] = 0.;
		} else {
			const double E = (a - 2.) / sqrt(a * (a - 3.));
			const double vphi = sqrt(gsl_pow_2(E) * r * (r - 2.) - gsl_pow_2(r - 2.)) / (r * E);
			const double tp5 = vphi * sin(periapsis), tp6 = -vphi * cos(periapsis) * cos(inclination);
			const double xp5 = tp5 * cos(ascending_node) - tp6 * sin(ascending_node), xp6 = tp5 * sin(ascending_node) + tp6 * cos(ascending_node), xp7 = -vphi * cos(periapsis) * sin(inclination);
			position_[5] = (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * cos(observer_inclination) + xp7 * sin(observer_inclination);
			position_[6] = xp6 * cos(observer_rotation) - xp5 * sin(observer_rotation);
			position_[7] = xp7 * cos(observer_inclination) - (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * sin(observer_inclination);
		}
		CartesianToSpherical(position_);
		return metric_->NormalizeTimelikeGeodesic(position_);
	}
	int Particle::InitializeSchwarzschildKeplerianApocenter(double a, double e, double inclination, double periapsis, double ascending_node, double observer_inclination, double observer_rotation) {
		integrator_ = metric_->GetIntegrator(time_, coordinate_, GEODESIC);
		position_[0] = 0.;
		double r = a * (1. + e); // apocenter
		const double tp1 = r * cos(periapsis), tp2 = r * sin(periapsis) * cos(inclination);
		const double xp1 = tp1 * cos(ascending_node) - tp2 * sin(ascending_node), xp2 = tp2 * cos(ascending_node) + tp1 * sin(ascending_node), xp3 = r * sin(periapsis) * sin(inclination);
		position_[1] = (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * cos(observer_inclination) + xp3 * sin(observer_inclination);
		position_[2] = xp2 * cos(observer_rotation) - xp1 * sin(observer_rotation);
		position_[3] = xp3 * cos(observer_inclination) - (xp1 * cos(observer_rotation) + xp2 * sin(observer_rotation)) * sin(observer_inclination);
		if (fixed_) {
			position_[5] = 0.;
			position_[6] = 0.;
			position_[7] = 0.;
		} else {
			const double E = (a - 2.) / sqrt(a * (a - 3.));
			const double vphi = sqrt(gsl_pow_2(E) * r * (r - 2.) - gsl_pow_2(r - 2.)) / (r * E);
			const double tp5 = -vphi * sin(periapsis), tp6 = vphi * cos(periapsis) * cos(inclination);
			const double xp5 = tp5 * cos(ascending_node) - tp6 * sin(ascending_node), xp6 = tp5 * sin(ascending_node) + tp6 * cos(ascending_node), xp7 = vphi * cos(periapsis) * sin(inclination);
			position_[5] = (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * cos(observer_inclination) + xp7 * sin(observer_inclination);
			position_[6] = xp6 * cos(observer_rotation) - xp5 * sin(observer_rotation);
			position_[7] = xp7 * cos(observer_inclination) - (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * sin(observer_inclination);
		}
		CartesianToSpherical(position_);
		return metric_->NormalizeTimelikeGeodesic(position_);
	}
	int Particle::InitializeKeplerianHarmonic(double a, double e, double inclination, double periapsis, double ascending_node, double true_anomaly, double observer_inclination, double observer_rotation) {
		InitializeKeplerian(a, e, inclination, periapsis, ascending_node, true_anomaly, observer_inclination, observer_rotation);
		position_[1] += 1.;
		return metric_->NormalizeTimelikeGeodesic(position_);
	}
	int Particle::InitializeCircular(double r, double phi, double v_phi_ratio) {
		integrator_ = metric_->GetIntegrator(time_, coordinate_, CIRCULAR);
		position_[0] = 0.;
		position_[1] = r;
		position_[2] = M_PI_2;
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
		return metric_->NormalizeTimelikeGeodesic(position_);
	}
	int Particle::InitializeHelical(double r, double theta, double phi, double v_r, double v_phi) {
		integrator_ = metric_->GetIntegrator(time_, coordinate_, HELICAL);
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
		return metric_->NormalizeTimelikeGeodesic(position_);
	}
	int Particle::IntegratorApply(double *t, double t1, double *h) {
		if (fixed_) {
			*t = t1;
			return GSL_SUCCESS;
		}
		return integrator_->Apply(t, t1, h, position_);
	}
	int Particle::IntegratorApplyFixedStep(double *t, const double h) {
		if (fixed_) {
			*t += h;
			return GSL_SUCCESS;
		}
		return integrator_->ApplyFixedStep(t, h, position_);
	}
	int Particle::IntegratorReset() {
		return integrator_->Reset();
	}
	int Particle::Hit(const double current[], const double last[]) {
		return false;
	}
	double Particle::Redshift(const double photon[], TimeSystem photon_time) {
		return metric_->Redshift(position_, photon, time_, photon_time);
	}
	int Particle::MetricTensor(gsl_matrix *metric) {
		return metric_->MetricTensor(position_, metric);
	}
	double Particle::DotProduct(const double x[], const double y[], const size_t dimension) {
		return metric_->DotProduct(position_, x, y, dimension);
	}
	int Particle::LocalInertialFrame(gsl_matrix *coordinate) {
		return metric_->LocalInertialFrame(position_, time_, coordinate);
	}
	double Particle::Energy() {
		return metric_->Energy(position_, time_, coordinate_);
	}
	double Particle::AngularMomentum() {
		return metric_->AngularMomentum(position_, time_, coordinate_);
	}
	double Particle::CarterConstant() {
		return metric_->CarterConstant(position_, 1., time_, coordinate_);
	}
	Star::Star(shared_ptr<Metric> metric, TimeSystem time, DynamicalSystem coordinate, double radius, bool fixed) : Particle(metric, time, coordinate, fixed), radius_(radius), radius_square_(gsl_pow_2(radius)) {}
	int Star::Hit(const double current[], const double last[]) {
		double a2 = metric_->DistanceSquare(position_, current, 3);
		if (a2 <= radius_square_)
			return 1;
		double b2 = metric_->DistanceSquare(position_, last, 3), c2 = metric_->DistanceSquare(current, last, 3);
		if (a2 + c2 > b2 && b2 + c2 > a2 && 2 * a2 * c2 - gsl_pow_2(a2 - b2 + c2) <= 4 * c2 * radius_square_)
			return 1; // if min distance between current and last < radius, return 1;
		return 0;
	}
	HotSpot::HotSpot(shared_ptr<Metric> metric, TimeSystem time, DynamicalSystem coordinate, double spectral_index, double luminosity, double luminosity_mu, double luminosity_sigma, bool fixed) : Particle(metric, time, coordinate, fixed), spectral_index_(spectral_index), luminosity_(luminosity), luminosity_mu_(luminosity_mu), luminosity_sigma_(luminosity_sigma) {}
	double HotSpot::Luminosity(double t) {
		double intrinsic_time = time_ == T ? position_[0] : t;
		return luminosity_ * exp(-0.5 * gsl_pow_2((intrinsic_time - luminosity_mu_) / luminosity_sigma_));
	}
	double HotSpot::SpectralDensity(double t, double redshift) {
		double intrinsic_time = time_ == T ? position_[0] : t;
		return luminosity_ * pow(redshift, spectral_index_) * exp(-0.5 * gsl_pow_2((intrinsic_time - luminosity_mu_) / luminosity_sigma_));
	}
	int Disk::Hit(const double current[], const double last[]) {
		if (!OppositeSign(current[2], last[2]))
			return 0;
		if (abs(current[2]) < M_PI_4)
			return 0;
		return 1;
	}
} // namespace SBody
