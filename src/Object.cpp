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

namespace SBody {
	std::vector<Object *> Object::object_list_;
	Object::Object(std::shared_ptr<Metric> metric) : metric_(metric) {
		object_list_.push_back(this);
	}
	Star::Star(std::shared_ptr<Metric> metric, time_system time, coordinate_system coordinate, double radius, bool fixed) : Object(metric), time_(time), coordinate_(coordinate), fixed_(fixed), radius_(radius), radius_square_(gsl_pow_2(radius)) {
		for (int i = 0; i < 8; ++i)
			position_[i] = 0;
	}
	int Star::Position(double *position) {
		std::copy(position_, position_ + 8, position);
		return GSL_SUCCESS;
	}
	int Star::InitializeKeplerian(double a, double e, double inclination, double periapsis, double ascending_node, double true_anomaly, double observer_inclination, double observer_rotation) {
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
			const double vphi = sqrt((1 - e * e) * a) / r, vr = GSL_SIGN(PhiDifference(true_anomaly)) * sqrt(std::max(0., 2. / r - 1. / a - vphi * vphi));
			const double tp5 = vphi * sin(periapsis) - vr * cos(periapsis), tp6 = -(vphi * cos(periapsis) + vr * sin(periapsis)) * cos(inclination);
			const double xp5 = tp5 * cos(ascending_node) - tp6 * sin(ascending_node), xp6 = tp5 * sin(ascending_node) + tp6 * cos(ascending_node), xp7 = -(vphi * cos(periapsis) + vr * sin(periapsis)) * sin(inclination);
			position_[5] = (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * cos(observer_inclination) + xp7 * sin(observer_inclination);
			position_[6] = xp6 * cos(observer_rotation) - xp5 * sin(observer_rotation);
			position_[7] = xp7 * cos(observer_inclination) - (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * sin(observer_inclination);
		}
		CartesianToSpherical(position_);
		return metric_->NormalizeTimelikeGeodesic(position_);
	}
	int Star::InitializeGeodesic(double orbital_radius, double inclination, double periapsis, double ascending_node, double v_r, double v_phi, double observer_inclination, double observer_rotation) {
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
	int Star::InitializeSchwarzschildKeplerianPericenter(double a, double e, double inclination, double periapsis, double ascending_node, double observer_inclination, double observer_rotation) {
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
	int Star::InitializeSchwarzschildKeplerianApocenter(double a, double e, double inclination, double periapsis, double ascending_node, double observer_inclination, double observer_rotation) {
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
			// const double vphi = sqrt(gsl_pow_2(E) * r * (r - 2.) - gsl_pow_2(r - 2.)) / (r * E);
			const double vphi = sqrt((1 - e * e) * a) / r;
			const double tp5 = -vphi * sin(periapsis), tp6 = vphi * cos(periapsis) * cos(inclination);
			const double xp5 = tp5 * cos(ascending_node) - tp6 * sin(ascending_node), xp6 = tp5 * sin(ascending_node) + tp6 * cos(ascending_node), xp7 = vphi * cos(periapsis) * sin(inclination);
			position_[5] = (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * cos(observer_inclination) + xp7 * sin(observer_inclination);
			position_[6] = xp6 * cos(observer_rotation) - xp5 * sin(observer_rotation);
			position_[7] = xp7 * cos(observer_inclination) - (xp5 * cos(observer_rotation) + xp6 * sin(observer_rotation)) * sin(observer_inclination);
		}
		CartesianToSpherical(position_);
		return metric_->NormalizeTimelikeGeodesic(position_);
	}
	int Star::InitializeKeplerianHarmonic(double a, double e, double inclination, double periapsis, double ascending_node, double true_anomaly, double observer_inclination, double observer_rotation) {
		InitializeKeplerian(a, e, inclination, periapsis, ascending_node, true_anomaly, observer_inclination, observer_rotation);
		position_[1] += 1.;
		return metric_->NormalizeTimelikeGeodesic(position_);
	}
	int Star::InitializeCircular(double r, double phi, double v_phi_ratio) {
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
	int Star::InitializeHelical(double r, double theta, double phi, double v_r, double v_phi) {
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
	int Star::IntegratorApply(double *t, double t1, double *h) {
		if (fixed_) {
			*t = t1;
			return GSL_SUCCESS;
		}
		return integrator_->Apply(t, t1, h, position_);
	}
	int Star::IntegratorApplyFixedStep(double *t, const double h) {
		if (fixed_) {
			*t += h;
			return GSL_SUCCESS;
		}
		return integrator_->ApplyFixedStep(t, h, position_);
	}
	int Star::IntegratorReset() {
		return integrator_->Reset();
	}
	int Star::Hit(const double current[], const double last[]) {
		double a2 = metric_->DistanceSquare(position_, current, 3);
		if (a2 <= radius_square_)
			return 1;
		double b2 = metric_->DistanceSquare(position_, last, 3), c2 = metric_->DistanceSquare(current, last, 3);
		if (a2 + c2 > b2 && b2 + c2 > a2 && 2 * a2 * c2 - gsl_pow_2(a2 - b2 + c2) <= 4 * c2 * radius_square_)
			return 1; // if min distance between current and last < radius, return 1;
		return 0;
	}
	double Star::Redshift(const double photon[]) {
		const double u[4] = {1., position_[5], position_[6], position_[7]}, v[4] = {1., photon[5], photon[6], photon[7]};
		return -metric_->DotProduct(position_, u, v, 4) / (position_[4] * photon[4]);
	}
	double Star::RedshiftTau(const double photon[]) {
		const double u[4] = {1., position_[5], position_[6], position_[7]};
		return -metric_->DotProduct(position_, u, photon + 4, 4) / position_[4];
	}
	int Star::MetricTensor(gsl_matrix *metric) {
		return metric_->MetricTensor(position_, metric);
	}
	double Star::DotProduct(const double x[], const double y[], const size_t dimension) {
		return metric_->DotProduct(position_, x, y, dimension);
	}
	int Star::LocalInertialFrame(gsl_matrix *coordinate) {
		return metric_->LocalInertialFrame(position_, coordinate);
	}
	double Star::Energy() {
		return metric_->Energy(position_, coordinate_);
	}
	double Star::AngularMomentum() {
		return metric_->AngularMomentum(position_, coordinate_);
	}
	double Star::CarterConstant() {
		return metric_->CarterConstant(position_, 1., coordinate_);
	}
} // namespace SBody
