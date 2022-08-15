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

#include <gsl/gsl_math.h>

#include "Metric.h"
#include "Utility.h"

namespace SBody {
	std::vector<Object *> Object::object_list_;
	Star::Star(double radius, const double position[], bool fixed) : fixed_(fixed), radius_(radius), radius_square_(gsl_pow_2(radius)) {
		std::copy(position, position + 8, pos);
	}
	int Star::Hit(const double current[], const double last[]) {
		double a2 = metric::ds2(pos, current, 3);
		if (a2 <= radius_square_)
			return 1;
		double b2 = metric::ds2(pos, last, 3), c2 = metric::ds2(current, last, 3);
		if (a2 + c2 > b2 && b2 + c2 > a2 && 2 * a2 * c2 - gsl_pow_2(a2 - b2 + c2) <= 4 * c2 * radius_square_)
			return 1; // if min distance between current and last < radius, return 1;
		return 0;
	}
	double Star::Redshift(const double photon[]) {
		const double u[4] = {1., pos[5], pos[6], pos[7]}, v[4] = {1., photon[5], photon[6], photon[7]};
		return -metric::dot(pos, u, v, 4) / (pos[4] * photon[4]);
	}
	double Star::RedshiftTau(const double photon[]) {
		const double u[4] = {1., pos[5], pos[6], pos[7]};
		return -metric::dot(pos, u, photon + 4, 4) / pos[4];
	}
	Disk::Disk(double inner_radius, double outer_radius) : inner_radius_(inner_radius), outer_radius_(outer_radius) {}
	int Disk::Hit(const double current[], const double last[]) {
		if (inner_radius_ <= current[1] && current[1] <= outer_radius_ && OppositeSign(current[2] - M_PI_2, last[2] - M_PI_2))
			return 1;
		// if min distance between current and last < radius, return 1;
		return 0;
	}
	double Disk::Redshift(const double ph[]) {
		return 1.;
	}
	ThickDisk::ThickDisk(double inner_radius, double outer_radius, double half_angle) : Disk(inner_radius, outer_radius), half_angle_(half_angle) {}
	int ThickDisk::Hit(const double current[], const double last[]) { // TODO: need update
		if (inner_radius_ <= current[1] && current[1] <= outer_radius_ && OppositeSign(current[2] - M_PI_2, last[2] - M_PI_2))
			return 1;
		// if min distance between current and last < radius, return 1;
		return 0;
	}
	double ThickDisk::Redshift(const double ph[]) {
		return 1.;
	}
	Torus::Torus(double major_radius, double minor_radius) : major_radius_(major_radius), minor_radius_(minor_radius) {}
	int Torus::Hit(const double current[], const double last[]) {
		const double rc[4] = {0, major_radius_, M_PI_2, current[3]};
		if (metric::ds2(rc, current, 3) <= minor_radius_) // FIXME:not accurate if minorRadius not << majorRadius
			return 1;
		// if min distance between current and last < radius, return 1;
		return 0;
	}
	double Torus::Redshift(const double ph[]) {
		return 1.;
	}
} // namespace SBody
