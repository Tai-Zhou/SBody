#include "Object.h"

#include <gsl/gsl_math.h>

#include "Metric.h"
#include "Utility.h"

namespace SBody {
	namespace Object {
		std::vector<object *> objectList;
		star::star(double radius, const double position[], int fixed) : fixed(fixed), radius(radius), r2(gsl_pow_2(radius)) {
			for (int i = 0; i < 8; ++i)
				pos[i] = position[i];
		}
		int star::hit(const double current[], const double last[]) {
			double a2 = Metric::ds2(pos, current, 3);
			if (a2 <= r2)
				return 1;
			double b2 = Metric::ds2(pos, last, 3), c2 = Metric::ds2(current, last, 3);
			if (a2 + c2 > b2 && b2 + c2 > a2 && 2 * a2 * c2 - gsl_pow_2(a2 - b2 + c2) <= 4 * c2 * r2)
				return 1; // if min distance between current and last < radius, return 1;
			return 0;
		}
		double star::frequency(const double ph[]) {
			const double u[4] = {1, pos[5], pos[6], pos[7]}, v[4] = {1, -ph[5], -ph[6], -ph[7]};
			return -Metric::dot(pos, u, v, 4) / (pos[4] * ph[4]);
		}
		disk::disk(double innerRadius, double outerRadius) : innerRadius(innerRadius), outerRadius(outerRadius) {}
		int disk::hit(const double current[], const double last[]) {
			if (innerRadius <= current[1] && current[1] <= outerRadius && oppositeSign(current[2] - M_PI_2, last[2] - M_PI_2))
				return 1;
			// if min distance between current and last < radius, return 1;
			return 0;
		}
		double disk::frequency(const double ph[]) {
			return 1.;
		}
		thickDisk::thickDisk(double innerRadius, double outerRadius, double halfAngle) : disk(innerRadius, outerRadius), halfAngle(halfAngle) {}
		int thickDisk::hit(const double current[], const double last[]) { //TODO: need update
			if (innerRadius <= current[1] && current[1] <= outerRadius && oppositeSign(current[2] - M_PI_2, last[2] - M_PI_2))
				return 1;
			// if min distance between current and last < radius, return 1;
			return 0;
		}
		double thickDisk::frequency(const double ph[]) {
			return 1.;
		}
		torus::torus(double majorRadius, double minorRadius) : majorRadius(majorRadius), minorRadius(minorRadius) {}
		int torus::hit(const double current[], const double last[]) {
			const double rc[4] = {0, majorRadius, M_PI_2, current[3]};
			if (Metric::ds2(rc, current, 3) <= minorRadius) //FIXME:not accurate if minorRadius not << majorRadius
				return 1;
			// if min distance between current and last < radius, return 1;
			return 0;
		}
		double torus::frequency(const double ph[]) {
			return 1.;
		}
	} // namespace Object
} // namespace SBody
