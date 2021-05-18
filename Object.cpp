#include "Object.h"

#include <gsl/gsl_math.h>

#include "Metric.h"
#include "Utility.h"

namespace SBody {
	namespace Object {
		star::star(double radius, const double position[], int fixed) : fixed(fixed), radius(radius) {
			for (int i = 0; i < 8; ++i)
				pos[i] = position[i];
		}
		int star::hit(const double current[], const double last[]) {
			if (Metric::ds2(pos, current) <= radius)
				return 1;
			// if min distance between current and last < radius, return 1;
			return 0;
		}
		int disk::hit(const double current[], const double last[]) {
			if (innerRadius <= current[1] && current[1] <= outerRadius && oppositeSign(current[2] - M_PI_2, last[2] - M_PI_2))
				return 1;
			// if min distance between current and last < radius, return 1;
			return 0;
		}
		thickDisk::thickDisk(double innerRadius, double outerRadius, double halfAngle) : disk(innerRadius, outerRadius), halfAngle(halfAngle) {}
		int torus::hit(const double current[], const double last[]) {
			const double rc[4] = {0, majorRadius, M_PI_2, current[3]};
			if (Metric::ds2(rc, current) <= minorRadius) //FIXME:not accurate if minorRadius not << majorRadius
				return 1;
			// if min distance between current and last < radius, return 1;
			return 0;
		}
	} // namespace Object
} // namespace SBody
