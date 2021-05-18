#include "Object.h"

#include "Metric.h"

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
	} // namespace Object
} // namespace SBody
