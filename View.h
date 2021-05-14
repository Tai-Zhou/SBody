#ifndef _VIEW_H
#define _VIEW_H

#include <array>
#include <vector>

#include "Utility.h"

namespace SBody {
	class view {
	  private:
		const size_t pixel;
		const double viewAngle;
		const double r;
		const double theta;
		const double phi;
		const double tFinal;
		const size_t duration;
		const size_t frame;
		std::vector<std::vector<std::array<double, 3>>> traces;

	  public:
		view(const size_t pixel, const double viewAngle, const double r, const double theta, const double phi, const double tFinal, const size_t duration = 300, const size_t frame = 30);
		void traceBack(size_t NSK, struct source *params);
	};
} // namespace SBody

#endif
