#ifndef _VIEW_H
#define _VIEW_H

#include <array>
#include <vector>

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
		std::vector<std::array<double, 9>> initials;

	  public:
		view(size_t pixel, double viewAngle, double r, double theta, double phi, double tFinal, size_t duration = 300, size_t frame = 30);
		void traceBack();
	};
} // namespace SBody

#endif
