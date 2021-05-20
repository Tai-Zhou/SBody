#ifndef _VIEW_H
#define _VIEW_H

#include <array>
#include <vector>

namespace SBody {
	class view {
	  protected:
		const double viewAngle;
		const double r;
		const double theta;
		const double tFinal;
		const size_t duration;
		const size_t frame;

	  public:
		view(double viewAngle, double r, double theta, double tFinal, size_t duration = 300, size_t frame = 30);
		void traceBack();
	};
	class camera : public view {
	  protected:
		const size_t pixel;
		std::vector<std::array<double, 9>> initials;

	  public:
		camera(size_t pixel, double viewAngle, double r, double theta, double tFinal, size_t duration = 300, size_t frame = 30);
		void traceBack();
	};
} // namespace SBody

#endif
