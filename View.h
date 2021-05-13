#ifndef _VIEW_H
#define _VIEW_H

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
		std::vector<std::vector<std::tuple<double, double, double>>> traces;

	  public:
		view(const size_t _pixel, const double _viewAngle, const double _r, const double _theta, const double _phi, const double _tFinal, const size_t _duration = 300, const size_t _frame = 30);
		void traceBack(size_t metric, struct source *param);
	};
} // namespace SBody

#endif
