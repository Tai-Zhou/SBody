#ifndef _VIEW_H
#define _VIEW_H

#include <array>
#include <string>
#include <vector>

#include "Object.h"

namespace SBody {
	class view {
	  protected:
		const double r;
		const double theta;
		const double sinto;
		const double costo;
		const double tFinal;
		std::vector<std::vector<double>> screen;

	  public:
		view(double r, double theta);
		int traceBack(Object::star &s, int rayNO);
		int shadow(int n);
		int save(std::string fileName = "view");
	};
	class camera : public view {
	  protected:
		const size_t pixel;
		const double viewAngle;
		std::vector<std::array<double, 9>> initials;

	  public:
		camera(size_t pixel, double viewAngle, double r, double theta);
		int initialize();
		int traceBack();
	};
} // namespace SBody

#endif
