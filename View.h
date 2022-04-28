#ifndef _VIEW_H
#define _VIEW_H

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "IO.h"
#include "Object.h"

namespace SBody {
	class view {
	  protected:
		const double r;
		const double theta;
		const double sinto;
		const double costo;
		const double tFinal;
		std::unique_ptr<IO::file> output;

	  public:
		view(double r, double theta, std::string fileName);
		int traceBack(Object::star &s, int rayNO);
		int shadow(int n);
	};
	class camera : public view {
	  protected:
		const size_t pixel;
		const double viewAngle;
		std::vector<std::array<double, 9>> initials;
		std::vector<std::vector<double>> screen;

	  public:
		camera(size_t pixel, double viewAngle, double r, double theta, std::string fileName);
		int traceBack();
		int lens();
		int save();
	};
} // namespace SBody

#endif
