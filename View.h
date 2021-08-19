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
		const size_t duration;
		const size_t frame;
		std::vector<std::vector<double>> screen;

	  public:
		view(double r, double theta, double tFinal, size_t duration = 300, size_t frame = 30);
		void traceBack(Object::star &s, int rayNO);
		void save(std::string fileName = "view");
	};
	class camera : public view {
	  protected:
		const size_t pixel;
		const double viewAngle;
		std::vector<std::array<double, 9>> initials;

	  public:
		camera(size_t pixel, double viewAngle, double r, double theta, double tFinal, size_t duration = 300, size_t frame = 30);
		void initialize();
		void traceBack();
	};
} // namespace SBody

#endif
