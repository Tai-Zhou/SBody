#include "Utility.h"

#include <cmath>

#include <gsl/gsl_math.h>

#include "Constant.h"

namespace SBody {
	// Dot product of vector x·y, or x·x if y == nullptr
	double dot(const double x[], const double y[], int dimension) {
		if (dimension == 3) {
			if (y == NULL)
				return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
			return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
		}
		double sum = 0;
		while (dimension--)
			sum += x[dimension] * y[dimension];
		return sum;
	}
	// Length of vector x, with 3 dimensions set by default
	double norm(const double x[], int dimension) {
		if (dimension == 3)
			return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
		double sum = 0;
		while (dimension--)
			sum += x[dimension] * x[dimension];
		return sqrt(sum);
	}
	// Cross product of vector x✖y, stored in z
	void cross(const double x[], const double y[], double z[]) {
		z[0] = x[1] * y[2] - x[2] * y[1];
		z[1] = x[2] * y[0] - x[0] * y[2];
		z[2] = x[0] * y[1] - x[1] * y[0];
	}
	double _0x(const double x) {
		return x > 0 ? x : 0;
	}
	double _0x1(const double x) {
		if (x < 0)
			return 0;
		return x < 1 ? x : 1;
	}
}
