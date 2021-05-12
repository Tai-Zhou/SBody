#ifndef _UTILITY_H
#define _UTILITY_H

namespace SBody {
	// Dot product of vector x·y, or x·x if y == nullptr
	double dot(const double x[], const double y[] = nullptr, int dimension = 3);

	// Length of vector x, with 3 dimensions set by default
	double norm(const double x[], int dimension = 3);

	// Cross product of vector x \times y, stored in z
	void cross(const double x[], const double y[], double z[]);
	double _0x(const double x);
	double _0x1(const double x);
} // namespace SBody

#endif
