#ifndef SBODY_UTILITY_H
#define SBODY_UTILITY_H

#include <string>

#include <gsl/gsl_odeiv2.h>

namespace SBody {
	extern double absAcc;
	extern double relAcc;
	constexpr double epsilon = 1e-10;
	class integrator {
	  private:
		const int coordinate;
		const gsl_odeiv2_step_type *type;
		gsl_odeiv2_control *control;
		gsl_odeiv2_evolve *evolve;
		gsl_odeiv2_step *step;
		gsl_odeiv2_system system;

	  public:
		integrator(int (*function)(double, const double *, double *, void *), int (*jacobian)(double, const double *, double *, double *, void *), int coordinate, void *params = nullptr, const gsl_odeiv2_step_type *type = gsl_odeiv2_step_rk8pd);
		~integrator();
		int apply(double *t, double t1, double *h, double *y);
		int apply_fixed(double *t, const double h, double *y);
		int reset();
		int checkCoordinate(double *y);
	};

	// Dot product of vector x·y, or x·x if y == nullptr
	double dot(const double x[], const double y[] = nullptr, size_t dimension = 3);

	// Length of vector x, with 3 dimensions set by default
	double norm(const double x[], size_t dimension = 3);

	// Cross product of vector x \times y, stored in z
	void cross(const double x[], const double y[], double z[]);

	// return sign of x
	int sign(int x);
	double sign(double x);

	// return 1 if x, y have opposite signs
	int oppositeSign(double x, double y);

	// return x in [0, 2*pi)
	double mod2Pi(double x);

	double _0x(double x);
	double _0x1(double x);
} // namespace SBody

#endif
