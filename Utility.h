#ifndef _UTILITY_H
#define _UTILITY_H

#include <gsl/gsl_odeiv2.h>

namespace SBody {
	extern double absAcc;
	extern double relAcc;
	struct integrator {
		const gsl_odeiv2_step_type *type;
		gsl_odeiv2_control *const control;
		gsl_odeiv2_evolve *const evolve;
		gsl_odeiv2_step *const step;
		gsl_odeiv2_system system;
		integrator(const int _metric, void *_params, const gsl_odeiv2_step_type *_type = gsl_odeiv2_step_rk8pd);
		int apply(double *t, double t1, double *h, double *y);
	};

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
