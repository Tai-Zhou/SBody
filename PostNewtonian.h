#ifndef _POSTNEWTONIAN_H
#define _POSTNEWTONIAN_H

namespace postnewtonian {
	extern const int dimension;
	extern int PN;
	int c2s(const double x[], double r[]);
	int s2c(const double r[], double x[]);
	int function(double t, const double y[], double dydt[], void *params);
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
	double energy(const double y[], void *params);
	double angularMomentum(const double y[], void *params);
} // namespace postnewtonian

#endif
