#ifndef _NEWTON_H
#define _NEWTON_H

namespace newton {
	extern const int dimension;
	int c2s(const double x[], double r[]);
	int s2c(const double r[], double x[]);
	double energy(const double y[]);
	double angularMomentum(const double y[]);
	int function(double t, const double y[], double dydt[], void *params);
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
} // namespace newton

#endif
