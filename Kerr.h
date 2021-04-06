#ifndef _KERR_H
#define _KERR_H

namespace kerr {
	extern const int dimension;
	int c2s(const double x[], double r[]);
	int s2c(const double r[], double x[]);
	int function(double t, const double y[], double dydt[], void *params);
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
	double energy(const double y[], void *params);
	double angularMomentum(const double y[], void *params);
	namespace particle {
		int normalization(double y[], void *params);
	} // namespace particle
	namespace light {
		int normalization(double y[], void *params);
	} // namespace light
} // namespace kerr

#endif
