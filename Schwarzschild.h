#ifndef _SCHWARZSCHILD_H
#define _SCHWARZSCHILD_H

namespace schwarzschild {
	extern const int dimension;
	int c2s(const double x[], double r[]);
	int s2c(const double r[], double x[]);
	int function(double t, const double y[], double dydt[], void *params);
	int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
	double energy(const double y[]);
	double angularMomentum(const double y[]);
	namespace particle {
		int normalization(double y[], void *params = nullptr);
	} // namespace particle
	namespace light {
		int normalization(double y[], void *params = nullptr);
	} // namespace light
} // namespace schwarzschild

#endif
