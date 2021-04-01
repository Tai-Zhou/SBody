#include "Schwarzschild.h"

#include <cmath>

#include <gsl/gsl_errno.h>

#include "Constant.h"
#include "Utility.h"

namespace kerr {
	namespace particle {
		int function(double t, const double y[], double dydt[], void *params);
		int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
	} // namespace particle
	namespace light {
		int function(double t, const double y[], double dydt[], void *params);
		int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
	} // namespace light
} // namespace kerr
