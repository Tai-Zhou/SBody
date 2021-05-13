#include "Utility.h"

#include <cmath>

#include <gsl/gsl_math.h>

#include "Constant.h"
#include "Metric.h"

namespace SBody {
	double absAcc = 1e-15, relAcc = 1e-15;
	integrator::integrator(const size_t NSK, void *const params, const gsl_odeiv2_step_type *type) : type(type), control(gsl_odeiv2_control_y_new(absAcc, relAcc)), evolve(gsl_odeiv2_evolve_alloc(8)), step(gsl_odeiv2_step_alloc(type, 8)) {
		system = {Metric::function[NSK], Metric::jacobian[NSK], 8, params};
	}
	int integrator::apply(double *t, const double t1, double *h, double *y) {
		return gsl_odeiv2_evolve_apply(evolve, control, step, &system, t, t1, h, y);
	}
	double dot(const double x[], const double y[], size_t dimension) {
		if (dimension == 3) {
			if (y == nullptr)
				return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
			return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
		}
		double sum = 0;
		while (dimension--)
			sum += x[dimension] * y[dimension];
		return sum;
	}
	double norm(const double x[], size_t dimension) {
		if (dimension == 3)
			return sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
		double sum = 0;
		while (dimension--)
			sum += x[dimension] * x[dimension];
		return sqrt(sum);
	}
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
} // namespace SBody
