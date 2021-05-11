#include "Metric.h"

#include <cmath>

#include <gsl/gsl_math.h>

#include "Utility.h"

namespace SBody {
	namespace Metric {
		int c2s(const double x[], const double v[], double r[], double w[]) {
			// x = {x, y, z}
			// v = {v_x, v_y, v_z}
			// r = {r, \theta, \phi}
			// w = {v_r, v_\theta, v_\phi}
			r[0] = norm(x);
			if (r[0] == 0)
				return 1;
			r[1] = acos(x[2] / r[0]);
			w[0] = dot(x, v) / r[0];
			double normXY = norm(x, 2);
			if (normXY == 0) {
				double normVXY = norm(v, 2);
				if (normVXY == 0) {
					r[2] = 0;
					w[1] = 0;
				}
				else {
					if (v[1] >= 0)
						r[2] = acos(v[0] / normVXY);
					else
						r[2] = 2 * M_PI - acos(v[0] / normVXY);
					if (x[2] >= 0)
						w[1] = normVXY / r[0];
					else
						w[1] = -normVXY / r[0];
				}
				w[2] = 0;
			}
			else {
				if (x[1] >= 0)
					r[2] = acos(x[0] / normXY);
				else
					r[2] = 2 * M_PI - acos(x[0] / normXY);
				w[1] = (-v[2] + x[2] / r[0] * w[0]) / normXY;
				w[2] = (v[1] * x[0] - v[0] * x[1]) / gsl_pow_2(normXY);
			}
			return 0;
		}
		int s2c(const double r[], const double w[], double x[], double v[]) {
			// r = {r, \theta, \phi}
			// w = {v_r, v_\theta, v_\phi}
			// x = {x, y, z}
			// v = {v_x, v_y, v_z}
			x[0] = r[0] * sin(r[1]) * cos(r[2]);
			x[1] = r[0] * sin(r[1]) * sin(r[2]);
			x[2] = r[0] * cos(r[1]);
			v[0] = w[0] * sin(r[1]) * cos(r[2]) + r[0] * cos(r[1]) * cos(r[2]) * w[1] - r[0] * sin(r[1]) * sin(r[2]) * w[2];
			v[1] = w[0] * sin(r[1]) * sin(r[2]) + r[0] * cos(r[1]) * sin(r[2]) * w[1] + r[0] * sin(r[1]) * cos(r[2]) * w[2];
			v[2] = w[0] * cos(r[1]) - r[0] * sin(r[1]) * w[1];
			return 0;
		}
	} // namespace Metric
} // namespace SBody
