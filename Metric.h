#ifndef _METRIC_H
#define _METRIC_H

namespace SBody {
	namespace Metric {
		// from cartesian to spherical
		int c2s(const double x[], const double v[], double r[], double w[]);
		// from spherical to cartesian
		int s2c(const double r[], const double w[], double x[], double v[]);
	} // namespace Metric
} // namespace SBody

#endif
