#include "View.h"

#include <array>
#include <cmath>
#include <tuple>

#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>

#include "Constant.h"
#include "Metric.h"

using namespace std;

namespace SBody {
	view::view(const size_t pixel, const double viewAngle, const double r, const double theta, const double phi, const double tFinal, const size_t duration, const size_t frame)
		: pixel(pixel),
		  viewAngle(viewAngle),
		  r(r),
		  theta(theta),
		  phi(phi),
		  tFinal(tFinal),
		  duration(duration),
		  frame(frame) {
		traces = vector<vector<tuple<double, double, double>>>(pixel * pixel, {make_tuple(r, theta, phi)});
	}
	void view::traceBack(size_t NSK, source *params) {
		const double sint = sin(theta), tana_pix = 2. * tan(0.5 * viewAngle) / pixel, tFrame = tFinal / duration;
		double ph[8] = {0., r, theta, phi, 1., -1., 0., 0.}, tStep = 0, t;
		int status = 0;
		integrator integ(NSK, params);
		if (theta < Constant::epsilon || M_PI - theta < Constant::epsilon) {
			for (int i = 0; i < pixel; ++i)
				for (int j = 0; j < pixel; ++j) {
					const double k = gsl_hypot(i - 0.5 * pixel, j - 0.5 * pixel);
					ph[6] = (theta < M_PI_2 ? 1 : -1) * tana_pix * k;
					if (k < Constant::epsilon)
						ph[3] = 0;
					else {
						if (2 * i <= pixel)
							ph[3] = (theta < M_PI_2 ? 1. : -1.) * acos((j - 0.5 * pixel) / k);
						else
							ph[3] = (theta < M_PI_2 ? 1. : -1.) * (2. * M_PI - acos((j - 0.5 * pixel) / k));
					}
					if (NSK == 1) {
						Metric::Schwarzschild::lightNormalization(ph, params);
						for (t = 0; t < tFinal; t += tFrame) {
							for (int k = 0; k < 4; ++k)
								ph[k] += tFrame * ph[k + 4];
							traces[i * pixel + j].push_back(make_tuple(ph[1], ph[2], ph[3]));
						}
					}
					else if (NSK == 2) {
						Metric::Kerr::lightNormalization(ph, params);
						for (t = 0; t < tFinal; t += tFrame) {
							for (int k = 0; k < 4; ++k)
								ph[k] += tFrame * ph[k + 4];
							traces[i * pixel + j].push_back(make_tuple(ph[1], ph[2], ph[3]));
						}
					}
				}
		}
		else
			for (int i = 0; i < pixel; ++i) {
				ph[6] = tana_pix * (i - 0.5 * pixel);
				for (int j = 0; j < pixel; ++j) {
					ph[7] = tana_pix * (j - 0.5 * pixel) / sint;
					if (NSK == 1) {
						Metric::Schwarzschild::lightNormalization(ph, params);
						for (t = 0; t < tFinal; t += tFrame) {
							for (int k = 0; k < 4; ++k)
								ph[k] += tFrame * ph[k + 4];
							traces[i * pixel + j].push_back(make_tuple(ph[1], ph[2], ph[3]));
						}
					}
					else if (NSK == 2) {
						Metric::Kerr::lightNormalization(ph, params);
						for (t = 0; t < tFinal; t += tFrame) {
							for (int k = 0; k < 4; ++k)
								ph[k] += tFrame * ph[k + 4];
							traces[i * pixel + j].push_back(make_tuple(ph[1], ph[2], ph[3]));
						}
					}
				}
			}
		/*if (NSK == 0) {
		}
		else if (NSK == 1) {
		}
		else if (NSK == 2) {
		}
		else if (NSK == 3) {
		}*/
	}
} // namespace SBody
