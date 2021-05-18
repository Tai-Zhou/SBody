#include "View.h"

#include <cmath>
#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>

#include "Constant.h"
#include "IO.h"
#include "Metric.h"

using namespace std;

namespace SBody {
	view::view(size_t pixel, double viewAngle, double r, double theta, double phi, double tFinal, size_t duration, size_t frame)
		: pixel(pixel),
		  viewAngle(viewAngle),
		  r(r),
		  theta(theta),
		  phi(phi),
		  tFinal(tFinal),
		  duration(duration),
		  frame(frame) {
		initials = vector<array<double, 9>>(pixel * pixel);
	}
	void view::traceBack(size_t NSK, source *params) {
		const double mass = params->mass, sint = sin(theta), tana_pix = 2. * tan(0.5 * viewAngle) / (r * pixel), tFrame = tFinal / duration;
		double ph[8] = {0., r, theta, phi, 1., -1., 0., 0.}, last[3], t = 0, tStep = 0, h = 1e-3;
		int status = 0;
		integrator integ(params);
		vector<vector<int>> screen(pixel, vector<int>(pixel));
		if (theta < Constant::epsilon || M_PI - theta < Constant::epsilon) {
			for (int i = 0; i < pixel; ++i)
				for (int j = 0; j < pixel; ++j) {
					const double k = gsl_hypot(i - 0.5 * pixel + 0.5, j - 0.5 * pixel + 0.5);
					ph[0] = 0;
					ph[1] = r;
					ph[2] = theta;
					ph[4] = 1;
					ph[5] = -1;
					ph[6] = (theta < M_PI_2 ? 1 : -1) * tana_pix * k;
					ph[7] = 0;
					if (k < Constant::epsilon)
						ph[3] = 0;
					else {
						if (2 * i <= pixel)
							ph[3] = (theta < M_PI_2 ? 1. : -1.) * acos((j - 0.5 * pixel + 0.5) / k);
						else
							ph[3] = (theta < M_PI_2 ? 1. : -1.) * (2. * M_PI - acos((j - 0.5 * pixel + 0.5) / k));
					}
					if (NSK == 1)
						Metric::Schwarzschild::lightNormalization(ph, params);
					else if (NSK == 2)
						Metric::Kerr::lightNormalization(ph, params);
					else if (NSK == 3)
						Metric::KerrH::lightNormalization(ph, params);
					while (status == 0) {
						status = integ.apply(&t, GSL_POSINF, &h, ph);
						if (ph[1] < 100 * mass) {
							initials[i * pixel + j] = {ph[0], ph[1], ph[2], ph[3], ph[4], ph[5], ph[6], ph[7], t};
							break;
						}
					}
				}
		}
		else
			for (int i = 0; i < pixel; ++i)
				for (int j = 0; j < pixel; ++j) {
					ph[0] = 0;
					ph[1] = r;
					ph[2] = theta;
					ph[3] = phi;
					ph[4] = 1;
					ph[5] = -1;
					ph[6] = tana_pix * (i - 0.5 * pixel + 0.5);
					ph[7] = tana_pix * (j - 0.5 * pixel + 0.5) / sint;
					if (NSK == 1)
						Metric::Schwarzschild::lightNormalization(ph, params);
					else if (NSK == 2)
						Metric::Kerr::lightNormalization(ph, params);
					else if (NSK == 3)
						Metric::KerrH::lightNormalization(ph, params);
					while (status == 0) {
						status = integ.apply(&t, GSL_POSINF, &h, ph);
						if (ph[1] < 100 * mass) {
							initials[i * pixel + j] = {ph[0], ph[1], ph[2], ph[3], ph[4], ph[5], ph[6], ph[7], t};
							break;
						}
					}
				}
		for (int i = 0; i < pixel; ++i)
			for (int j = 0; j < pixel; ++j) {
				t = 0;
				tStep = 0;
				h = 1e-3;
				while (tStep < tFinal) {
					tStep += tFrame;
					while (status == 0 && t < tStep) {
						last[0] = ph[1];
						last[1] = ph[2];
						last[2] = ph[3];
						status = integ.apply(&t, tStep, &h, ph);
						if (ph[1] > 6 * mass && ph[1] < 12 * mass && (ph[2] - M_PI_2) * (last[1] - M_PI_2) < 0)
							screen[i][j] = 1;
					}
				}
			}
		IO::NumPySave(screen, "screen");
	}
} // namespace SBody
