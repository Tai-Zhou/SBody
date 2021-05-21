#include "View.h"

#include <cmath>
#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>

#include "Constant.h"
#include "IO.h"
#include "Metric.h"
#include "Object.h"
#include "Utility.h"

using namespace std;

namespace SBody {
	view::view(double viewAngle, double r, double theta, double tFinal, size_t duration, size_t frame)
		: viewAngle(viewAngle), r(r), theta(theta), tFinal(tFinal), duration(duration), frame(frame) {}
	void view::traceBack() {}
	camera::camera(size_t pixel, double viewAngle, double r, double theta, double tFinal, size_t duration, size_t frame) : view(viewAngle, r, theta, tFinal, duration, frame), pixel(pixel) {
		initials = vector<array<double, 9>>(pixel * pixel);
		screen = vector<vector<int>>(pixel, vector<int>(pixel));
		const double sint = sin(theta), tana_pix = 2. * tan(0.5 * viewAngle) / (r * pixel);
		double ph[8] = {0., r, theta, 0, 1., -1., 0., 0.}, t = 0, h = 1e-3;
		int status = 0;
		integrator integ(Metric::function, Metric::jacobian, 0);
		if (theta < epsilon || M_PI - theta < epsilon) {
			for (int i = 0; i < pixel; ++i)
				for (int j = 0; j < pixel; ++j) {
					integ.reset();
					t = 0;
					h = 1e-3;
					const double k = gsl_hypot(i - 0.5 * pixel + 0.5, j - 0.5 * pixel + 0.5);
					ph[0] = 0;
					ph[1] = r;
					ph[2] = theta < M_PI_2 ? 0 : M_PI;
					ph[4] = 1;
					ph[5] = -1;
					ph[6] = sign(M_PI_2 - theta) * tana_pix * k;
					ph[7] = 0;
					if (k < epsilon)
						ph[3] = 0;
					else {
						if (2 * i <= pixel)
							ph[3] = (theta < M_PI_2 ? 1. : -1.) * acos((j - 0.5 * pixel + 0.5) / k);
						else
							ph[3] = (theta < M_PI_2 ? 1. : -1.) * (2. * M_PI - acos((j - 0.5 * pixel + 0.5) / k));
					}
					Metric::lightNormalization(ph);
					while (status == 0 && t < r + 100) {
						status = integ.apply(&t, r + 100, &h, ph);
						if (ph[1] < 100 * Metric::m) {
							initials[i * pixel + j] = {ph[0], ph[1], ph[2], ph[3], ph[4], ph[5], ph[6], ph[7], t};
							break;
						}
					}
				}
		}
		else
			for (int i = 0; i < pixel; ++i)
				for (int j = 0; j < pixel; ++j) {
					integ.reset();
					t = 0;
					h = 1e-3;
					ph[0] = 0;
					ph[1] = r;
					ph[2] = theta;
					ph[3] = 0;
					ph[4] = 1;
					ph[5] = -1;
					ph[6] = tana_pix * (i - 0.5 * pixel + 0.5);
					ph[7] = tana_pix * (j - 0.5 * pixel + 0.5) / sint;
					Metric::lightNormalization(ph);
					while (status == 0 && t < r + 100. * Metric::m) {
						status = integ.apply(&t, r + 100. * Metric::m, &h, ph);
						if (ph[1] < 100 * Metric::m) {
							initials[i * pixel + j] = {ph[0], ph[1], ph[2], ph[3], ph[4], ph[5], ph[6], ph[7], t};
							break;
						}
					}
				}
	}
	void camera::traceBack() {
		double ph[8] = {0., r, theta, 0, 1., -1., 0., 0.}, last[8], t = 0, h = 1e-3;
		int status = 0;
		integrator integ(Metric::function, Metric::jacobian, 0);
		for (int i = 0; i < pixel; ++i)
			for (int j = 0; j < pixel; ++j) {
				integ.reset();
				h = 1e-3;
				ph[0] = initials[i * pixel + j][0];
				ph[1] = initials[i * pixel + j][1];
				ph[2] = initials[i * pixel + j][2];
				ph[3] = initials[i * pixel + j][3];
				ph[4] = initials[i * pixel + j][4];
				ph[5] = initials[i * pixel + j][5];
				ph[6] = initials[i * pixel + j][6];
				ph[7] = initials[i * pixel + j][7];
				t = 0; // TODO: should add t_0 (initials[8]) here
				while (status == 0 && t < 1000. * Metric::m) {
					last[1] = ph[1];
					last[2] = ph[2];
					last[3] = ph[3];
					status = integ.apply(&t, 1000. * Metric::m, &h, ph);
					for (auto objP : Object::objectList)
						if (objP->hit(ph, last))
							screen[i][j] = 1;
					if (screen[i][j])
						break;
				}
			}
	}
	void camera::save() {
		IO::NumPySave(screen, "screen");
	}
} // namespace SBody
