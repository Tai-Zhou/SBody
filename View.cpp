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
	void view::traceBack() {
		const double sint = sin(theta), tana_pix = 2. * tan(0.5 * viewAngle) / (r * pixel), tFrame = 30 * Constant::R_sun;
		double ph[8] = {0., r, theta, phi, 1., -1., 0., 0.}, last[8], t = 0, tStep = 0, h = 1e-3;
		int status = 0;
		integrator integ(Metric::function, Metric::jacobian, 0);
		vector<vector<int>> screen(pixel, vector<int>(pixel));
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
					while (status == 0 && t < tFinal) {
						status = integ.apply(&t, tFinal, &h, ph);
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
					ph[3] = phi;
					ph[4] = 1;
					ph[5] = -1;
					ph[6] = tana_pix * (i - 0.5 * pixel + 0.5);
					ph[7] = tana_pix * (j - 0.5 * pixel + 0.5) / sint;
					Metric::lightNormalization(ph);
					while (status == 0 && t < tFinal) {
						status = integ.apply(&t, tFinal, &h, ph);
						if (ph[1] < 100 * Metric::m) {
							initials[i * pixel + j] = {ph[0], ph[1], ph[2], ph[3], ph[4], ph[5], ph[6], ph[7], t};
							break;
						}
					}
				}
		double pos[8] = {0., 6. * Metric::m, M_PI_2, M_PI};
		Object::star fStar(5 * Constant::R_sun, pos, 1);
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
				t = initials[i * pixel + j][8];
				while (status == 0 && t < tFinal) {
					last[1] = ph[1];
					last[2] = ph[2];
					last[3] = ph[3];
					status = integ.apply(&t, tFinal, &h, ph);
					if (fStar.hit(ph, last)) {
						screen[i][j] = 1;
						break;
					}
				}
			}
		IO::NumPySave(screen, "screen");
	}
} // namespace SBody
