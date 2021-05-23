#include "View.h"

#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <omp.h>

#include "Constant.h"
#include "IO.h"
#include "Metric.h"
#include "Object.h"
#include "Utility.h"

using namespace std;

namespace SBody {
	view::view(double r, double theta, double tFinal, size_t duration, size_t frame) : r(r), theta(theta), tFinal(tFinal), duration(duration), frame(frame) {}
	void view::traceBack(Object::star &s) {
		const double sints = sin(s.pos[2]), costs = cos(s.pos[2]), sinto = sin(theta), costo = cos(theta), rs = s.pos[1], phis = s.pos[3], t1 = r + 100. * Metric::m;
		double alpha0 = GSL_POSINF, alpha1 = rs * sints * sin(phis), beta0 = GSL_POSINF, beta1 = rs * costs * sinto - rs * sints * cos(phis) * costo, ph[8], rp = rs, phip = phis;
		while (abs(alpha1 - alpha0) > abs(1e-3 * alpha1) || abs(beta1 - beta0) > abs(1e-3 * beta1)) {
			alpha0 = alpha1;
			beta0 = beta1;
			ph[0] = 0.;
			ph[1] = r;
			ph[2] = theta;
			ph[4] = 1.;
			ph[5] = -1.;
			if (sinto < epsilon) {
				const double k = gsl_hypot(alpha1, beta1);
				if (alpha1 >= 0)
					ph[3] = (theta < M_PI_2 ? 1. : -1.) * acos(-beta1 / k);
				else
					ph[3] = M_2PI - (theta < M_PI_2 ? 1. : -1.) * acos(-beta1 / k);
				ph[6] = sign(M_PI_2 - theta) * k / gsl_pow_2(r);
				ph[7] = 0.;
			}
			else {
				ph[3] = 0.,
				ph[6] = -beta1 / gsl_pow_2(r);
				ph[7] = alpha1 / (gsl_pow_2(r) * sinto);
			}
			double lastt, t = 0, h = 1e-3;
			Metric::lightNormalization(ph, 1.);
			integrator integ(Metric::function, Metric::jacobian, 0);
			int status = 0;
			while (status == 0 && t < t1) {
				lastt = ph[2];
				status = integ.apply(&t, t1, &h, ph);
				if (oppositeSign(ph[2] - s.pos[2], lastt - s.pos[2])) {
					const double dt = (s.pos[2] - ph[2]) / ph[6];
					ph[1] += dt * ph[5];
					ph[2] = s.pos[2];
					ph[3] += dt * ph[7];
					const double dphisint = (ph[3] - phis) * sints;
					rp += rs - ph[1] * cos(dphisint);
					phip -= asin(ph[1] * sin(dphisint) / (rp * sints));
					alpha1 = rp * sints * sin(phip);
					beta1 = rp * costs * sinto - rp * sints * cos(phip) * costo;
					break;
				}
			}
		}
		screen.push_back({alpha1, beta1, s.frequency(ph)});
	}
	void view::save(string fileName) {
		IO::NumPy<double> output(fileName);
		output.save(screen);
	}
	camera::camera(size_t pixel, double viewAngle, double r, double theta, double tFinal, size_t duration, size_t frame) : view(r, theta, tFinal, duration, frame), pixel(pixel), viewAngle(viewAngle) {
		screen = vector<vector<double>>(pixel, vector<double>(pixel));
		initials = vector<array<double, 9>>(pixel * pixel);
	}
	void camera::initialize() {
		const double sint = sin(theta), tana_pix = 2. * tan(0.5 * viewAngle) / (r * pixel), t1 = r + 100. * Metric::m;
		if (theta < epsilon || M_PI - theta < epsilon) {
#pragma omp parallel for
			for (int p = pixel * pixel - 1; p >= 0; --p) {
				int i = p / pixel;
				int j = p - i * pixel;
				const double k = gsl_hypot(i - 0.5 * pixel + 0.5, j - 0.5 * pixel + 0.5);
				double ph[8] = {0., r, theta < M_PI_2 ? 0. : M_PI, 0., 1., -1., sign(M_PI_2 - theta) * tana_pix * k, 0.}, t = 0., h = 1e-3;
				int status = 0;
				integrator integ(Metric::function, Metric::jacobian, 0);
				if (k >= epsilon) {
					if (2 * i <= pixel)
						ph[3] = (theta < M_PI_2 ? 1. : -1.) * acos((j - 0.5 * pixel + 0.5) / k);
					else
						ph[3] = (theta < M_PI_2 ? 1. : -1.) * (2. * M_PI - acos((j - 0.5 * pixel + 0.5) / k));
				}
				Metric::lightNormalization(ph, 1.);
				while (status == 0 && t < t1) {
					status = integ.apply(&t, t1, &h, ph);
					if (ph[1] < 100 * Metric::m) {
						initials[p] = {ph[0], ph[1], ph[2], ph[3], ph[4], ph[5], ph[6], ph[7], t};
						break;
					}
				}
			}
		}
		else
#pragma omp parallel for
			for (int p = pixel * pixel - 1; p >= 0; --p) {
				int i = p / pixel;
				int j = p - i * pixel;
				double ph[8] = {0., r, theta, 0., 1., -1., tana_pix * (i - 0.5 * pixel + 0.5), tana_pix * (j - 0.5 * pixel + 0.5) / sint}, t = 0., h = 1e-3;
				int status = 0;
				integrator integ(Metric::function, Metric::jacobian, 0);
				Metric::lightNormalization(ph, 1.);
				while (status == 0 && t < t1) {
					status = integ.apply(&t, t1, &h, ph);
					if (ph[1] < 100 * Metric::m) {
						initials[p] = {ph[0], ph[1], ph[2], ph[3], ph[4], ph[5], ph[6], ph[7], t};
						break;
					}
				}
			}
	}
	void camera::traceBack() {
		const double t1 = 1000. * Metric::m;
#pragma omp parallel for
		for (int p = pixel * pixel - 1; p >= 0; --p) {
			int i = p / pixel;
			int j = p - i * pixel;
			double ph[8], last[8], t = 0., h = 1e-3; // TODO: should add t_0 (initials[8]) here
			int status = 0;
			integrator integ(Metric::function, Metric::jacobian, 0);
			ph[0] = initials[p][0];
			ph[1] = initials[p][1];
			ph[2] = initials[p][2];
			ph[3] = initials[p][3];
			ph[4] = initials[p][4];
			ph[5] = initials[p][5];
			ph[6] = initials[p][6];
			ph[7] = initials[p][7];
			while (status == 0 && t < t1) {
				last[1] = ph[1];
				last[2] = ph[2];
				last[3] = ph[3];
				status = integ.apply(&t, t1, &h, ph);
				for (auto objP : Object::objectList)
					if (objP->hit(ph, last))
						screen[i][j] = objP->frequency(ph); //FIXME: if multi objects
				if (screen[i][j] > epsilon)
					break;
			}
		}
	}
} // namespace SBody
