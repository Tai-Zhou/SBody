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
	view::view(double r, double theta, double tFinal, size_t duration, size_t frame) : r(r), theta(theta), sinto(sin(theta)), costo(cos(theta)), tFinal(r + 2e4 * Metric::m), duration(duration), frame(frame) {}
	void view::traceBack(Object::star &s, int rayNO) {
		//Metric::setMetric(2, 4e6, 1);
		const double rs = s.pos[1], thetas = s.pos[2], phis = s.pos[3], sints = sin(thetas), costs = cos(thetas), sinps = sin(phis), cosps = cos(phis);
		double alpha0 = GSL_POSINF, alpha1 = rs * sints * sin(phis), beta0 = GSL_POSINF, beta1 = rs * costs * sinto - rs * sints * cosps * costo, ph[9], cosph, last[9], h;
		vector<double> qdq(12);
		vector<vector<double>> rec;
		while (abs(alpha1 - alpha0) > epsilon * (1. + gsl_hypot(alpha1, beta1)) || abs(beta1 - beta0) > epsilon * (1. + gsl_hypot(alpha1, beta1))) {
			rec.clear();
			alpha0 = alpha1;
			beta0 = beta1;
			ph[0] = 0.;
			ph[1] = r;
			ph[4] = 1.;
			ph[5] = -1.;
			ph[8] = 0.;
			if (sinto < epsilon) {
				const double k = gsl_hypot(alpha1, beta1);
				if (M_PI_2 > theta) {
					ph[2] = 1e-15;
					ph[3] = alpha1 >= 0. ? acos(-beta1 / k) : M_2PI - acos(-beta1 / k);
					ph[6] = k / gsl_pow_2(r);
				}
				else {
					ph[2] = M_PI - 1e-15;
					ph[3] = alpha1 >= 0. ? -acos(-beta1 / k) : M_2PI + acos(-beta1 / k);
					ph[6] = -k / gsl_pow_2(r);
				}
				ph[7] = 0.;
			}
			else {
				ph[2] = theta;
				ph[3] = 0.,
				ph[6] = -beta1 / gsl_pow_2(r);
				ph[7] = alpha1 / (gsl_pow_2(r) * sinto);
			}
			h = 1e-3;
			Metric::lightNormalization(ph, 1.);
			const int kerrh = 1;
			if (kerrh)
				Metric::KerrH::qdq2qp(ph);
			integrator integ(kerrh ? Metric::KerrH::function : Metric::function, kerrh ? Metric::KerrH::jacobian : Metric::jacobian, 0);
			int status = 0, fixed = 0;
			while (status == 0 && ph[8] < tFinal) {
				copy(ph, ph + 9, last);
				if (fixed)
					status = integ.apply_fixed(ph + 8, h, ph);
				else
					status = integ.apply(ph + 8, tFinal, &h, ph);
				cosph = (rs * sints * cosps - ph[1] * sin(ph[2]) * cos(ph[3])) * sinto + (rs * costs - ph[1] * cos(ph[2])) * costo;
				if (cosph >= 0) {
					if (cosph > rs * epsilon) {
						copy(last, last + 9, ph);
						h *= 0.5;
						integ.reset();
						fixed = 1;
					}
					else {
						alpha1 += rs * sints * sinps - ph[1] * sin(ph[2]) * sin(ph[3]); //TODO:OPT!
						beta1 += rs * costs * sinto - rs * sints * cosps * costo - ph[1] * cos(ph[2]) * sinto + ph[1] * sin(ph[2]) * cos(ph[3]) * costo;
						break;
					}
				}
				else {
					copy(ph, ph + 8, qdq.begin());
					if (kerrh)
						Metric::KerrH::qp2qdq(qdq.data());
					qdq[8] = ph[8] / Constant::s;
					qdq[9] = Metric::energy(qdq.data());
					qdq[10] = Metric::angularMomentum(qdq.data());
					qdq[11] = Metric::carterL(qdq.data());
					rec.push_back(qdq);
				}
			}
		}
		screen.push_back({alpha1, beta1, s.frequency(qdq.data()), ph[8] / Constant::s});
		IO::NumPy<double> output("Trace " + to_string(rayNO));
		output.save(rec);
		//Metric::setMetric(1, 4e6, 1);
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
