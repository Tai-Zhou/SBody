#include "View.h"

#include <cmath>
#include <iostream>

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
	view::view(double r, double theta) : r(r), theta(theta), sinto(sin(theta)), costo(cos(theta)), tFinal(r + 2e4 * Metric::m) {}
	int view::traceBack(Object::star &s, int rayNO) { //FIXME:!!!!
		const double rs = s.pos[1], thetas = s.pos[2], phis = s.pos[3], sints = sin(thetas), costs = cos(thetas), sinps = sin(phis), cosps = cos(phis);
		double alpha0 = GSL_POSINF, alpha1 = (rs + Metric::m) * sints * sin(phis), beta0 = GSL_POSINF, beta1 = (rs + Metric::m) * (costs * sinto - sints * cosps * costo), ph[10], cosph, last[10], h;
		//Metric::a *= -1.;
		vector<double> qdq(12);
		vector<vector<double>> rec;
		const int Hamiltonian = 0, record = 0;
		while (gsl_hypot(alpha1 - alpha0, beta1 - beta0) > epsilon * (1. + gsl_hypot(alpha1, beta1))) {
			rec.clear();
			alpha0 = alpha1;
			beta0 = beta1;
			ph[0] = 0.;
			ph[1] = r;
			ph[4] = 1.;
			ph[5] = -1.;
			ph[8] = 0.;
			ph[9] = 0.;
			if (sinto < epsilon) {
				const double k = gsl_hypot(alpha1, beta1);
				if (theta < M_PI_2) {
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
			if (Hamiltonian)
				Metric::qdq2qp(ph);
			integrator integ(Metric::function, Metric::jacobian, 2);
			int status = 0, fixed = 0;
			while (status <= 0) {
				if (status == -1) {
					h *= 0.5;
					fixed = 0;
				}
				copy(ph, ph + 10, last);
				if (fixed)
					status = integ.apply_fixed(ph + 9, h, ph);
				else
					status = integ.apply(ph + 9, tFinal, &h, ph);
				cosph = (rs * sints * cosps - ph[1] * sign(ph[2]) * sin(ph[2]) * cos(ph[3])) * sinto + (rs * costs - ph[1] * sign(ph[2]) * cos(ph[2])) * costo;
				if (cosph > rs * relAcc) {
					copy(last, last + 10, ph);
					h *= 0.5;
					fixed = 1;
					integ.reset();
				}
				else {
					if (ph[9] > tFinal * 1e-8) {
						ph[8] += ph[9];
						ph[9] = 0;
					}
					copy(ph, ph + 8, qdq.begin());
					if (Hamiltonian)
						Metric::qp2qdq(qdq.data());
					qdq[8] = (ph[8] + ph[9]) / Constant::s;
					qdq[9] = Metric::energy(qdq.data());
					qdq[10] = Metric::angularMomentum(qdq.data());
					qdq[11] = Metric::carter(qdq.data(), 0.);
					if (record)
						rec.push_back(qdq);
					if (cosph >= 0) {
						alpha1 += rs * sints * sinps - ph[1] * sign(ph[2]) * sin(ph[2]) * sin(ph[3]); //TODO:OPT!
						beta1 += (rs * costs - ph[1] * sign(ph[2]) * cos(ph[2])) * sinto - (rs * sints * cosps - ph[1] * sign(ph[2]) * sin(ph[2]) * cos(ph[3])) * costo;
						break;
					}
				}
			}
			if (status > 0) {
				cerr << "[!] view::traceBack status = " << status << endl;
				return status;
			}
		}
		//Metric::a *= -1.;
		screen.push_back({alpha1, beta1, s.frequency(qdq.data()), qdq[8]});
		if (record) {
			IO::NumPy<double> output("Trace " + to_string(rayNO));
			output.save(rec);
		}
		return 0;
	}
	int view::shadow(int n) {
		const double interval = M_2PI / n;
		vector<vector<double>> rec(n);
		Metric::a *= -1.;
		Metric::l *= -1.;
		double h, rin = 2. * Metric::m, rout = 10. * Metric::m, rmid = 6. * Metric::m, ph[10];
		for (int i = 0; i < n; ++i) {
			IO::progressBar.set_progress(100. * i / n);
			const double angle = i * interval, sina = sin(angle), cosa = cos(angle);
			integrator integ(Metric::KerrTaubNUT::functionTau, Metric::jacobian, 2);
			int status = 0;
			rin -= 2. * interval * rmid;
			rout += 1.5 * interval * rmid;
			while (rout - rin > epsilon * (rin + rout)) {
				rmid = 0.5 * (rin + rout);
				ph[0] = 0.;
				ph[1] = r;
				ph[4] = 1.;
				ph[5] = -1.;
				ph[8] = 0.;
				ph[9] = 0.;
				if (sinto < epsilon) {
					if (theta < M_PI_2) {
						ph[2] = 1e-15;
						ph[3] = M_PI_2 + angle;
						ph[6] = rmid / gsl_pow_2(r);
					}
					else {
						ph[2] = M_PI - 1e-15;
						ph[3] = M_PI_2 - angle;
						ph[6] = -rmid / gsl_pow_2(r);
					}
					ph[7] = 0.;
				}
				else {
					ph[2] = theta;
					ph[3] = 0.,
					ph[6] = -rmid * sina / gsl_pow_2(r);
					ph[7] = rmid * cosa / (gsl_pow_2(r) * sinto);
				}
				h = 1.;
				Metric::lightNormalization(ph, 1.);
				while (status <= 0 && ph[8] + ph[9] < tFinal) {
					status = integ.apply(ph + 9, tFinal, &h, ph);
					if (ph[9] > tFinal * 1e-8) {
						ph[8] += ph[9];
						ph[9] = 0;
					}
					if (ph[4] >= 1e6 || ph[5] >= epsilon)
						break;
				}
				if (status > 0) {
					cerr << "[!] view::shadow status = " << status << endl;
					return status;
				}
				if (ph[5] >= epsilon)
					rout = rmid;
				else
					rin = rmid;
				integ.reset();
			}
			rec[i] = {rmid * cosa, rmid * sina};
		}
		Metric::a *= -1.;
		Metric::l *= -1.;
		IO::NumPy<double> output("shadow " + to_string(Metric::a / Metric::m) + "," + to_string(Metric::l / Metric::m));
		IO::progressBar.set_progress(100.);
		IO::progressBar.set_option(indicators::option::PrefixText{"Complete"});
		IO::progressBar.mark_as_completed();
		indicators::show_console_cursor(true);
		return output.save(rec);
	}
	int view::save(string fileName) {
		IO::NumPy<double> output(fileName);
		return output.save(screen);
	}
	camera::camera(size_t pixel, double viewAngle, double r, double theta) : view(r, theta), pixel(pixel), viewAngle(viewAngle) {
		screen = vector<vector<double>>(pixel, vector<double>(pixel));
		initials = vector<array<double, 9>>(pixel * pixel);
	}
	int camera::initialize() {
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
				if (status > 0)
					cerr << "[!] camera::initialize status = " << status << endl;
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
				if (status > 0)
					cerr << "[!] camera::initialize status = " << status << endl;
			}
		return 0;
	}
	int camera::traceBack() {
		const double t1 = 1000. * Metric::m;
#pragma omp parallel for
		for (int p = pixel * pixel - 1; p >= 0; --p) {
			int i = p / pixel;
			int j = p - i * pixel;
			double ph[8], last[8], t = 0., h = 1e-3; // TODO: should add t_0 (initials[8]) here
			int status = 0;
			integrator integ(Metric::function, Metric::jacobian, 0);
			copy(initials[p].data(), initials[p].data() + 8, ph);
			while (status == 0 && t < t1) {
				copy(ph, ph + 3, last);
				status = integ.apply(&t, t1, &h, ph);
				for (auto objP : Object::objectList)
					if (objP->hit(ph, last))
						screen[i][j] = objP->frequency(ph); //FIXME: if multi objects
				if (screen[i][j] > epsilon)
					break;
			}
			if (status > 0)
				cerr << "[!] camera::traceBack status = " << status << endl;
		}
		return 0;
	}
} // namespace SBody
