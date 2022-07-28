/**
 * @file View.cpp
 * @author Tai Zhou
 * @brief This is the file setting the view of a distant observer.
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "View.h"

#include <cmath>
#include <iostream>

#include <fmt/core.h>
#include <gsl/gsl_math.h>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "IO.h"
#include "Metric.h"
#include "Object.h"
#include "Unit.h"
#include "Utility.h"

using namespace std;

namespace SBody {
	view::view(double r, double theta, string filename) : r(r), theta(theta), sinto(sin(theta)), costo(cos(theta)), tFinal(-r - 2e4 * 1.), output(make_unique<NumPy>(filename, vector<int>({28}))) {}
	int view::traceStar(Object::star &s, int rayNO) { // FIXME:!!!!
		const double rs = s.pos[1], thetas = s.pos[2], phis = s.pos[3], sints = sin(thetas), costs = cos(thetas), sinps = sin(phis), cosps = cos(phis);
		double alpha0 = GSL_POSINF, alpha1 = (rs + 1.) * sints * sin(phis), beta0 = GSL_POSINF, beta1 = (rs + 1.) * (costs * sinto - sints * cosps * costo), ph[10], cosph, last[10], h;
#ifdef RECORD_TRACE
		vector<double> qdq(12);
		IO::NumPy rec("Trace " + to_string(rayNO), 12);
#endif
#ifdef VIEW_TAU
		integrator integ(Metric::functionTau, Metric::jacobian, 2);
#elif VIEW_HAMILTONIAN
		integrator integ(Metric::functionHamiltonian, Metric::jacobian, 2);
#else
		integrator integ(Metric::function, Metric::jacobian, 2);
#endif
		while (gsl_hypot(alpha1 - alpha0, beta1 - beta0) > epsilon * (1. + gsl_hypot(alpha1, beta1))) {
			alpha0 = alpha1;
			beta0 = beta1;
			ph[0] = 0.;
			ph[1] = r;
			ph[4] = 1.;
			ph[5] = 1.;
			ph[8] = 0.;
			ph[9] = 0.;
			if (sinto < epsilon) {
				const double k = gsl_hypot(alpha1, beta1);
				if (theta < M_PI_2) {
					ph[2] = 1e-15;
					ph[3] = alpha1 >= 0. ? acos(-beta1 / k) : M_2PI - acos(-beta1 / k);
					ph[6] = -k / gsl_pow_2(r);
				} else {
					ph[2] = M_PI - 1e-15;
					ph[3] = alpha1 >= 0. ? -acos(-beta1 / k) : M_2PI + acos(-beta1 / k);
					ph[6] = k / gsl_pow_2(r);
				}
				ph[7] = 0.;
			} else {
				ph[2] = theta;
				ph[3] = 0.,
				ph[6] = beta1 / gsl_pow_2(r);
				ph[7] = -alpha1 / (gsl_pow_2(r) * sinto);
			}
			Metric::lightNormalization(ph, 1.);
#ifdef VIEW_HAMILTONIAN
			Metric::qdq2qp(ph);
#endif
			int status = 0, fixed = 0;
			h = -1.;
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
				cosph = (rs * sints * cosps - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * cos(ph[3])) * sinto + (rs * costs - ph[1] * GSL_SIGN(ph[2]) * cos(ph[2])) * costo;
				if (cosph > rs * relAcc) {
					copy(last, last + 10, ph);
					h *= 0.3;
					fixed = 1;
					integ.reset();
				} else {
					if (ph[9] < tFinal * 1e-8) {
						ph[8] += ph[9];
						ph[9] = 0;
					}
#ifdef RECORD_TRACE
					copy(ph, ph + 8, qdq.begin());
#ifdef VIEW_HAMILTONIAN
					Metric::qp2qdq(qdq.data());
#endif
					qdq[8] = ph[8] + ph[9];
					qdq[9] = Metric::energy(qdq.data());
					qdq[10] = Metric::angularMomentum(qdq.data());
					qdq[11] = Metric::carter(qdq.data(), 0.);
					rec.save(qdq);
#endif
					if (cosph >= 0) {
						alpha1 += rs * sints * sinps - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * sin(ph[3]); // TODO:OPT!
						beta1 += (rs * costs - ph[1] * GSL_SIGN(ph[2]) * cos(ph[2])) * sinto - (rs * sints * cosps - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * cos(ph[3])) * costo;
#ifdef VIEW_HAMILTONIAN
						Metric::qp2qdq(ph);
#endif
						break;
					}
				}
			}
			if (status > 0) {
				fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
				return status;
			}
		}
		double diffAlphaBeta[4][2] = {{1.e-6, 0.}, {-1.e-6, 0.}, {0., 1.e-6}, {0., -1.e-6}};
		double res[4][8];
		for (int i = 0; i < 4; ++i) {
			alpha0 = alpha1 + diffAlphaBeta[i][0];
			beta0 = beta1 + diffAlphaBeta[i][1];
			ph[0] = 0.;
			ph[1] = r;
			ph[4] = 1.;
			ph[5] = 1.;
			ph[8] = 0.;
			ph[9] = 0.;
			if (sinto < epsilon) {
				const double k = gsl_hypot(alpha0, beta0);
				if (theta < M_PI_2) {
					ph[2] = 1e-15;
					ph[3] = alpha0 >= 0. ? acos(-beta0 / k) : M_2PI - acos(-beta0 / k);
					ph[6] = -k / gsl_pow_2(r);
				} else {
					ph[2] = M_PI - 1e-15;
					ph[3] = alpha0 >= 0. ? -acos(-beta0 / k) : M_2PI + acos(-beta0 / k);
					ph[6] = k / gsl_pow_2(r);
				}
				ph[7] = 0.;
			} else {
				ph[2] = theta;
				ph[3] = 0.,
				ph[6] = beta0 / gsl_pow_2(r);
				ph[7] = -alpha0 / (gsl_pow_2(r) * sinto);
			}
			Metric::lightNormalization(ph, 1.);
#ifdef VIEW_HAMILTONIAN
			Metric::qdq2qp(ph);
#endif
			int status = 0, fixed = 0;
			h = -1.;
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
				cosph = (rs * sints * cosps - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * cos(ph[3])) * sinto + (rs * costs - ph[1] * GSL_SIGN(ph[2]) * cos(ph[2])) * costo;
				if (cosph > rs * relAcc) {
					copy(last, last + 10, ph);
					h *= 0.3;
					fixed = 1;
					integ.reset();
				} else {
					if (ph[9] < tFinal * 1e-8) {
						ph[8] += ph[9];
						ph[9] = 0;
					}
#ifdef RECORD_TRACE
					copy(ph, ph + 8, qdq.begin());
#ifdef VIEW_HAMILTONIAN
					Metric::qp2qdq(qdq.data());
#endif
					qdq[8] = ph[8] + ph[9];
					qdq[9] = Metric::energy(qdq.data());
					qdq[10] = Metric::angularMomentum(qdq.data());
					qdq[11] = Metric::carter(qdq.data(), 0.);
					rec.save(qdq);
#endif
					if (cosph >= 0) {
#ifdef VIEW_HAMILTONIAN
						Metric::qp2qdq(ph);
#endif
						copy(ph, ph + 8, res[i]);
						break;
					}
				}
			}
			if (status > 0) {
				fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
				return status;
			}
		}
#ifdef VIEW_TAU
		output->save({alpha1, beta1, s.frequencyTau(ph), (ph[8] + ph[9]) / Unit::s});
#else
		output->save({alpha1, beta1, s.frequency(ph), (ph[8] + ph[9]) / Unit::s, res[0][1], res[0][2], res[0][3], res[0][5], res[0][6], res[0][7], res[1][1], res[1][2], res[1][3], res[1][5], res[1][6], res[1][7], res[2][1], res[2][2], res[2][3], res[2][5], res[2][6], res[2][7], res[3][1], res[3][2], res[3][3], res[3][5], res[3][6], res[3][7]});
#endif
		return 0;
	}
	int view::shadow(int n) {
		const double interval = M_2PI / n;
		NumPy rec("shadow " + to_string(Metric::a) + "," + to_string(Metric::l), {2});
		double h, rin = 2., rout = 10., rmid = 6., ph[10];
#ifdef VIEW_TAU
		integrator integ(Metric::functionTau, Metric::jacobian, 2);
#else
		integrator integ(Metric::function, Metric::jacobian, 2);
#endif
		indicators::BlockProgressBar shadowProgressBar{
			indicators::option::ShowElapsedTime{true},
			indicators::option::ShowRemainingTime{true},
			indicators::option::PrefixText{"? Shadow"},
			indicators::option::ForegroundColor{indicators::Color(4)},
			indicators::option::MaxProgress{n},
			indicators::option::FontStyles{vector<indicators::FontStyle>{indicators::FontStyle::bold}}};
		int progressBarIndex = ProgressBar::bars_.push_back(shadowProgressBar);
		for (int i = 0; i < n; ++i) {
			const double angle = i * interval, sina = sin(angle), cosa = cos(angle);
			int status = 0;
			while (rout - rin > epsilon * (rin + rout)) {
				rmid = 0.5 * (rin + rout);
				ph[0] = 0.;
				ph[1] = r;
				ph[4] = 1.;
				ph[5] = 1.;
				ph[8] = 0.;
				ph[9] = 0.;
				if (sinto < epsilon) {
					if (theta < M_PI_2) {
						ph[2] = 1e-15;
						ph[3] = M_PI_2 + angle;
						ph[6] = -rmid / gsl_pow_2(r);
					} else {
						ph[2] = M_PI - 1e-15;
						ph[3] = M_PI_2 - angle;
						ph[6] = rmid / gsl_pow_2(r);
					}
					ph[7] = 0.;
				} else {
					ph[2] = theta;
					ph[3] = 0.,
					ph[6] = rmid * sina / gsl_pow_2(r);
					ph[7] = -rmid * cosa / (gsl_pow_2(r) * sinto);
				}
				h = -1.;
				Metric::lightNormalization(ph, 1.);
				while (status <= 0 && ph[8] + ph[9] > tFinal) {
					status = integ.apply(ph + 9, tFinal, &h, ph);
					if (ph[9] < tFinal * 1e-8) {
						ph[8] += ph[9];
						ph[9] = 0.;
					}
					if (ph[4] >= 1e6 || ph[5] <= 0 || ph[1] <= 0)
						break;
				}
				if (status > 0) {
					fmt::print(stderr, "[!] view::shadow status = {}\n", status);
					return status;
				}
				if (ph[5] <= 0)
					rout = rmid;
				else
					rin = rmid;
				integ.reset();
			}
			rin -= 2. * interval * rmid;
			rout += 2. * interval * rmid;
			rec.save({rmid * cosa, rmid * sina});
			if (ProgressBar::display_)
				ProgressBar::bars_[progressBarIndex].tick();
		}
		if (ProgressBar::display_)
			ProgressBar::SetComplete(progressBarIndex, "! Shadow");
		return 0;
	}
	camera::camera(size_t pixel, double halfAngle, double r, double theta, string fileName) : view(r, theta, fileName), pixel(pixel), halfAngle(halfAngle) {
		screen = vector<vector<double>>(pixel, vector<double>(pixel));
		initials = vector<array<double, 10>>(pixel * pixel);
		const double sint = sin(theta), tana_pix = 2. * tan(halfAngle) / (r * pixel), t1 = r + 100.;
		if (theta < epsilon || M_PI - theta < epsilon) {
			for (int i = 0; i < pixel; ++i)
				for (int j = 0; j < pixel; ++j) {
					const double k = gsl_hypot(i - 0.5 * pixel + 0.5, j - 0.5 * pixel + 0.5);
					initials[i * pixel + j] = {0., r, theta < M_PI_2 ? 1e-15 : M_PI - 1e-15, 0., 1., 1., GSL_SIGN(theta - M_PI_2) * tana_pix * k, 0., 0.};
					if (k >= epsilon) {
						if (2 * i <= pixel)
							initials[i * pixel + j][3] = (theta < M_PI_2 ? 1. : -1.) * acos((j - 0.5 * pixel + 0.5) / k);
						else
							initials[i * pixel + j][3] = (theta < M_PI_2 ? 1. : -1.) * (2. * M_PI - acos((j - 0.5 * pixel + 0.5) / k));
					}
				}
		} else
			for (int i = 0; i < pixel; ++i)
				for (int j = 0; j < pixel; ++j)
					initials[i * pixel + j] = {0., r, theta, 0., 1., 1., -tana_pix * (i - 0.5 * pixel + 0.5), -tana_pix * (j - 0.5 * pixel + 0.5) / sint, 0.};
#pragma omp parallel for
		for (int p = pixel * pixel - 1; p >= 0; --p) {
			int status = 0;
			initials[p][9] = -1.;
#ifdef VIEW_TAU
			integrator integ(Metric::functionTau, Metric::jacobian, 2);
#else
			integrator integ(Metric::functionHamiltonian, Metric::jacobian, 2);
#endif
			Metric::lightNormalization(initials[p].data(), 1.);
			Metric::qdq2qp(initials[p].data());
			while (status <= 0 && initials[p][8] > t1 && initials[p][1] > 100)
				status = integ.apply(initials[p].data() + 8, t1, initials[p].data() + 9, initials[p].data());
			if (status > 0)
				fmt::print(stderr, "[!] camera::initialize status = {}\n", status);
		}
	}
	int camera::traceStar() {
		const double t1 = -1000.;
#pragma omp parallel for
		for (int p = pixel * pixel - 1; p >= 0; --p) {
			int i = p / pixel;
			int j = p - i * pixel;
			double ph[10], last[10]; // TODO: should add t_0 (initials[8]) here
			int status = 0;
			copy(initials[p].begin(), initials[p].end(), ph);
#ifdef VIEW_TAU
			integrator integ(Metric::functionTau, Metric::jacobian, 2);
#else
			integrator integ(Metric::function, Metric::jacobian, 2);
#endif
			while (status <= 0 && ph[8] > t1) {
				copy(ph, ph + 10, last);
				status = integ.apply(ph + 8, t1, ph + 9, ph);
				for (auto objP : Object::objectList)
					if (objP->hit(ph, last))
						screen[i][j] = objP->frequency(ph); // FIXME: if multi objects
				if (screen[i][j] > epsilon)
					break;
			}
			if (status > 0)
				fmt::print(stderr, "[!] camera::traceStar status = {}\n", status);
		}
		return 0;
	}
	int camera::lens() {
		const double t1 = -1000. * r, pixelPerAngle = 0.5 * pixel / halfAngle;
		NumPy rec("lens", {2});
		indicators::BlockProgressBar lensProgressBar{
			indicators::option::ShowElapsedTime{true},
			indicators::option::ShowRemainingTime{true},
			indicators::option::ForegroundColor{indicators::Color(5)},
			indicators::option::PrefixText{"? lens"},
			indicators::option::MaxProgress{pixel * pixel},
			indicators::option::FontStyles{vector<indicators::FontStyle>{indicators::FontStyle::bold}}};
		int progressBarIndex = ProgressBar::bars_.push_back(lensProgressBar);
		for (int i = 0; i < pixel; ++i)
			for (int j = 0; j < pixel; ++j) {
				double ph[10], phc[10];
				int status = 0;
				copy(initials[i * pixel + j].begin(), initials[i * pixel + j].end(), ph);
#ifdef VIEW_TAU
				integrator integ(Metric::functionTau, Metric::jacobian, 2);
#else
				integrator integ(Metric::functionHamiltonian, Metric::jacobian, 2);
#endif
				while (status <= 0 && ph[8] > t1 && ph[1] > 3. && ph[1] < 3.e2)
					status = integ.apply(ph + 8, t1, ph + 9, ph);
				if (status > 0)
					fmt::print(stderr, "[!] camera::lens status = {}\n", status);
				if (ph[1] <= 3.)
					rec.save({NAN, NAN});
				else {
					Metric::qp2qdq(ph);
					if (ph[2] < 0)
						ph[2] += M_PI;
					Metric::s2c(ph, phc);
					rec.save({phc[6] * pixelPerAngle, (phc[7] * sinto - phc[5] * costo) * pixelPerAngle});
				}
				if (ProgressBar::display_)
					ProgressBar::bars_[progressBarIndex].tick();
			}
		if (ProgressBar::display_)
			ProgressBar::SetComplete(progressBarIndex, "! lens");
		return 0;
	}
	int camera::save() {
		for (const vector<double> &line : screen)
			output->save(line);
		return 0;
	}
} // namespace SBody
