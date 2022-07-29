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
	View::View(double r, double theta, string file_name) : r_(r), theta_(theta), sin_theta_observer_(sin(theta)), cos_theta_observer_(cos(theta)), t_final_(-r - 2e4 * 1.), output_(make_unique<NumPy>(file_name, vector<int>({28}))) {}
	int View::TraceStar(Star &star, int ray_number) { // FIXME:!!!!
		const double rs = star.pos[1], thetas = star.pos[2], phis = star.pos[3], sints = sin(thetas), costs = cos(thetas), sinps = sin(phis), cosps = cos(phis);
		double alpha0 = GSL_POSINF, alpha1 = (rs + 1.) * sints * sin(phis), beta0 = GSL_POSINF, beta1 = (rs + 1.) * (costs * sin_theta_observer_ - sints * cosps * cos_theta_observer_), ph[10], cosph, last[10], h;
#ifdef RECORD_TRACE
		vector<double> qdq(12);
		IO::NumPy rec("Trace " + to_string(ray_number), 12);
#endif
#ifdef VIEW_TAU
		Integrator integrator(Metric::functionTau, Metric::jacobian, 2);
#elif VIEW_HAMILTONIAN
		Integrator integrator(Metric::functionHamiltonian, Metric::jacobian, 2);
#else
		Integrator integrator(Metric::function, Metric::jacobian, 2);
#endif
		while (gsl_hypot(alpha1 - alpha0, beta1 - beta0) > epsilon * (1. + gsl_hypot(alpha1, beta1))) {
			alpha0 = alpha1;
			beta0 = beta1;
			ph[0] = 0.;
			ph[1] = r_;
			ph[4] = 1.;
			ph[5] = 1.;
			ph[8] = 0.;
			ph[9] = 0.;
			if (sin_theta_observer_ < epsilon) {
				const double k = gsl_hypot(alpha1, beta1);
				if (theta_ < M_PI_2) {
					ph[2] = 1e-15;
					ph[3] = alpha1 >= 0. ? acos(-beta1 / k) : M_2PI - acos(-beta1 / k);
					ph[6] = -k / gsl_pow_2(r_);
				} else {
					ph[2] = M_PI - 1e-15;
					ph[3] = alpha1 >= 0. ? -acos(-beta1 / k) : M_2PI + acos(-beta1 / k);
					ph[6] = k / gsl_pow_2(r_);
				}
				ph[7] = 0.;
			} else {
				ph[2] = theta_;
				ph[3] = 0.,
				ph[6] = beta1 / gsl_pow_2(r_);
				ph[7] = -alpha1 / (gsl_pow_2(r_) * sin_theta_observer_);
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
					status = integrator.ApplyFixedStep(ph + 9, h, ph);
				else
					status = integrator.Apply(ph + 9, t_final_, &h, ph);
				cosph = (rs * sints * cosps - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * cos(ph[3])) * sin_theta_observer_ + (rs * costs - ph[1] * GSL_SIGN(ph[2]) * cos(ph[2])) * cos_theta_observer_;
				if (cosph > rs * relAcc) {
					copy(last, last + 10, ph);
					h *= 0.3;
					fixed = 1;
					integrator.Reset();
				} else {
					if (ph[9] < t_final_ * 1e-8) {
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
						beta1 += (rs * costs - ph[1] * GSL_SIGN(ph[2]) * cos(ph[2])) * sin_theta_observer_ - (rs * sints * cosps - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * cos(ph[3])) * cos_theta_observer_;
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
			ph[1] = r_;
			ph[4] = 1.;
			ph[5] = 1.;
			ph[8] = 0.;
			ph[9] = 0.;
			if (sin_theta_observer_ < epsilon) {
				const double k = gsl_hypot(alpha0, beta0);
				if (theta_ < M_PI_2) {
					ph[2] = 1e-15;
					ph[3] = alpha0 >= 0. ? acos(-beta0 / k) : M_2PI - acos(-beta0 / k);
					ph[6] = -k / gsl_pow_2(r_);
				} else {
					ph[2] = M_PI - 1e-15;
					ph[3] = alpha0 >= 0. ? -acos(-beta0 / k) : M_2PI + acos(-beta0 / k);
					ph[6] = k / gsl_pow_2(r_);
				}
				ph[7] = 0.;
			} else {
				ph[2] = theta_;
				ph[3] = 0.,
				ph[6] = beta0 / gsl_pow_2(r_);
				ph[7] = -alpha0 / (gsl_pow_2(r_) * sin_theta_observer_);
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
					status = integrator.ApplyFixedStep(ph + 9, h, ph);
				else
					status = integrator.Apply(ph + 9, t_final_, &h, ph);
				cosph = (rs * sints * cosps - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * cos(ph[3])) * sin_theta_observer_ + (rs * costs - ph[1] * GSL_SIGN(ph[2]) * cos(ph[2])) * cos_theta_observer_;
				if (cosph > rs * relAcc) {
					copy(last, last + 10, ph);
					h *= 0.3;
					fixed = 1;
					integrator.Reset();
				} else {
					if (ph[9] < t_final_ * 1e-8) {
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
		output->save({alpha1, beta1, star.FrequencyTau(ph), (ph[8] + ph[9]) / Unit::star});
#else
		output_->save({alpha1, beta1, star.Frequency(ph), (ph[8] + ph[9]) / Unit::s, res[0][1], res[0][2], res[0][3], res[0][5], res[0][6], res[0][7], res[1][1], res[1][2], res[1][3], res[1][5], res[1][6], res[1][7], res[2][1], res[2][2], res[2][3], res[2][5], res[2][6], res[2][7], res[3][1], res[3][2], res[3][3], res[3][5], res[3][6], res[3][7]});
#endif
		return 0;
	}
	int View::Shadow(int n) {
		const double interval = M_2PI / n;
		NumPy rec("shadow " + to_string(Metric::a) + "," + to_string(Metric::l), {2});
		double h, rin = 2., rout = 10., rmid = 6., ph[10];
#ifdef VIEW_TAU
		Integrator integrator(Metric::functionTau, Metric::jacobian, 2);
#else
		Integrator integrator(Metric::function, Metric::jacobian, 2);
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
				ph[1] = r_;
				ph[4] = 1.;
				ph[5] = 1.;
				ph[8] = 0.;
				ph[9] = 0.;
				if (sin_theta_observer_ < epsilon) {
					if (theta_ < M_PI_2) {
						ph[2] = 1e-15;
						ph[3] = M_PI_2 + angle;
						ph[6] = -rmid / gsl_pow_2(r_);
					} else {
						ph[2] = M_PI - 1e-15;
						ph[3] = M_PI_2 - angle;
						ph[6] = rmid / gsl_pow_2(r_);
					}
					ph[7] = 0.;
				} else {
					ph[2] = theta_;
					ph[3] = 0.,
					ph[6] = rmid * sina / gsl_pow_2(r_);
					ph[7] = -rmid * cosa / (gsl_pow_2(r_) * sin_theta_observer_);
				}
				h = -1.;
				Metric::lightNormalization(ph, 1.);
				while (status <= 0 && ph[8] + ph[9] > t_final_) {
					status = integrator.Apply(ph + 9, t_final_, &h, ph);
					if (ph[9] < t_final_ * 1e-8) {
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
				integrator.Reset();
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
	Camera::Camera(size_t pixel, double half_angle, double r, double theta, string file_name) : View(r, theta, file_name), pixel_(pixel), half_angle_(half_angle) {
		screen_ = vector<vector<double>>(pixel, vector<double>(pixel));
		initials_ = vector<array<double, 10>>(pixel * pixel);
		const double sint = sin(theta), tana_pix = 2. * tan(half_angle_) / (r * pixel), t1 = r + 100.;
		if (theta < epsilon || M_PI - theta < epsilon) {
			for (int i = 0; i < pixel; ++i)
				for (int j = 0; j < pixel; ++j) {
					const double k = gsl_hypot(i - 0.5 * pixel + 0.5, j - 0.5 * pixel + 0.5);
					initials_[i * pixel + j] = {0., r, theta < M_PI_2 ? 1e-15 : M_PI - 1e-15, 0., 1., 1., GSL_SIGN(theta - M_PI_2) * tana_pix * k, 0., 0.};
					if (k >= epsilon) {
						if (2 * i <= pixel)
							initials_[i * pixel + j][3] = (theta < M_PI_2 ? 1. : -1.) * acos((j - 0.5 * pixel + 0.5) / k);
						else
							initials_[i * pixel + j][3] = (theta < M_PI_2 ? 1. : -1.) * (2. * M_PI - acos((j - 0.5 * pixel + 0.5) / k));
					}
				}
		} else
			for (int i = 0; i < pixel; ++i)
				for (int j = 0; j < pixel; ++j)
					initials_[i * pixel + j] = {0., r, theta, 0., 1., 1., -tana_pix * (i - 0.5 * pixel + 0.5), -tana_pix * (j - 0.5 * pixel + 0.5) / sint, 0.};
#pragma omp parallel for
		for (int p = pixel * pixel - 1; p >= 0; --p) {
			int status = 0;
			initials_[p][9] = -1.;
#ifdef VIEW_TAU
			Integrator integrator(Metric::functionTau, Metric::jacobian, 2);
#else
			Integrator integrator(Metric::functionHamiltonian, Metric::jacobian, 2);
#endif
			Metric::lightNormalization(initials_[p].data(), 1.);
			Metric::qdq2qp(initials_[p].data());
			while (status <= 0 && initials_[p][8] > t1 && initials_[p][1] > 100)
				status = integrator.Apply(initials_[p].data() + 8, t1, initials_[p].data() + 9, initials_[p].data());
			if (status > 0)
				fmt::print(stderr, "[!] camera::initialize status = {}\n", status);
		}
	}
	int Camera::TraceStar() {
		const double t1 = -1000.;
#pragma omp parallel for
		for (int p = pixel_ * pixel_ - 1; p >= 0; --p) {
			int i = p / pixel_;
			int j = p - i * pixel_;
			double ph[10], last[10]; // TODO: should add t_0 (initials[8]) here
			int status = 0;
			copy(initials_[p].begin(), initials_[p].end(), ph);
#ifdef VIEW_TAU
			Integrator integrator(Metric::functionTau, Metric::jacobian, 2);
#else
			Integrator integrator(Metric::function, Metric::jacobian, 2);
#endif
			while (status <= 0 && ph[8] > t1) {
				copy(ph, ph + 10, last);
				status = integrator.Apply(ph + 8, t1, ph + 9, ph);
				for (auto objP : Object::object_list_)
					if (objP->Hit(ph, last))
						screen_[i][j] = objP->Frequency(ph); // FIXME: if multi objects
				if (screen_[i][j] > epsilon)
					break;
			}
			if (status > 0)
				fmt::print(stderr, "[!] camera::traceStar status = {}\n", status);
		}
		return 0;
	}
	int Camera::Lens() {
		const double t1 = -1000. * r_, pixelPerAngle = 0.5 * pixel_ / half_angle_;
		NumPy rec("lens", {2});
		indicators::BlockProgressBar lensProgressBar{
			indicators::option::ShowElapsedTime{true},
			indicators::option::ShowRemainingTime{true},
			indicators::option::ForegroundColor{indicators::Color(5)},
			indicators::option::PrefixText{"? lens"},
			indicators::option::MaxProgress{pixel_ * pixel_},
			indicators::option::FontStyles{vector<indicators::FontStyle>{indicators::FontStyle::bold}}};
		int progressBarIndex = ProgressBar::bars_.push_back(lensProgressBar);
		for (int i = 0; i < pixel_; ++i)
			for (int j = 0; j < pixel_; ++j) {
				double ph[10], phc[10];
				int status = 0;
				copy(initials_[i * pixel_ + j].begin(), initials_[i * pixel_ + j].end(), ph);
#ifdef VIEW_TAU
				Integrator integrator(Metric::functionTau, Metric::jacobian, 2);
#else
				Integrator integrator(Metric::functionHamiltonian, Metric::jacobian, 2);
#endif
				while (status <= 0 && ph[8] > t1 && ph[1] > 3. && ph[1] < 3.e2)
					status = integrator.Apply(ph + 8, t1, ph + 9, ph);
				if (status > 0)
					fmt::print(stderr, "[!] camera::lens status = {}\n", status);
				if (ph[1] <= 3.)
					rec.save({NAN, NAN});
				else {
					Metric::qp2qdq(ph);
					if (ph[2] < 0)
						ph[2] += M_PI;
					Metric::s2c(ph, phc);
					rec.save({phc[6] * pixelPerAngle, (phc[7] * sin_theta_observer_ - phc[5] * cos_theta_observer_) * pixelPerAngle});
				}
				if (ProgressBar::display_)
					ProgressBar::bars_[progressBarIndex].tick();
			}
		if (ProgressBar::display_)
			ProgressBar::SetComplete(progressBarIndex, "! lens");
		return 0;
	}
	int Camera::Save() {
		for (const vector<double> &line : screen_)
			output_->save(line);
		return 0;
	}
} // namespace SBody
