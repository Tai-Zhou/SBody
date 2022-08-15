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
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

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
	View::View(double r, double theta, string file_name) : r_(r), theta_(theta), sin_theta_observer_(sin(theta)), cos_theta_observer_(cos(theta)), t_final_(-r - 2e4 * 1.), output_(make_unique<NumPy>(file_name, vector<int>({5}))) {}
	int View::TraceStar(Star &star, int ray_number) { // FIXME:!!!!
		const double rs = star.pos[1], thetas = star.pos[2], phis = star.pos[3], sints = sin(thetas), costs = cos(thetas), sinps = sin(phis), cosps = cos(phis);
		double alpha0 = GSL_POSINF, alpha1 = (rs + 1.) * sints * sin(phis), beta0 = GSL_POSINF, beta1 = (rs + 1.) * (costs * sin_theta_observer_ - sints * cosps * cos_theta_observer_), ph[10], cosph, last[10], h;
#ifdef RECORD_TRACE
		vector<double> qdq(12);
		IO::NumPy rec("Trace " + to_string(ray_number), 12);
#endif
#ifdef VIEW_TAU
		Integrator integrator(metric::functionTau, metric::jacobian, 2);
#elif VIEW_HAMILTONIAN
		Integrator integrator(metric::functionHamiltonian, metric::jacobian, 2);
#else
		Integrator integrator(metric::function, metric::jacobian, 2);
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
			metric::lightNormalization(ph, 1.);
#ifdef VIEW_HAMILTONIAN
			metric::qdq2qp(ph);
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
					metric::qp2qdq(qdq.data());
#endif
					qdq[8] = ph[8] + ph[9];
					qdq[9] = metric::energy(qdq.data());
					qdq[10] = metric::angularMomentum(qdq.data());
					qdq[11] = metric::carter(qdq.data(), 0.);
					rec.save(qdq);
#endif
					if (cosph >= 0) {
						alpha1 += rs * sints * sinps - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * sin(ph[3]); // TODO:OPT!
						beta1 += (rs * costs - ph[1] * GSL_SIGN(ph[2]) * cos(ph[2])) * sin_theta_observer_ - (rs * sints * cosps - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * cos(ph[3])) * cos_theta_observer_;
#ifdef VIEW_HAMILTONIAN
						metric::qp2qdq(ph);
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
		int n = 100;
		double coordinate[16], ph2[10], rec2[100][8], recph[8], interval = M_2PI / n;
		metric::Schwarzschild::LocalInertialFrame(star.pos, coordinate);
		const double ph_reg[4] = {1. / ph[4], ph[5] / ph[4], ph[6] / ph[4], ph[7] / ph[4]};
		double ph_coor[4] = {
			metric::Schwarzschild::dot(star.pos, ph_reg, coordinate, 4),
			metric::Schwarzschild::dot(star.pos, ph_reg, coordinate + 4, 4),
			metric::Schwarzschild::dot(star.pos, ph_reg, coordinate + 8, 4),
			metric::Schwarzschild::dot(star.pos, ph_reg, coordinate + 12, 4)};
		for (int i = 1; i < 4; ++i)
			ph_coor[i] /= ph_coor[0];
		const double ph_coor_sph[2] = {acos(ph_coor[3]), atan2(ph_coor[2], ph_coor[1])};
		for (int i = 0; i < n; ++i) {
			const double angle = i * interval, sin_angle = sin(angle), cos_angle = cos(angle);
			const double ph2_coor[3] = {
				cos_angle * sin_epsilon,
				sin_angle * sin_epsilon,
				cos_epsilon};
			const double ph2_coor_theta[3] = {
				ph2_coor[0] * cos(ph_coor_sph[0]) + ph2_coor[2] * sin(ph_coor_sph[0]),
				ph2_coor[1],
				ph2_coor[2] * cos(ph_coor_sph[0]) - ph2_coor[0] * sin(ph_coor_sph[0])};
			// const double ph2_coor_phi[3] = {
			//	ph2_coor_theta[0] * cos(ph_coor_sph[1]) - ph2_coor_theta[1] * sin(ph_coor_sph[1]),
			//	ph2_coor_theta[1] * cos(ph_coor_sph[1]) + ph2_coor_theta[0] * sin(ph_coor_sph[1]),
			//	ph2_coor_theta[2]};
			copy(ph, ph + 4, ph2);
			ph2[4] = 1.;
			ph2[5] = ph2_coor_theta[0] * cos(ph_coor_sph[1]) - ph2_coor_theta[1] * sin(ph_coor_sph[1]);
			ph2[6] = ph2_coor_theta[1] * cos(ph_coor_sph[1]) + ph2_coor_theta[0] * sin(ph_coor_sph[1]);
			ph2[7] = ph2_coor_theta[2];
			gsl_vector_view ph2_view = gsl_vector_view_array(ph2 + 4, 4);
			// coordinate * gmunu * ph2 = ph_coor_phi
			gsl_matrix *coordinate_gmunu = gsl_matrix_alloc(4, 4);
			gsl_matrix_set_zero(coordinate_gmunu);
			double gmunu[16];
			metric::Schwarzschild::gmunu(star.pos, gmunu);
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j)
					for (int k = 0; k < 4; ++k)
						coordinate_gmunu->data[i * 4 + j] += coordinate[i * 4 + k] * gmunu[k * 4 + j];
			gsl_permutation *perm = gsl_permutation_alloc(4);
			int signum;
			gsl_linalg_LU_decomp(coordinate_gmunu, perm, &signum);
			gsl_linalg_LU_svx(coordinate_gmunu, perm, &ph2_view.vector);
			// ph2[4] = 1. / (coordinate[0] + coordinate[4] * ph2_coor_phi[0] + coordinate[8] * ph2_coor_phi[1] + coordinate[12] * ph2_coor_phi[2]);
			// ph2[5] = ph2[4] * (coordinate[1] + coordinate[5] * ph2_coor_phi[0] + coordinate[9] * ph2_coor_phi[1] + coordinate[13] * ph2_coor_phi[2]);
			// ph2[6] = ph2[4] * (coordinate[2] + coordinate[6] * ph2_coor_phi[0] + coordinate[10] * ph2_coor_phi[1] + coordinate[14] * ph2_coor_phi[2]);
			// ph2[7] = ph2[4] * (coordinate[3] + coordinate[7] * ph2_coor_phi[0] + coordinate[11] * ph2_coor_phi[1] + coordinate[15] * ph2_coor_phi[2]);
			ph2[8] = 0.;
			ph2[9] = 0.;
			metric::lightNormalization(ph2, 1.);
#ifdef VIEW_HAMILTONIAN
			metric::qdq2qp(ph2);
#endif
			h = 1.;
			int status = 0;
			integrator.Reset();
			while (status <= 0 && ph2[1] < 1000.) {
				status = integrator.Apply(ph2 + 9, -t_final_, &h, ph2);
				if (ph2[9] > t_final_ * -1e-8) {
					ph2[8] += ph[9];
					ph2[9] = 0;
				}
			}
			if (status > 0) {
				fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
				return status;
			}
			copy(ph2, ph2 + 8, rec2[i]);
#ifdef VIEW_HAMILTONIAN
			metric::qp2qdq(rec2[i]);
#endif
		}
		copy(ph, ph + 8, ph2);
		ph2[8] = 0.;
		ph2[9] = 0.;
		metric::lightNormalization(ph2, 1.);
#ifdef VIEW_HAMILTONIAN
		metric::qdq2qp(ph2);
#endif
		h = 1.;
		int status = 0;
		integrator.Reset();
		while (status <= 0 && ph2[1] < 1e15) {
			status = integrator.Apply(ph2 + 9, -t_final_, &h, ph2);
			if (ph2[9] > t_final_ * -1e-8) {
				ph2[8] += ph[9];
				ph2[9] = 0;
			}
		}
		if (status > 0) {
			fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
			return status;
		}
		copy(ph2, ph2 + 8, recph);
#ifdef VIEW_HAMILTONIAN
		metric::qp2qdq(recph);
#endif
		const double sin_theta = sin(recph[2]), cos_theta = cos(recph[2]), sin_phi = sin(recph[3]), cos_phi = cos(recph[3]);
		double vphph[3] = {recph[5] * sin_theta * cos_phi + recph[1] * (cos_theta * cos_phi * recph[6] - sin_theta * sin_phi * recph[7]),
						   recph[5] * sin_theta * sin_phi + recph[1] * (cos_theta * sin_phi * recph[6] + sin_theta * cos_phi * recph[7]),
						   recph[5] * cos_theta - recph[1] * sin_theta * recph[6]};
		const double vphph_norm = Norm(vphph);
		for (int j = 0; j < 3; ++j)
			vphph[j] /= vphph_norm;
		double final_angle[100];
		for (int i = 0; i < n; ++i) {
			const double sin_theta = sin(rec2[i][2]), cos_theta = cos(rec2[i][2]), sin_phi = sin(rec2[i][3]), cos_phi = cos(rec2[i][3]);
			double vph[3] = {rec2[i][5] * sin_theta * cos_phi + rec2[i][1] * (cos_theta * cos_phi * rec2[i][6] - sin_theta * sin_phi * rec2[i][7]),
							 rec2[i][5] * sin_theta * sin_phi + rec2[i][1] * (cos_theta * sin_phi * rec2[i][6] + sin_theta * cos_phi * rec2[i][7]),
							 rec2[i][5] * cos_theta - rec2[i][1] * sin_theta * rec2[i][6]};
			const double vph_norm = Norm(vph);
			for (int j = 0; j < 3; ++j)
				vph[j] /= vph_norm;
			double cross_product[3];
			Cross(vph, vphph, cross_product);
			final_angle[i] = Norm(cross_product);
		}
		double area = final_angle[0] * final_angle[99];
		for (int i = 1; i < n; ++i)
			area += final_angle[i] * final_angle[i - 1];
		area /= n;
#ifdef VIEW_TAU
		output->save({alpha1, beta1, star.RedshiftTau(ph), (ph[8] + ph[9]) / Unit::s, area / gsl_pow_2(epsilon)});
#else
		output_->save({alpha1, beta1, star.Redshift(ph), (ph[8] + ph[9]) / Unit::s, area / gsl_pow_2(epsilon)});
#endif
		return 0;
	}
	int View::Shadow(int n) {
		const double interval = M_2PI / n;
		NumPy rec("shadow " + to_string(metric::a) + "," + to_string(metric::l), {2});
		double h, rin = 2., rout = 10., rmid = 6., ph[10];
#ifdef VIEW_TAU
		Integrator integrator(metric::functionTau, metric::jacobian, 2);
#else
		Integrator integrator(metric::function, metric::jacobian, 2);
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
				metric::lightNormalization(ph, 1.);
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
			Integrator integrator(metric::functionTau, metric::jacobian, 2);
#else
			Integrator integrator(metric::functionHamiltonian, metric::jacobian, 2);
#endif
			metric::lightNormalization(initials_[p].data(), 1.);
			metric::qdq2qp(initials_[p].data());
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
			Integrator integrator(metric::functionTau, metric::jacobian, 2);
#else
			Integrator integrator(metric::function, metric::jacobian, 2);
#endif
			while (status <= 0 && ph[8] > t1) {
				copy(ph, ph + 10, last);
				status = integrator.Apply(ph + 8, t1, ph + 9, ph);
				for (auto objP : Object::object_list_)
					if (objP->Hit(ph, last))
						screen_[i][j] = objP->Redshift(ph); // FIXME: if multi objects
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
				Integrator integrator(metric::functionTau, metric::jacobian, 2);
#else
				Integrator integrator(metric::functionHamiltonian, metric::jacobian, 2);
#endif
				while (status <= 0 && ph[8] > t1 && ph[1] > 3. && ph[1] < 3.e2)
					status = integrator.Apply(ph + 8, t1, ph + 9, ph);
				if (status > 0)
					fmt::print(stderr, "[!] camera::lens status = {}\n", status);
				if (ph[1] <= 3.)
					rec.save({NAN, NAN});
				else {
					metric::qp2qdq(ph);
					if (ph[2] < 0)
						ph[2] += M_PI;
					metric::s2c(ph, phc);
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
