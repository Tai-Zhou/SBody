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
#include <gsl/gsl_blas.h>
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
	View::View(shared_ptr<Metric> metric, double r, double theta, string file_name) : metric_(move(metric)), r_(r), theta_(theta), sin_theta_observer_(sin(theta)), cos_theta_observer_(cos(theta)), t_final_(-r - 2e4 * 1.), output_(make_unique<NumPy>(file_name, vector<int>({5}))) {}
	int View::InitializePhoton(double *photon, double alpha, double beta) {
		photon[0] = 0.;
		photon[1] = r_;
		photon[5] = 1.;
		photon[8] = 0.;
		photon[9] = 0.;
		if (sin_theta_observer_ < epsilon) {
			const double k = gsl_hypot(alpha, beta);
			if (theta_ < M_PI_2) {
				photon[2] = 1e-15;
				photon[3] = ModBy2Pi(atan2(alpha, -beta));
				photon[6] = -k / gsl_pow_2(r_);
			} else {
				photon[2] = M_PI - 1e-15;
				photon[3] = ModBy2Pi(atan2(alpha, beta));
				photon[6] = k / gsl_pow_2(r_);
			}
			photon[7] = 0.;
		} else {
			photon[2] = theta_;
			photon[3] = 0.,
			photon[6] = beta / gsl_pow_2(r_);
			photon[7] = -alpha / (gsl_pow_2(r_) * sin_theta_observer_);
		}
		return metric_->NormalizeNullGeodesic(photon, 1.);
	}
	int View::TraceStar(Star &star, int ray_number) { // FIXME:!!!!
		const double r_star = star.position_[1], sin_theta_star = sin(star.position_[2]), cos_theta_star = cos(star.position_[2]), sin_phi_star = sin(star.position_[3]), cos_phi_star = cos(star.position_[3]);
		double alpha0 = GSL_POSINF, alpha1 = (r_star + 1.) * sin_theta_star * sin_phi_star, beta0 = GSL_POSINF, beta1 = (r_star + 1.) * (cos_theta_star * sin_theta_observer_ - sin_theta_star * cos_phi_star * cos_theta_observer_), ph[10], cosph, last[10], h;
#ifdef RECORD_TRACE
		vector<double> qdq(12);
		IO::NumPy rec("Trace " + to_string(ray_number), 12);
#endif
		Integrator &&integrator = metric_->GetIntegrator(2);
		while (gsl_hypot(alpha1 - alpha0, beta1 - beta0) > epsilon * (1. + gsl_hypot(alpha1, beta1))) {
			alpha0 = alpha1;
			beta0 = beta1;
			InitializePhoton(ph, alpha0, beta0);
			metric_->LagrangianToHamiltonian(ph);
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
				cosph = (r_star * sin_theta_star * cos_phi_star - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * cos(ph[3])) * sin_theta_observer_ + (r_star * cos_theta_star - ph[1] * GSL_SIGN(ph[2]) * cos(ph[2])) * cos_theta_observer_;
				if (cosph > r_star * relative_accuracy) {
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
						alpha1 += r_star * sin_theta_star * sin_phi_star - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * sin(ph[3]); // TODO:OPT!
						beta1 += (r_star * cos_theta_star - ph[1] * GSL_SIGN(ph[2]) * cos(ph[2])) * sin_theta_observer_ - (r_star * sin_theta_star * cos_phi_star - ph[1] * GSL_SIGN(ph[2]) * sin(ph[2]) * cos(ph[3])) * cos_theta_observer_;
#ifdef VIEW_HAMILTONIAN
						metric_->HamiltonianToLagrangian(ph);
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
		double ph2[9], rec2[sample_number][8], recph[8], interval = M_2PI / sample_number, time_limit;
		copy(ph, ph + 8, ph2);
		ph2[8] = 0.;
		metric_->NormalizeNullGeodesic(ph2, 1.);
#ifdef VIEW_HAMILTONIAN
		metric_->LagrangianToHamiltonian(ph2);
#endif
		h = 1.;
		int status = 0;
		integrator.Reset();
		while (status <= 0 && ph2[1] < 1.e3)
			status = integrator.Apply(ph2 + 8, -t_final_, &h, ph2);
		if (status > 0) {
			fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
			return status;
		}
		copy(ph2, ph2 + 8, recph);
		time_limit = ph2[8];
#ifdef VIEW_HAMILTONIAN
		metric_->HamiltonianToLagrangian(recph);
#endif
		gsl_matrix *coordinate = gsl_matrix_alloc(4, 4), *gmunu = gsl_matrix_alloc(4, 4), *coordinate_gmunu = gsl_matrix_alloc(4, 4);
		gsl_permutation *perm = gsl_permutation_alloc(4);
		star.GetMetricTensor(gmunu);
		star.LocalInertialFrame(coordinate->data);
		const double ph_reg[4] = {1. / ph[4], ph[5] / ph[4], ph[6] / ph[4], ph[7] / ph[4]};
		double ph_coor[4] = {
			star.DotProduct(ph_reg, coordinate->data, 4),
			star.DotProduct(ph_reg, coordinate->data + 4, 4),
			star.DotProduct(ph_reg, coordinate->data + 8, 4),
			star.DotProduct(ph_reg, coordinate->data + 12, 4)};
		for (int i = 1; i < 4; ++i)
			ph_coor[i] /= ph_coor[0];
		const double ph_coor_sph[2] = {acos(ph_coor[3]), atan2(ph_coor[2], ph_coor[1])};
		for (int i = 0; i < sample_number; ++i) {
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
			gsl_matrix_set_zero(coordinate_gmunu);
			gsl_blas_dsymm(CblasLeft, CblasUpper, 1., coordinate, gmunu, 0., coordinate_gmunu);
			int signum;
			gsl_linalg_LU_decomp(coordinate_gmunu, perm, &signum);
			gsl_linalg_LU_svx(coordinate_gmunu, perm, &ph2_view.vector);
			// ph2[4] = 1. / (coordinate[0] + coordinate[4] * ph2_coor_phi[0] + coordinate[8] * ph2_coor_phi[1] + coordinate[12] * ph2_coor_phi[2]);
			// ph2[5] = ph2[4] * (coordinate[1] + coordinate[5] * ph2_coor_phi[0] + coordinate[9] * ph2_coor_phi[1] + coordinate[13] * ph2_coor_phi[2]);
			// ph2[6] = ph2[4] * (coordinate[2] + coordinate[6] * ph2_coor_phi[0] + coordinate[10] * ph2_coor_phi[1] + coordinate[14] * ph2_coor_phi[2]);
			// ph2[7] = ph2[4] * (coordinate[3] + coordinate[7] * ph2_coor_phi[0] + coordinate[11] * ph2_coor_phi[1] + coordinate[15] * ph2_coor_phi[2]);
			ph2[8] = 0.;
			metric_->NormalizeNullGeodesic(ph2, 1.);
#ifdef VIEW_HAMILTONIAN
			metric_->LagrangianToHamiltonian(ph2);
#endif
			h = 1.;
			int status = 0;
			integrator.Reset();
			while (status <= 0 && ph2[8] < time_limit)
				status = integrator.Apply(ph2 + 8, time_limit, &h, ph2);
			if (status > 0) {
				fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
				return status;
			}
			copy(ph2, ph2 + 8, rec2[i]);
#ifdef VIEW_HAMILTONIAN
			metric_->HamiltonianToLagrangian(rec2[i]);
#endif
		}
		const double sin_theta = sin(recph[2]), cos_theta = cos(recph[2]), sin_phi = sin(recph[3]), cos_phi = cos(recph[3]);
		double vphph[3] = {recph[5] * sin_theta * cos_phi + recph[1] * (cos_theta * cos_phi * recph[6] - sin_theta * sin_phi * recph[7]),
						   recph[5] * sin_theta * sin_phi + recph[1] * (cos_theta * sin_phi * recph[6] + sin_theta * cos_phi * recph[7]),
						   recph[5] * cos_theta - recph[1] * sin_theta * recph[6]};
		const double vphph_norm = Norm(vphph);
		for (int j = 0; j < 3; ++j)
			vphph[j] /= vphph_norm;
		double rec3[100][3];
		for (int i = 0; i < sample_number; ++i) {
			const double sin_theta = sin(rec2[i][2]), cos_theta = cos(rec2[i][2]), sin_phi = sin(rec2[i][3]), cos_phi = cos(rec2[i][3]);
			rec3[i][0] = rec2[i][5] * sin_theta * cos_phi + rec2[i][1] * (cos_theta * cos_phi * rec2[i][6] - sin_theta * sin_phi * rec2[i][7]);
			rec3[i][1] = rec2[i][5] * sin_theta * sin_phi + rec2[i][1] * (cos_theta * sin_phi * rec2[i][6] + sin_theta * cos_phi * rec2[i][7]);
			rec3[i][2] = rec2[i][5] * cos_theta - rec2[i][1] * sin_theta * rec2[i][6];
			const double vph_norm = Norm(rec3[i]);
			for (int j = 0; j < 3; ++j)
				rec3[i][j] /= vph_norm;
		}
		double cross_product[3];
		Cross(rec3[0], rec3[sample_number - 1], cross_product);
		double area = Dot(cross_product, vphph);
		for (int i = 1; i < sample_number; ++i) {
			Cross(rec3[i], rec3[i - 1], cross_product);
			area += Dot(cross_product, vphph);
		}
#ifdef VIEW_TAU
		output->Save({alpha1, beta1, star.RedshiftTau(ph), (ph[8] + ph[9]) / Unit::s, abs(area) / (M_2PI * gsl_pow_2(epsilon))});
#else
		output_->Save({alpha1, beta1, star.Redshift(ph), (ph[8] + ph[9]) / Unit::s, abs(area) / (M_2PI * gsl_pow_2(epsilon))});
#endif
		gsl_matrix_free(coordinate);
		gsl_matrix_free(gmunu);
		gsl_matrix_free(coordinate_gmunu);
		gsl_permutation_free(perm);
		return 0;
	}
	int View::Shadow() {
		const double interval = M_2PI / sample_number;
		NumPy rec("shadow", {2});
		double h, rin = 2., rout = 10., rmid = 6., photon[10];
		Integrator &&integrator = metric_->GetIntegrator(2);
		indicators::BlockProgressBar shadowProgressBar{
			indicators::option::ShowElapsedTime{true},
			indicators::option::ShowRemainingTime{true},
			indicators::option::PrefixText{"? Shadow"},
			indicators::option::ForegroundColor{indicators::Color(4)},
			indicators::option::MaxProgress{sample_number},
			indicators::option::FontStyles{vector<indicators::FontStyle>{indicators::FontStyle::bold}}};
		int progressBarIndex = ProgressBar::bars_.push_back(shadowProgressBar);
		for (int i = 0; i < sample_number; ++i) {
			const double angle = i * interval, sina = sin(angle), cosa = cos(angle);
			int status = 0;
			while (rout - rin > epsilon * (rin + rout)) {
				rmid = 0.5 * (rin + rout);
				InitializePhoton(photon, rmid * cosa, rmid * sina);
				h = -1.;
				while (status <= 0 && photon[8] + photon[9] > t_final_) {
					status = integrator.Apply(photon + 9, t_final_, &h, photon);
					if (photon[9] < t_final_ * 1e-8) {
						photon[8] += photon[9];
						photon[9] = 0.;
					}
					if (photon[4] >= 1e6 || photon[5] <= 0 || photon[1] <= 0)
						break;
				}
				if (status > 0) {
					fmt::print(stderr, "[!] view::shadow status = {}\n", status);
					return status;
				}
				if (photon[5] <= 0)
					rout = rmid;
				else
					rin = rmid;
				integrator.Reset();
			}
			rin -= 2. * interval * rmid;
			rout += 2. * interval * rmid;
			rec.Save({rmid * cosa, rmid * sina});
			if (ProgressBar::display_)
				ProgressBar::bars_[progressBarIndex].tick();
		}
		if (ProgressBar::display_)
			ProgressBar::SetComplete(progressBarIndex, "! Shadow");
		return 0;
	}
	Camera::Camera(shared_ptr<Metric> metric, size_t pixel, double half_angle, double r, double theta, string file_name) : View(move(metric), r, theta, file_name), pixel_(pixel), half_angle_(half_angle) {
		screen_ = vector<vector<double>>(pixel, vector<double>(pixel));
		initials_ = vector<array<double, 10>>(pixel * pixel);
		const double pixel_size = 2. * half_angle * r / pixel, t1 = r + 100.;
		for (int i = 0; i < pixel; ++i)
			for (int j = 0; j < pixel; ++j)
				InitializePhoton(initials_[i * pixel + j].data(), pixel_size * (i - 0.5 * pixel + 0.5), pixel_size * (j - 0.5 * pixel + 0.5));
#pragma omp parallel for
		for (int p = pixel * pixel - 1; p >= 0; --p) {
			int status = 0;
			initials_[p][9] = -1.;
			Integrator &&integrator = metric_->GetIntegrator(2);
			metric_->NormalizeNullGeodesic(initials_[p].data(), 1.);
			metric_->LagrangianToHamiltonian(initials_[p].data());
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
			double ph[10], last[10];
			int status = 0;
			copy(initials_[p].begin(), initials_[p].end(), ph);
			Integrator &&integrator = metric_->GetIntegrator(2);
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
				Integrator &&integrator = metric_->GetIntegrator(2);
				while (status <= 0 && ph[8] > t1 && ph[1] > 3. && ph[1] < 3.e2)
					status = integrator.Apply(ph + 8, t1, ph + 9, ph);
				if (status > 0)
					fmt::print(stderr, "[!] camera::lens status = {}\n", status);
				if (ph[1] <= 3.)
					rec.Save({NAN, NAN});
				else {
					metric_->HamiltonianToLagrangian(ph);
					if (ph[2] < 0)
						ph[2] += M_PI;
					SphericalToCartesian(ph, phc, 8);
					rec.Save({phc[6] * pixelPerAngle, (phc[7] * sin_theta_observer_ - phc[5] * cos_theta_observer_) * pixelPerAngle});
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
			output_->Save(line);
		return 0;
	}
} // namespace SBody
