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
#include <memory>

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
	View::View(shared_ptr<Metric> metric, double r, double theta, string file_name) : metric_(move(metric)), r_(r), theta_(theta), sin_theta_observer_(sin(theta)), cos_theta_observer_(cos(theta)), t_final_(-r - 2e4 * 1.), output_(make_unique<NumPy>(file_name, vector<int>({8}))) {
		bar_ = make_unique<indicators::BlockProgressBar>(indicators::option::ShowElapsedTime{true},
														 indicators::option::ShowRemainingTime{true},
														 indicators::option::ForegroundColor{indicators::Color(4)},
														 indicators::option::FontStyles{vector<indicators::FontStyle>{indicators::FontStyle::bold}});
	}
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
		double alpha0 = GSL_POSINF, alpha1 = (r_star + 1.) * sin_theta_star * sin_phi_star, beta0 = GSL_POSINF, beta1 = (r_star + 1.) * (cos_theta_star * sin_theta_observer_ - sin_theta_star * cos_phi_star * cos_theta_observer_), cosph, h;
		GslBlock collector;
		gsl_vector *photon = collector.VectorAlloc(10), *last = collector.VectorAlloc(10);
#ifdef RECORD_TRACE
		vector<double> qdq(12);
		IO::NumPy rec("Trace " + to_string(ray_number), 12);
#endif
		Integrator &&integrator = metric_->GetIntegrator(2);
		while (gsl_hypot(alpha1 - alpha0, beta1 - beta0) > epsilon * (1. + gsl_hypot(alpha1, beta1))) {
			alpha0 = alpha1;
			beta0 = beta1;
			InitializePhoton(photon->data, alpha0, beta0);
			metric_->LagrangianToHamiltonian(photon->data);
			int status = 0, fixed = 0;
			h = -1.;
			while (status <= 0) {
				if (status == -1) {
					h *= 0.5;
					fixed = 0;
				}
				gsl_vector_memcpy(last, photon);
				if (fixed)
					status = integrator.ApplyFixedStep(gsl_vector_ptr(photon, 9), h, photon->data);
				else
					status = integrator.Apply(gsl_vector_ptr(photon, 9), t_final_, &h, photon->data);
				cosph = (r_star * sin_theta_star * cos_phi_star - gsl_vector_get(photon, 1) * GSL_SIGN(gsl_vector_get(photon, 2)) * sin(gsl_vector_get(photon, 2)) * cos(gsl_vector_get(photon, 3))) * sin_theta_observer_ + (r_star * cos_theta_star - gsl_vector_get(photon, 1) * GSL_SIGN(gsl_vector_get(photon, 2)) * cos(gsl_vector_get(photon, 2))) * cos_theta_observer_;
				if (cosph > r_star * relative_accuracy) {
					gsl_vector_memcpy(photon, last);
					h *= 0.3;
					fixed = 1;
					integrator.Reset();
				} else {
					if (gsl_vector_get(photon, 9) < t_final_ * 1e-8) {
						gsl_vector_set(photon, 8, gsl_vector_get(photon, 8) + gsl_vector_get(photon, 9));
						gsl_vector_set(photon, 9, 0);
					}
#ifdef RECORD_TRACE
					copy(photon, photon + 8, qdq.begin());
					metric_->HamiltonianToLagrangian(qdq.data());
					qdq[8] = photon[8] + photon[9];
					qdq[9] = metric::energy(qdq.data());
					qdq[10] = metric::angularMomentum(qdq.data());
					qdq[11] = metric::carter(qdq.data(), 0.);
					rec.save(qdq);
#endif
					if (cosph >= 0) {
						alpha1 += r_star * sin_theta_star * sin_phi_star - gsl_vector_get(photon, 1) * GSL_SIGN(gsl_vector_get(photon, 2)) * sin(gsl_vector_get(photon, 2)) * sin(gsl_vector_get(photon, 3)); // TODO:OPT!
						beta1 += (r_star * cos_theta_star - gsl_vector_get(photon, 1) * GSL_SIGN(gsl_vector_get(photon, 2)) * cos(gsl_vector_get(photon, 2))) * sin_theta_observer_ - (r_star * sin_theta_star * cos_phi_star - gsl_vector_get(photon, 1) * GSL_SIGN(gsl_vector_get(photon, 2)) * sin(gsl_vector_get(photon, 2)) * cos(gsl_vector_get(photon, 3))) * cos_theta_observer_;
						metric_->HamiltonianToLagrangian(photon->data);
						break;
					}
				}
			}
			if (status > 0) {
				fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
				return status;
			}
		}
		double ph2[9], rec2[sample_number][3], rec3[sample_number][3], recph[3], interval = M_2PI / sample_number, time_limit = 0;
		gsl_vector_set(photon, 0, 0.);
		copy(photon->data, photon->data + 8, ph2);
		metric_->NormalizeNullGeodesic(ph2, 1.);
		metric_->LagrangianToHamiltonian(ph2);
		h = 1.;
		int status = 0, signum;
		integrator.Reset();
		while (status <= 0 && ph2[1] < 1.e3)
			status = integrator.Apply(&time_limit, -t_final_, &h, ph2);
		if (status > 0) {
			fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
			return status;
		}
		metric_->HamiltonianToLagrangian(ph2);
		SphericalToCartesian(ph2);
		copy(ph2 + 5, ph2 + 8, recph);
		const double recph_norm = Norm(recph);
		for (int i = 0; i < 3; ++i)
			recph[i] /= recph_norm;
		gsl_matrix *coordinate = collector.MatrixAlloc(4, 4), *gmunu = collector.MatrixAlloc(4, 4), *coordinate_gmunu = collector.MatrixAlloc(4, 4), *coordinate_static = collector.MatrixCalloc(4, 4);
		gsl_permutation *permutation = collector.PermutationAlloc(4);
		star.MetricTensor(gmunu);
		star.LocalInertialFrame(coordinate);
		gsl_blas_dsymm(CblasRight, CblasUpper, 1., gmunu, coordinate, 0., coordinate_gmunu);
		gsl_linalg_LU_decomp(coordinate_gmunu, permutation, &signum);
		gsl_matrix_set(coordinate_static, 0, 0, sqrt(-1. / gmunu->data[0]));
		gsl_matrix_set(coordinate_static, 1, 1, sqrt(1. / gmunu->data[5]));
		gsl_matrix_set(coordinate_static, 2, 2, sqrt(1. / gmunu->data[10]));
		gsl_matrix_set(coordinate_static, 3, 3, sqrt(1. / gmunu->data[15]));
		const double ph_reg[4] = {1., gsl_vector_get(photon, 5), gsl_vector_get(photon, 6), gsl_vector_get(photon, 7)};
		double ph_static[4] = {
			star.DotProduct(ph_reg, coordinate_static->data, 4),
			star.DotProduct(ph_reg, coordinate_static->data + 4, 4),
			star.DotProduct(ph_reg, coordinate_static->data + 8, 4),
			star.DotProduct(ph_reg, coordinate_static->data + 12, 4)};
		for (int i = 1; i < 4; ++i)
			ph_static[i] /= ph_static[0];
		double ph_coor[4] = {
			star.DotProduct(ph_reg, coordinate->data, 4),
			star.DotProduct(ph_reg, coordinate->data + 4, 4),
			star.DotProduct(ph_reg, coordinate->data + 8, 4),
			star.DotProduct(ph_reg, coordinate->data + 12, 4)};
		for (int i = 1; i < 4; ++i)
			ph_coor[i] /= ph_coor[0];
		CartesianToSpherical(ph_coor, 4);
		for (int i = 0; i < sample_number; ++i) {
			const double angle = i * interval, sin_angle = sin(angle), cos_angle = cos(angle);
			copy(photon->data, photon->data + 4, ph2);
			ph2[4] = 1.;
			ph2[5] = cos_angle * sin_epsilon;
			ph2[6] = sin_angle * sin_epsilon;
			ph2[7] = cos_epsilon;
			RotateAroundAxis(ph2 + 5, 1, ph_coor[2]);
			RotateAroundAxis(ph2 + 5, 2, ph_coor[3]);
			gsl_vector_view ph2_view = gsl_vector_view_array(ph2 + 4, 4);
			gsl_linalg_LU_svx(coordinate_gmunu, permutation, &ph2_view.vector);
			ph2[8] = 0.;
			metric_->NormalizeNullGeodesic(ph2);
			rec3[i][0] = star.DotProduct(ph2 + 4, coordinate_static->data + 4, 4);
			rec3[i][1] = star.DotProduct(ph2 + 4, coordinate_static->data + 8, 4);
			rec3[i][2] = star.DotProduct(ph2 + 4, coordinate_static->data + 12, 4);
			const double vph_norm = star.DotProduct(ph2 + 4, coordinate_static->data, 4);
			for (int j = 0; j < 3; ++j)
				rec3[i][j] /= vph_norm;
			metric_->LagrangianToHamiltonian(ph2);
			h = 1.;
			int status = 0;
			integrator.Reset();
			while (status <= 0 && ph2[8] < time_limit)
				status = integrator.Apply(ph2 + 8, time_limit, &h, ph2);
			if (status > 0) {
				fmt::print(stderr, "[!] view::traceStar status = {}\n", status);
				return status;
			}
			metric_->HamiltonianToLagrangian(ph2);
			SphericalToCartesian(ph2);
			copy(ph2 + 5, ph2 + 8, rec2[i]);
			const double rec2_norm = Norm(rec2[i]); // TODO: == ph2[4]?
			for (int j = 0; j < 3; ++j)
				rec2[i][j] /= rec2_norm;
		}
		if (ray_number == 2500) {
			NumPy cone_record("cone_record", {3});
			cone_record.Save(recph + 5, 3);
			for (int i = 0; i < sample_number; ++i)
				cone_record.Save(rec2[i], 3);
		}
		double area = DotCross(recph + 5, rec2[0], rec2[sample_number - 1]);
		double area2 = DotCross(ph_static + 1, rec3[0], rec3[sample_number - 1]);
		for (int i = 1; i < sample_number; ++i) {
			area += DotCross(recph + 5, rec2[i], rec2[i - 1]);
			area2 += DotCross(ph_static + 1, rec3[i], rec3[i - 1]);
		}
#ifdef VIEW_TAU
		output->Save({alpha1, beta1, star.RedshiftTau(photon), (photon[8] + photon[9]) / Unit::s, abs(area) / (M_2PI * gsl_pow_2(epsilon)), sqrt(-1. / gmunu->data[0]), abs(area2) / (M_2PI * gsl_pow_2(epsilon))});
#else
		output_->Save({alpha1, beta1, star.Redshift(photon->data), (gsl_vector_get(photon, 8) + gsl_vector_get(photon, 9)) / Unit::s, abs(area) / (M_2PI * gsl_pow_2(epsilon)), sqrt(-1. / gmunu->data[0]), abs(area2) / (M_2PI * gsl_pow_2(epsilon))});
#endif
		return 0;
	}
	int View::Shadow() {
		const double interval = M_2PI / sample_number;
		NumPy rec("shadow", {2});
		double h, rin = 2., rout = 10., rmid = 6., photon[10];
		Integrator &&integrator = metric_->GetIntegrator(2);
		bar_->set_option(indicators::option::MaxProgress(sample_number));
		bar_->set_option(indicators::option::PrefixText("? Shadow"));
		int progressBarIndex = ProgressBar::bars_.push_back(*bar_);
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
		for (size_t i = 0; i < pixel; ++i)
			for (size_t j = 0; j < pixel; ++j)
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
		bar_->set_option(indicators::option::MaxProgress(pixel_ * pixel_));
		bar_->set_option(indicators::option::PrefixText("? Lens"));
		int progressBarIndex = ProgressBar::bars_.push_back(*bar_);
		for (size_t i = 0; i < pixel_; ++i)
			for (size_t j = 0; j < pixel_; ++j) {
				double ph[10];
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
					SphericalToCartesian(ph);
					rec.Save({ph[6] * pixelPerAngle, (ph[7] * sin_theta_observer_ - ph[5] * cos_theta_observer_) * pixelPerAngle});
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
