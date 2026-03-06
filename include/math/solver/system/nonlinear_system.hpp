/**
 * @file nonlinear_system.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_NONLINEAR_SYSTEM_HPP
#define SBODY_NONLINEAR_SYSTEM_HPP

#include <cstddef>
#include <numeric>
#include <string>
#include <vector>

#include <boost/array.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fmt/core.h>

#include "infra/log.hpp"
#include "math/geometry.hpp"
#include "math/scalar.hpp"

namespace sbody {
	template <typename T, std::size_t N>
	int CoordinateOrthogonalization(const boost::numeric::ublas::bounded_vector<T, N> &x, boost::numeric::ublas::bounded_matrix<T, N, N> &coordinate) {
		namespace ublas = boost::numeric::ublas;
		std::vector<ublas::matrix_column<ublas::bounded_matrix<T, N, N>>> coordinate_columns;
		coordinate_columns.reserve(N);
		for (int i = 0; i < N; ++i)
			coordinate_columns.emplace_back(ublas::matrix_column<ublas::bounded_matrix<T, N, N>>(coordinate, i));
		int column_idx = 1;
		coordinate_columns[0] = x / ublas::norm_2(x);
		for (int coordinate_idx = 0; coordinate_idx < N; ++coordinate_idx) {
			coordinate_columns[column_idx] = ublas::zero_vector<T>(N);
			coordinate_columns[column_idx][coordinate_idx] = 1.;
			for (int i = 0; i < column_idx; ++i) {
				T dot_product = ublas::inner_prod(coordinate_columns[i], coordinate_columns[column_idx]);
				coordinate_columns[column_idx] -= coordinate_columns[i] * dot_product;
			}
			T base_norm = ublas::norm_2(coordinate_columns[column_idx]);
			if (base_norm > boost::math::tools::root_epsilon<T>()) {
				coordinate_columns[column_idx] /= base_norm;
				if (++column_idx >= N)
					return StatusCode::SUCCESS;
			}
		}
		return StatusCode::FAILURE;
	}

	template <typename T>
	bool TrustRegionStep(T rho, T step_norm, T &trust_region_radius, T max_radius) {
		if (rho < 0.25)
			trust_region_radius *= 0.25;
		else if (rho > 0.75 && step_norm > 0.99 * trust_region_radius)
			trust_region_radius = std::min(2. * trust_region_radius, max_radius);
		return rho > 0.1;
	}

	template <typename T, std::size_t N, typename Parameter>
	class NonlinearSystemSolver {
	  protected:
		boost::numeric::ublas::bounded_vector<T, N> x_;
		boost::numeric::ublas::bounded_vector<T, N> f_;
		boost::numeric::ublas::bounded_vector<T, N> dx_;
		int (*function_)(const boost::numeric::ublas::bounded_vector<T, N> &, boost::numeric::ublas::bounded_vector<T, N> &, Parameter &);
		Parameter &function_parameters_;

	  public:
		NonlinearSystemSolver(int (*function)(const boost::numeric::ublas::bounded_vector<T, N> &, boost::numeric::ublas::bounded_vector<T, N> &, Parameter &), Parameter &params) : function_(function), function_parameters_(params) {}
		virtual int Iterate() = 0;
		int Solve(T absolute_accuracy, int max_iteration = 256) {
			for (; max_iteration > 0; --max_iteration) {
				if (boost::numeric::ublas::norm_1(f_) <= absolute_accuracy)
					return StatusCode::SUCCESS;
				try {
					if (int status = Iterate(); status != StatusCode::SUCCESS)
						return status;
				} catch (...) {
					return StatusCode::NUMERIC_ERROR;
				}
				if (std::any_of(x_.begin(), x_.end(), [](T x_i) { return !std::isfinite(x_i); }))
					return StatusCode::NUMERIC_ERROR;
			}
			return StatusCode::MAX_ITERATIONS_EXCEEDED;
		}
		int Solve(T absolute_accuracy, T relative_accuracy, int max_iteration = 256) {
			for (; max_iteration > 0; --max_iteration) {
				if (boost::numeric::ublas::norm_1(f_) <= absolute_accuracy) {
					int status = StatusCode::SUCCESS;
					for (int i = 0; i < N; ++i)
						if (dx_(i) > std::abs(x_(i)) * relative_accuracy + absolute_accuracy) {
							status = StatusCode::CONTINUE;
							break;
						}
					if (status == StatusCode::SUCCESS)
						return StatusCode::SUCCESS;
				}
				try {
					if (int status = Iterate(); status != StatusCode::SUCCESS)
						return status;
				} catch (...) {
					return StatusCode::NUMERIC_ERROR;
				}
				if (std::any_of(x_.begin(), x_.end(), [](T x) { return !std::isfinite(x); }))
					return StatusCode::NUMERIC_ERROR;
			}
			return StatusCode::MAX_ITERATIONS_EXCEEDED;
		}
		const boost::numeric::ublas::bounded_vector<T, N> &Root() {
			return x_;
		}
		const boost::numeric::ublas::bounded_vector<T, N> &Value() {
			return f_;
		}
		const boost::numeric::ublas::bounded_vector<T, N> &StepSize() {
			return dx_;
		}
		int Eval() {
			return function_(x_, f_, function_parameters_);
		}
		int Eval(const boost::numeric::ublas::bounded_vector<T, N> &x, boost::numeric::ublas::bounded_vector<T, N> &f) {
			return function_(x, f, function_parameters_);
		}
		int ScaleX() {
			int status;
			for (status = Eval(x_, f_); status == StatusCode::SCALING_REQUIRED; status = Eval(x_, f_))
				x_ *= f_(0);
			return status;
		}
		int ScaleX(bool &scaled) {
			int status;
			for (status = Eval(x_, f_); status == StatusCode::SCALING_REQUIRED; status = Eval(x_, f_)) {
				x_ *= f_(0);
				scaled = true;
			}
			return status;
		}
		int ScaleX(boost::numeric::ublas::bounded_vector<T, N> &x, boost::numeric::ublas::bounded_vector<T, N> &f) {
			int status;
			for (status = Eval(x, f); status == StatusCode::SCALING_REQUIRED; status = Eval(x, f))
				x *= f(0);
			return status;
		}
		int ScaleX(boost::numeric::ublas::bounded_vector<T, N> &x, boost::numeric::ublas::bounded_vector<T, N> &f, bool &scaled) {
			int status;
			for (status = Eval(x, f); status == StatusCode::SCALING_REQUIRED; status = Eval(x, f)) {
				x *= f(0);
				scaled = true;
			}
			return status;
		}
		int TwoSidedJacobian(const boost::numeric::ublas::bounded_vector<T, N> &x, const boost::numeric::ublas::bounded_vector<T, N> &f, T relative_epsilon, boost::numeric::ublas::bounded_matrix<T, N, N> &jacobian) {
			namespace ublas = boost::numeric::ublas;
			int status = StatusCode::SUCCESS;
			ublas::bounded_vector<T, N> x1(x), f_trial_plus, f_trial_minus;
			for (std::size_t j = 0; j < N; ++j) {
				const T x_j = x(j);
				// auto jacobian_column_j = ublas::matrix_column(jacobian, j);
				T dx, dx0 = std::abs(x_j) >= 1. ? relative_epsilon * x_j : std::copysign(relative_epsilon, x_j), dx_limit = std::abs(x_j) * boost::math::tools::epsilon<T>();
				for (dx = dx0; std::abs(dx) > dx_limit; dx *= 0.5) {
					x1(j) = x_j + dx;
					if (status = Eval(x1, f_trial_plus); status == StatusCode::SUCCESS)
						break;
				}
				if (status != StatusCode::SUCCESS) {
					for (dx = dx0; std::abs(dx) > dx_limit; dx *= 0.5) {
						x1(j) = x_j - dx;
						if (status = Eval(x1, f_trial_minus); status == StatusCode::SUCCESS)
							break;
					}
					if (status != StatusCode::SUCCESS)
						return StatusCode::FAILURE;
					ublas::matrix_column(jacobian, j) = (f - f_trial_minus) / dx;
				} else {
					x1(j) = x_j - dx;
					if (status = Eval(x1, f_trial_minus); status != StatusCode::SUCCESS)
						ublas::matrix_column(jacobian, j) = (f_trial_plus - f) / dx;
					else
						ublas::matrix_column(jacobian, j) = 0.5 * (f_trial_plus - f_trial_minus) / dx;
				}
				x1(j) = x_j;
			}
			return StatusCode::SUCCESS;
		}
		int FourPointJacobian(const boost::numeric::ublas::bounded_vector<T, N> &x, const boost::numeric::ublas::bounded_vector<T, N> &f, T relative_epsilon, boost::numeric::ublas::bounded_matrix<T, N, N> &jacobian) {
			return StatusCode::NOT_IMPLEMENTED;
			namespace ublas = boost::numeric::ublas;
			int status = StatusCode::SUCCESS;
			ublas::bounded_vector<T, N> x1(x), f_trial_plus, f_trial_minus;
			for (std::size_t j = 0; j < N; ++j) {
				const T x_j = x(j);
				// auto jacobian_column_j = ublas::matrix_column(jacobian, j);
				T dx, dx0 = std::abs(x_j) >= 1. ? relative_epsilon * x_j : std::copysign(relative_epsilon, x_j), dx_limit = std::abs(x_j) * boost::math::tools::epsilon<T>();
				for (dx = dx0; std::abs(dx) > dx_limit; dx *= 0.5) {
					x1(j) = x_j + dx;
					if (status = Eval(x1, f_trial_plus); status == StatusCode::SUCCESS)
						break;
				}
				if (status != StatusCode::SUCCESS) {
					for (dx = dx0; std::abs(dx) > dx_limit; dx *= 0.5) {
						x1(j) = x_j - dx;
						if (status = Eval(x1, f_trial_minus); status == StatusCode::SUCCESS)
							break;
					}
					if (status != StatusCode::SUCCESS)
						return StatusCode::FAILURE;
					ublas::matrix_column(jacobian, j) = (f - f_trial_minus) / dx;
				} else {
					x1(j) = x_j - dx;
					if (status = Eval(x1, f_trial_minus); status != StatusCode::SUCCESS)
						ublas::matrix_column(jacobian, j) = (f_trial_plus - f) / dx;
					else
						ublas::matrix_column(jacobian, j) = 0.5 * (f_trial_plus - f_trial_minus) / dx;
				}
				x1(j) = x_j;
			}
			return StatusCode::SUCCESS;
		}
		int TwoSidedDirectionalJacobian(const boost::numeric::ublas::bounded_vector<T, N> &x, const boost::numeric::ublas::bounded_vector<T, N> &direction, const boost::numeric::ublas::bounded_vector<T, N> &f, T relative_epsilon, boost::numeric::ublas::bounded_matrix<T, N, N> &jacobian) {
			namespace ublas = boost::numeric::ublas;
			int status = StatusCode::SUCCESS;
			ublas::bounded_vector<T, N> x1, f_trial_plus, f_trial_minus;
			auto coordinate = ublas::bounded_matrix<T, N, N>();
			if (status = CoordinateOrthogonalization(direction, coordinate); status != StatusCode::SUCCESS)
				return status;
			jacobian.clear();
			T x_norm = ublas::norm_2(x), dx_limit = x_norm * boost::math::tools::epsilon<T>();
			for (std::size_t j = 0; j < N; ++j) {
				auto coordinate_column_j = ublas::matrix_column(coordinate, j);
				T dot_product = ublas::inner_prod(x, coordinate_column_j);
				T dx, dx0 = relative_epsilon * std::copysign(std::max(T{1.}, x_norm), dot_product);
				for (dx = dx0; std::abs(dx) > dx_limit; dx *= 0.5) {
					x1 = x + dx * coordinate_column_j;
					if (status = Eval(x1, f_trial_plus); status == StatusCode::SUCCESS)
						break;
				}
				if (status != StatusCode::SUCCESS) {
					for (dx = dx0; std::abs(dx) > dx_limit; dx *= 0.5) {
						x1 = x - dx * coordinate_column_j;
						if (status = Eval(x1, f_trial_minus); status == StatusCode::SUCCESS)
							break;
					}
					if (status != StatusCode::SUCCESS)
						return StatusCode::FAILURE;
					jacobian += ublas::outer_prod((f - f_trial_minus) / dx, coordinate_column_j);
				} else {
					x1 = x - dx * coordinate_column_j;
					if (status = Eval(x1, f_trial_minus); status != StatusCode::SUCCESS) {
						jacobian += ublas::outer_prod((f_trial_plus - f) / dx, coordinate_column_j);
					} else
						jacobian += ublas::outer_prod((f_trial_plus - f_trial_minus) / (2. * dx), coordinate_column_j);
				}
			}
			return StatusCode::SUCCESS;
		}
		int Hessian(const boost::numeric::ublas::bounded_vector<T, N> &x, const boost::numeric::ublas::bounded_vector<T, N> &f, T relative_epsilon, boost::numeric::ublas::bounded_vector<T, N> &jacobian, boost::numeric::ublas::bounded_matrix<T, N, N> &hessian) { // TODO: Check!
			namespace ublas = boost::numeric::ublas;
			ublas::bounded_vector<T, N> x1;
			ublas::bounded_vector<T, N> f_trial;
			const T f_norm = ublas::norm_2(f);
			ublas::bounded_vector<T, N> dx;
			// f_{0, 0} = f0
			// f_{-1,0} = f0 - g0 * dx0 + 0.5 * G00 * dx0^2
			// f_{1, 0} = f0 + g0 * dx0 + 0.5 * G00 * dx0^2
			// f_{0,-1} = f0 - g1 * dx1 + 0.5 * G11 * dx1^2
			// f_{0, 1} = f0 + g1 * dx1 + 0.5 * G11 * dx1^2
			// f_{-1,-1}= f0 - g0 * dx0 - g1 * dx1 + 0.5 * G00 * dx0^2 + 0.5 * G11 * dx1^2 + G01 * dx0 * dx1
			// f_{1, 1} = f0 + g0 * dx0 + g1 * dx1 + 0.5 * G00 * dx0^2 + 0.5 * G11 * dx1^2 + G01 * dx0 * dx1
			// f_{-1,1} = f0 - g0 * dx0 + g1 * dx1 + 0.5 * G00 * dx0^2 + 0.5 * G11 * dx1^2 - G01 * dx0 * dx1
			// f_{1,-1} = f0 + g0 * dx0 - g1 * dx1 + 0.5 * G00 * dx0^2 + 0.5 * G11 * dx1^2 - G01 * dx0 * dx1

			x1 = x; // copy x into x1
			for (std::size_t i = 0; i < N; ++i) {
				T x_i = x(i);
				dx(i) = std::abs(x_i) >= 1. ? relative_epsilon * x_i : relative_epsilon * x_i;
				x1(i) = x_i + dx(i);
				if (int f_stat = Eval(x1, f_trial); f_stat != StatusCode::SUCCESS)
					return f_stat;
				T f_trial_norm = ublas::norm_2(f_trial);
				x1(i) = x_i - dx(i);
				if (int f_stat = Eval(x1, f_trial); f_stat != StatusCode::SUCCESS)
					return f_stat;
				x1(i) = x_i;
				jacobian(i) = (f_trial_norm - ublas::norm_2(f_trial)) / (2. * dx(i));
				hessian(i, i) = (f_trial_norm + ublas::norm_2(f_trial) - 2. * f_norm) / power2(dx(i));
			}
			for (std::size_t i = 1; i < N; ++i) {
				x1 = x;
				T x_i = x(i);
				for (std::size_t j = 0; j < i; ++j) {
					T df = 0., x_j = x(j);
					x1(i) = x_i + dx(i);
					x1(j) = x_j + dx(j);
					if (int f_stat = Eval(x1, f_trial); f_stat != StatusCode::SUCCESS) {
						return f_stat;
					}
					df += ublas::norm_2(f_trial);
					x1(j) = x_j - dx(j);
					if (int f_stat = Eval(x1, f_trial); f_stat != StatusCode::SUCCESS) {
						return f_stat;
					}
					df -= ublas::norm_2(f_trial);
					x1(i) = x_i - dx(i);
					if (int f_stat = Eval(x1, f_trial); f_stat != StatusCode::SUCCESS) {
						return f_stat;
					}
					df += ublas::norm_2(f_trial);
					x1(j) = x_j + dx(j);
					if (int f_stat = Eval(x1, f_trial); f_stat != StatusCode::SUCCESS)
						return f_stat;
					df -= ublas::norm_2(f_trial);
					df /= 4. * dx[i] * dx[j];
					hessian(i, j) = hessian(j, i) = df;
				}
			}
			return StatusCode::SUCCESS;
		}
	};

	template <typename T, std::size_t N, typename Parameter>
	class DirectionMultiFunctionSolver : public NonlinearSystemSolver<T, N, Parameter> {
	  protected:
		int directional_num_;
		T trust_radius_;
		T delta_angle_;
		T central_angle_;
		boost::numeric::ublas::bounded_vector<T, N> x_trial_;
		boost::numeric::ublas::bounded_vector<T, N> x_rec_;
		boost::numeric::ublas::bounded_vector<T, N> f_trial_;
		boost::numeric::ublas::bounded_vector<T, N> f_rec_;

	  public:
		DirectionMultiFunctionSolver(int (*function)(const boost::numeric::ublas::bounded_vector<T, N> &, boost::numeric::ublas::bounded_vector<T, N> &, Parameter &), Parameter &params) : NonlinearSystemSolver<T, N, Parameter>(function, params) {}
		int Set(const boost::numeric::ublas::bounded_vector<T, N> &x) {
			namespace ublas = boost::numeric::ublas;
			this->x_ = x;
			directional_num_ = 1;
			delta_angle_ = boost::math::constants::two_thirds_pi<T>();
			central_angle_ = 0.;
			if (int status = this->ScaleX(this->x_, this->f_); status != StatusCode::SUCCESS)
				return status;
			trust_radius_ = std::max(T{1.}, ublas::norm_2(this->x_) / ublas::norm_2(this->f_)) * boost::math::tools::root_epsilon<T>();
			return StatusCode::SUCCESS;
		}
		int Iterate() override {
			namespace ublas = boost::numeric::ublas;
			int status;
			const T f_norm = ublas::norm_2(this->f_), trust_radius_limit = ublas::norm_2(this->x_) * boost::math::tools::epsilon<T>() * 100;
			for (T effective_radius = std::max(trust_radius_limit, trust_radius_ * f_norm); directional_num_ < 64; effective_radius = std::max(trust_radius_limit, trust_radius_ * f_norm)) {
				T min_norm, delta_angle = delta_angle_, next_angle = central_angle_;
				T delta_x = effective_radius * std::cos(central_angle_), delta_y = effective_radius * std::sin(central_angle_);
				x_trial_(0) = this->x_(0) + delta_x;
				x_trial_(1) = this->x_(1) + delta_y;
				if (status = this->ScaleX(x_trial_, f_trial_); status != StatusCode::SUCCESS)
					return status;
				if (min_norm = ublas::norm_2(f_trial_); min_norm < f_norm) {
					next_angle = central_angle_;
					x_rec_ = x_trial_;
					f_rec_ = f_trial_;
				}
				for (int i = directional_num_; i > 0; --i) {
					for (int sign = 1; sign > -2; sign -= 2) {
						T angle = central_angle_ + sign * delta_angle;
						delta_x = trust_radius_ * std::cos(angle);
						delta_y = trust_radius_ * std::sin(angle);
						x_trial_(0) = this->x_(0) + delta_x;
						x_trial_(1) = this->x_(1) + delta_y;
						if (status = this->ScaleX(x_trial_, f_trial_); status != StatusCode::SUCCESS)
							return status;
						if (T f_dir_norm = ublas::norm_2(f_trial_); f_dir_norm < min_norm) {
							min_norm = f_dir_norm;
							next_angle = angle;
							x_rec_ = x_trial_;
							f_rec_ = f_trial_;
						}
					}
					if (delta_angle *= 2.; delta_angle > M_PI) {
						delta_angle *= 0.5;
						break;
					}
				}
				if (min_norm >= f_norm) {
					trust_radius_ *= 0.7;
					if (next_angle != central_angle_ + M_PI)
						delta_angle_ = std::max(delta_angle_, std::abs(next_angle - central_angle_)) * 0.5;
					else
						delta_angle_ = M_PI_4 * std::pow(0.5, directional_num_);
					directional_num_ += 2;
					central_angle_ = ModBy2Pi(next_angle);
					continue;
				}
				delta_x = 0.5 * effective_radius * std::cos(next_angle);
				delta_y = 0.5 * effective_radius * std::sin(next_angle);
				while (effective_radius > trust_radius_limit) {
					x_trial_(0) = this->x_(0) + delta_x;
					x_trial_(1) = this->x_(1) + delta_y;
					if (status = this->ScaleX(x_trial_, f_trial_); status != StatusCode::SUCCESS)
						return status;
					if (T f_mid_norm = ublas::norm_2(f_trial_); f_mid_norm < min_norm) {
						min_norm = f_mid_norm;
						x_rec_ = x_trial_;
						f_rec_ = f_trial_;
						delta_x *= 0.5;
						delta_y *= 0.5;
						trust_radius_ *= 0.5;
					} else
						break;
				}
				this->dx_ = x_rec_ - this->x_;
				this->x_ = x_rec_;
				this->f_ = f_rec_;
				if (directional_num_ > 1) {
					if (next_angle != central_angle_ + M_PI)
						delta_angle_ = std::max(delta_angle_, std::abs(next_angle - central_angle_));
					else
						delta_angle_ = boost::math::constants::half_pi<T>();
					directional_num_ = 1;
				} else {
					if (next_angle == central_angle_)
						delta_angle_ *= 0.5;
					trust_radius_ *= 1.5;
				}
				central_angle_ = ModBy2Pi(next_angle);
				return StatusCode::SUCCESS;
			}
			return StatusCode::MAX_ITERATIONS_EXCEEDED;
		}
	};
} // namespace sbody
#endif
