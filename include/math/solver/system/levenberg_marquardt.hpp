/**
 * @file levenberg_marquardt.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-06
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <numeric>

#include "infra/status.hpp"
#include "math/solver/system/nonlinear_system.hpp"

namespace sbody {
template <typename T, std::size_t N, typename Parameter>
class LevenbergMarquardtSolver : public NonlinearSystemSolver<T, N, Parameter> {
  protected:
    T lambda_;
    T trust_radius_;
    T epsilon_coefficient_;
    T f_norm_square_;
    static constexpr T epsilon_limit_ = boost::math::tools::root_epsilon<T>() * boost::math::tools::forth_root_epsilon<T>();
    bool directional_;
    boost::numeric::ublas::bounded_vector<T, N> x_trial_;
    boost::numeric::ublas::bounded_vector<T, N> f_norm_weight_;
    boost::numeric::ublas::bounded_vector<T, N> f_trial_;
    boost::numeric::ublas::bounded_vector<T, N> gradient_;
    boost::numeric::ublas::bounded_vector<T, N> B_lambda_1_dx;
    boost::numeric::ublas::bounded_matrix<T, N, N> jacobian_;
    boost::numeric::ublas::bounded_matrix<T, N, N> B_;
    boost::numeric::ublas::bounded_matrix<T, N, N> B_lambda_;
    boost::numeric::ublas::permutation_matrix<std::size_t> permutation_;

  public:
    LevenbergMarquardtSolver(int (*function)(const boost::numeric::ublas::bounded_vector<T, N> &,
                                             boost::numeric::ublas::bounded_vector<T, N> &, Parameter &),
                             Parameter &params)
        : NonlinearSystemSolver<T, N, Parameter>(function, params), permutation_(N) {}
    int Set(const boost::numeric::ublas::bounded_vector<T, N> &x, T trust_radius = 0.) {
        namespace ublas = boost::numeric::ublas;
        this->x_ = x;
        if (int status = this->ScaleX(); status != StatusCode::SUCCESS)
            return status;
        for (int i = 0; i < N; ++i)
            f_norm_weight_(i) = 1. / (boost::math::tools::root_epsilon<T>() + std::abs(this->f_(i)));
        f_norm_square_ = ublas::norm_2_square(this->f_); // ublas::norm_2(ublas::element_prod(this->f_, this->f_norm_weight_));
        lambda_ = 0.;
        trust_radius_ = trust_radius <= 0. ? ublas::norm_2(this->x_) : trust_radius;
        epsilon_coefficient_ = boost::math::tools::root_epsilon<T>();
        directional_ = false;
        return StatusCode::SUCCESS;
    }
    int Iterate() override {
        namespace ublas = boost::numeric::ublas;
        const T x_norm = ublas::norm_2(this->x_);
        const T trust_radius_limit = std::max(T{1.}, x_norm) * boost::math::tools::root_epsilon<T>();
        int status;
        while (trust_radius_ >= trust_radius_limit) {
            if (directional_)
                status = this->TwoSidedDirectionalJacobian(this->x_, this->dx_, this->f_, epsilon_coefficient_, jacobian_);
            else
                status = this->TwoSidedJacobian(this->x_, this->f_, epsilon_coefficient_, jacobian_);
            if (status != StatusCode::SUCCESS)
                return status;
            gradient_ = ublas::prod(ublas::trans(jacobian_), this->f_);
            B_ = ublas::prod(ublas::trans(jacobian_), jacobian_);
            std::iota(permutation_.begin(), permutation_.end(), 0);
            if (status = ublas::lu_factorize(jacobian_, permutation_); status != StatusCode::SUCCESS)
                return status;
            this->dx_ = this->f_;
            ublas::lu_substitute(jacobian_, permutation_, this->dx_);
            T dx_norm = ublas::norm_2(this->dx_);
            if (dx_norm <= trust_radius_)
                lambda_ = 0.;
            else {
                std::pair<T, T> lambda_root = boost::math::tools::newton_raphson_iterate(
                    [this, &status](T lambda) -> std::tuple<T, T> {
                        this->B_lambda_ = this->B_ + lambda * ublas::identity_matrix<T>(N);
                        std::iota(this->permutation_.begin(), this->permutation_.end(), 0);
                        if (status = ublas::lu_factorize(this->B_lambda_, this->permutation_); status != StatusCode::SUCCESS)
                            return status;
                        this->dx_ = this->gradient_;
                        ublas::lu_substitute(this->B_lambda_, this->permutation_, this->dx_);
                        this->B_lambda_1_dx = this->dx_;
                        ublas::lu_substitute(this->B_lambda_, this->permutation_, this->B_lambda_1_dx);
                        T dx_norm = ublas::norm_2(this->dx_);
                        return {dx_norm - trust_radius_, ublas::inner_prod(this->dx_, this->B_lambda_1_dx) / dx_norm};
                    },
                    lambda_, T{0.}, ublas::norm_2(gradient_) / trust_radius_, 32);
                lambda_ = 0.5 * (lambda_root.first + lambda_root.second);
            }
            x_trial_ = this->x_ - this->dx_;
            bool scaled = false;
            if (status = this->ScaleX(x_trial_, f_trial_, scaled); status != StatusCode::SUCCESS)
                return status;
            if (scaled) {
                this->dx_ = this->x_ - x_trial_;
                trust_radius_ = std::min(trust_radius_, ublas::norm_2(this->dx_));
            }
            T f_trial_norm_square =
                ublas::norm_2_square(f_trial_); // ublas::norm_2(ublas::element_prod(f_trial_, this->f_norm_weight_));
            T rho = (f_trial_norm_square - f_norm_square_) /
                    (-2. * ublas::inner_prod(gradient_, this->dx_) +
                     ublas::inner_prod(ublas::trans(this->dx_), ublas::prod(B_, this->dx_)));
            if (TrustRegionStep(rho, dx_norm, trust_radius_, x_norm)) {
                // x_trial_ is better than x_
                this->x_ = x_trial_;
                this->f_ = f_trial_;
                f_norm_square_ = f_trial_norm_square;
                return StatusCode::SUCCESS;
            }
            if (!directional_) {
                directional_ = true;
                continue;
            }
        }
        return StatusCode::TRUST_REGION_TOO_SMALL;
    }
    T TrustRadius() const { return trust_radius_; }
};
} // namespace sbody
