/**
 * @file kerr.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-06
 *
 * @copyright Copyright (c) 2026
 *
 */

#ifndef SBODY_KERR_HPP
#define SBODY_KERR_HPP

#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
// #include <boost/math/special_functions/next.hpp>
#include <boost/math/tools/roots.hpp>

#include "math/solver/system/gauss_newton.hpp"
#include "math/special_functions.hpp"
#include "metric/metric.hpp"

namespace sbody {
	template <typename T>
	struct KerrFastTraceParameters;

	template <typename T>
	class Kerr : public Metric<T> {
	  protected:
		const T a_, a2_, a4_;
		const T sqrt_delta_a2_;
		const T u_plus_1, u_plus, u_minus_1, u_minus, u_r;
		static int DeltaUMuPhi(const boost::numeric::ublas::bounded_vector<T, 2> &alpha_beta, boost::numeric::ublas::bounded_vector<T, 2> &delta_u_mu_phi, KerrFastTraceParameters<T> &params) {
			// The photon in the observer's frame has the tetrad velocity: [1, r / R, beta / R, -alpha / R], where R = std::sqrt(r^2 + alpha^2 + beta^2).
			T alpha, beta;
			if (params.line_of_sight_axis) {
				alpha = alpha_beta(0) * std::cos(alpha_beta(1));
				beta = alpha_beta(0) * std::sin(alpha_beta(1));
			} else {
				alpha = alpha_beta(0);
				beta = alpha_beta(1);
			}
			if (!std::isfinite(alpha) || !std::isfinite(beta))
				return StatusCode::NUMERIC_ERROR;
			boost::numeric::ublas::bounded_vector<T, 8> photon;
			T photon_time;
			int status;
			if (status = params.kerr->InitializePhoton(photon, photon_time, alpha, beta, params.r, params.r2, params.theta_obs, params.sin_theta_obs); status != StatusCode::SUCCESS)
				return status;
			const T a = params.kerr->a_, a2 = params.kerr->a2_;
			params.E = params.kerr->Energy(photon, TimeSystem::T, DynamicalSystem::LAGRANGIAN);
			const T E_1 = 1. / params.E;
			params.L = params.kerr->AngularMomentum(photon, TimeSystem::T, DynamicalSystem::LAGRANGIAN);
			const T l = params.L * E_1, l2 = power2(l);
			params.Q = params.kerr->CarterConstant(photon, 0., TimeSystem::T, DynamicalSystem::LAGRANGIAN);
			const T q2 = params.Q * power2(E_1);
			T I_u_0, I_u_1 = 0., I_u_plus_0, I_u_plus_1, I_u_minus_0, I_u_minus_1, I_u_2_0, I_u_2_1, I_u_4_0, I_u_4_1;
			if (UIntegral(a, a2, params.kerr->u_plus_1, params.kerr->u_minus_1, params.kerr->u_plus, params.kerr->u_minus, l, l2, q2, params.r, params.u_obs, params.u_obj, I_u_0, I_u_1, I_u_plus_0, I_u_plus_1, I_u_minus_0, I_u_minus_1, I_u_2_0, I_u_2_1, I_u_4_0, I_u_4_1) == StatusCode::NO_INTERSECTION_WITH_TARGET_SURFACE) {
				// photon failed to hit the r_obj surface, the impact params need to be smaller
				delta_u_mu_phi(0) = std::min(T{0.99999999}, I_u_0 - boost::math::tools::root_epsilon<T>());
				return StatusCode::SCALING_REQUIRED;
			}
			// M = q^2 + (a^2-q^2-l^2)*mu^2-a^2*mu^4
			T M_minus_plus[2];
			T delta_M, delta_M_plus;
			if (int root_num_M = polySolveQuadratic(-a2, a2 - l2 - q2, q2, M_minus_plus); root_num_M == 0)
				return StatusCode::FAILURE;
			delta_M = std::sqrt(power2(a2 - l2 - q2) + 4. * a2 * q2) / a2;
			if (M_minus_plus[1] > 1.) {
				M_minus_plus[1] = 1.;
				delta_M_plus = 0.;
			} else if (M_minus_plus[1] > 0.9) {
				delta_M_plus = 2. * l2 / (a2 + l2 + q2 + std::sqrt(power2(a2 + l2 + q2) - 4. * a2 * l2));
				M_minus_plus[1] = 1. - delta_M_plus;
			} else
				delta_M_plus = 1. - M_minus_plus[1];
			T mu_plus, mu_minus;
			T A, k, n;
			T I_mu_0, I_mu_full_turn;
			T I_t_mu_0, I_t_mu_full_turn;
			T I_phi_mu_0, I_phi_mu_full_turn;
			// Here the first turning point is determined by `beta` of the observer.
			int alpha_1, alpha_2;
			if (status = Mu0Integral(a, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, sign(beta), params.mu_obs, mu_plus, mu_minus, A, k, n, alpha_1, I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn); status != StatusCode::SUCCESS)
				return status;
			T mu_f_0, mu_f_1, t_mu, phi_mu_0, phi_mu_1;
			if (status = MuFIntegral(a, a2, l, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, mu_plus, mu_minus, I_u_0, I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn, A, k, n, alpha_1, alpha_2, mu_f_0, t_mu, phi_mu_0); status != StatusCode::SUCCESS)
				return status;
			const T phi_0 = phi_mu_0 + params.kerr->u_r * ((l * params.kerr->u_plus_1 + 2. * (a - l)) * I_u_plus_0 - (l * params.kerr->u_minus_1 + 2. * (a - l)) * I_u_minus_0);
			params.tau = -t_mu - I_u_4_0;
			params.t = params.tau - 2. * params.kerr->u_r * ((a * (a - l) + power2(params.kerr->u_plus_1)) * I_u_plus_0 - (a * (a - l) + power2(params.kerr->u_minus_1)) * I_u_minus_0 - 2. * std::sqrt(1. - a2) * I_u_2_0);
			params.u_dir = -1.;
			params.mu_dir = sign(mu_plus) * alpha_2;
			// delta_u_mu_phi(0) = std::cos(theta_f_rotate) - params.mu_obj;
			T theta_f, sin_theta_f, cos_theta_f;
			T iota_f, sin_iota_f, cos_iota_f;
			T phi_f, sin_phi_f, cos_phi_f;
			// RotateSphericalAxisWithIota(params.sin_theta_obs, params.mu_obs, params.sin_iota, params.cos_iota,squareRoot(1. - mu_f_0 * mu_f_0), mu_f_0, -std::sin(phi_0), std::cos(phi_0), theta_f, sin_theta_f, cos_theta_f, phi_f, sin_phi_f, cos_phi_f);
			if (params.line_of_sight_axis) {
				rotationAngleOfPerpendicularCoordiante(params.sin_theta_obs, params.mu_obs, squareRoot(1. - mu_f_0 * mu_f_0), mu_f_0, -std::sin(phi_0), std::cos(phi_0), phi_f, sin_phi_f, cos_phi_f, theta_f);
				if (std::abs(phiDifference(alpha_beta(1) - phi_f)) > boost::math::constants::half_pi<T>()) {
					phi_f += boost::math::constants::pi<T>();
					theta_f = boost::math::constants::two_pi<T>() - theta_f;
				}
				delta_u_mu_phi(0) = theta_f - params.theta_obs_obj;
				delta_u_mu_phi(1) = phiDifference(phi_f - params.iota);
			} else {
				rotateSphericalAxisWithIota(params.sin_theta_obs, params.mu_obs, params.sin_iota, params.cos_iota, squareRoot(1. - mu_f_0 * mu_f_0), mu_f_0, -std::sin(phi_0), std::cos(phi_0), theta_f, sin_theta_f, cos_theta_f, phi_f, sin_phi_f, cos_phi_f);
				delta_u_mu_phi(0) = cos_theta_f;
				delta_u_mu_phi(1) = phiDifference(phi_f - params.theta_obs_obj);
			}
			// delta_u_mu_phi(0) = cos_theta_f;
			// if (negative_direction && phi_f > boost::math::constants::pi<T>())
			// 	phi_f -= boost::math::constants::two_pi<T>();
			// delta_u_mu_phi(1) = phi_f - params.theta_observer_object;
			if (I_u_1 == 0.)
				return StatusCode::SUCCESS;
			if (status = MuFIntegral(a, a2, l, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, mu_plus, mu_minus, I_u_1, I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn, A, k, n, alpha_1, alpha_2, mu_f_1, t_mu, phi_mu_1); status != StatusCode::SUCCESS)
				return status;
			const T phi_1 = phi_mu_1 + params.kerr->u_r * ((l * params.kerr->u_plus_1 + 2. * (a - l)) * I_u_plus_1 - (l * params.kerr->u_minus_1 + 2. * (a - l)) * I_u_minus_1);
			const T angle_diff_0 = power2(delta_u_mu_phi(0)) + power2(delta_u_mu_phi(1));
			T delta_theta_1, delta_phi_1;
			// T delta_phi_1 = phi_f - params.theta_observer_object;
			if (params.line_of_sight_axis) {
				rotationAngleOfPerpendicularCoordiante(params.sin_theta_obs, params.mu_obs, squareRoot(1. - mu_f_1 * mu_f_1), mu_f_1, -std::sin(phi_1), std::cos(phi_1), phi_f, sin_phi_f, cos_phi_f, theta_f);
				if (std::abs(phiDifference(alpha_beta(1) - phi_f)) > boost::math::constants::half_pi<T>()) {
					phi_f += boost::math::constants::pi<T>();
					theta_f = boost::math::constants::two_pi<T>() - theta_f;
				}
				delta_theta_1 = theta_f - params.theta_obs_obj;
				delta_phi_1 = phiDifference(phi_f - params.iota);
			} else {
				rotateSphericalAxisWithIota(params.sin_theta_obs, params.mu_obs, params.sin_iota, params.cos_iota, squareRoot(1. - mu_f_1 * mu_f_1), mu_f_1, -std::sin(phi_1), std::cos(phi_1), theta_f, sin_theta_f, cos_theta_f, phi_f, sin_phi_f, cos_phi_f);
				delta_theta_1 = cos_theta_f;
				delta_phi_1 = phiDifference(phi_f - params.theta_obs_obj);
			}
			T angle_diff_1 = power2(delta_theta_1) + power2(delta_phi_1);
			if (angle_diff_0 <= angle_diff_1)
				return StatusCode::SUCCESS;
			params.tau = -t_mu - I_u_4_1;
			params.t = params.tau - 2. * params.kerr->u_r * ((a * (a - l) + power2(params.kerr->u_plus_1)) * I_u_plus_1 - (a * (a - l) + power2(params.kerr->u_minus_1)) * I_u_minus_1 - 2. * std::sqrt(1. - a2) * I_u_2_1);
			params.u_dir = 1.;
			params.mu_dir = sign(mu_plus) * alpha_2;
			delta_u_mu_phi(0) = delta_theta_1;
			delta_u_mu_phi(1) = delta_phi_1;
			return StatusCode::SUCCESS;
		}
		static int UIntegral(T a, T a2, T u_plus_1, T u_minus_1, T u_plus, T u_minus, T l, T l2, T q2, T r_obs, T u_obs, T u_obj, T &I_u_0, T &I_u_1, T &I_u_plus_0, T &I_u_plus_1, T &I_u_minus_0, T &I_u_minus_1, T &I_u_2_0, T &I_u_2_1, T &I_u_4_0, T &I_u_4_1) {
			// U=1+[a^2−q^2−l^2]u^2+2[(a−l)^2+q^2]u^3−a^2q^2u^4
			const T c = a2 - l2 - q2, d = 2. * (power2(a - l) + q2), e = -a2 * q2;
			T u_roots[4];
			if (std::abs(e) < absolute_accuracy) {
				// U=1+[a^2−l^2]u^2+2(a−l)^2u^3
				// This means the observer and the object are on the equatorial plane.
				if (a == l) { // U=1.
					I_u_0 = u_obj - u_obs;
					I_u_plus_0 = u_plus * log((u_obj - u_plus) / (u_obs - u_plus));
					I_u_minus_0 = u_minus * log((u_obj - u_minus) / (u_obs - u_minus));
					I_u_2_0 = log(u_obj * r_obs);
					I_u_4_0 = r_obs - 1. / u_obj;
					return StatusCode::SUCCESS;
				}
				const T d_1 = 1. / d, sqrt_d_1 = boost::math::constants::half_root_two<T>() / std::abs(a - l);
				if (int root_num_U = polySolveCubic(c * d_1, T{0.}, d_1, u_roots); root_num_U == 1) {
					const T f = -1. / (u_roots[0] * d), g = f / u_roots[0];
					I_u_0 = sqrt_d_1 * EllipticIntegral2Complex(0, u_obs, u_obj, T{1.}, T{0.}, f, g, T{1.}, -u_roots[0], T{1.});
					I_u_plus_0 = sqrt_d_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, T{1.}, -u_plus_1, f, g, T{1.}, -u_roots[0], T{1.});
					I_u_minus_0 = sqrt_d_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, T{1.}, -u_minus_1, f, g, T{1.}, -u_roots[0], T{1.});
					I_u_2_0 = sqrt_d_1 * EllipticIntegral2Complex(-2, u_obs, u_obj, T{0.}, T{1.}, f, g, T{1.}, -u_roots[0], T{1.});
					I_u_4_0 = sqrt_d_1 * EllipticIntegral2Complex(-4, u_obs, u_obj, T{0.}, T{1.}, f, g, T{1.}, -u_roots[0], T{1.});
					return StatusCode::SUCCESS;
				}
				if (u_roots[1] < u_obj) {
					I_u_0 = u_roots[1] / u_obj;
					return StatusCode::NO_INTERSECTION_WITH_TARGET_SURFACE;
				}
				if (u_roots[1] >= u_plus) { // photon falls into the BH
					I_u_0 = sqrt_d_1 * EllipticIntegral(0, u_obs, u_obj, T{1.}, T{0.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.});
					I_u_plus_0 = sqrt_d_1 * -EllipticIntegral(-2, u_obs, u_obj, T{1.}, -u_plus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.});
					I_u_minus_0 = sqrt_d_1 * -EllipticIntegral(-2, u_obs, u_obj, T{1.}, -u_minus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.});
					I_u_2_0 = sqrt_d_1 * EllipticIntegral(-2, u_obs, u_obj, T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.});
					I_u_4_0 = sqrt_d_1 * EllipticIntegral(-4, u_obs, u_obj, T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.});
					return StatusCode::SUCCESS;
				}
				aMinusPlusB(sqrt_d_1 * EllipticIntegral(0, u_obs, u_roots[1], T{1.}, T{0.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}), sqrt_d_1 * EllipticIntegral(0, u_obj, u_roots[1], T{1.}, T{0.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}), I_u_0, I_u_1);
				aMinusPlusB(sqrt_d_1 * -EllipticIntegral(-2, u_obs, u_roots[1], T{1.}, -u_plus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}), sqrt_d_1 * -EllipticIntegral(-2, u_obj, u_roots[1], T{1.}, -u_plus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}), I_u_plus_0, I_u_plus_1);
				aMinusPlusB(sqrt_d_1 * -EllipticIntegral(-2, u_obs, u_roots[1], T{1.}, -u_minus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}), sqrt_d_1 * -EllipticIntegral(-2, u_obj, u_roots[1], T{1.}, -u_minus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}), I_u_minus_0, I_u_minus_1);
				aMinusPlusB(sqrt_d_1 * EllipticIntegral(-2, u_obs, u_roots[1], T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}), sqrt_d_1 * EllipticIntegral(-2, u_obj, u_roots[1], T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}), I_u_2_0, I_u_2_1);
				aMinusPlusB(sqrt_d_1 * EllipticIntegral(-4, u_obs, u_roots[1], T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}), sqrt_d_1 * EllipticIntegral(-4, u_obj, u_roots[1], T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}), I_u_4_0, I_u_4_1);
				return StatusCode::SUCCESS;
			}
			const T e_1 = 1. / e, sqrt_e = std::sqrt(std::abs(e)), sqrt_e_1 = 1. / sqrt_e;
			const T c_sqrt_e = c * sqrt_e_1, d_e = d * e_1;
			if (int root_num_U = polySolveQuartic(d_e, c * e_1, T{0.}, e_1, u_roots); root_num_U == 0) { // q2 < 0., e > 0.
				const T h1 = boost::math::tools::halley_iterate(
					// Equation (21) in Dexter & Agol (2009) ApJ 696
					[c_sqrt_e, cubic_coeff = 2. * c_sqrt_e - sqrt_e * power2(d_e)](T h1) -> std::tuple<T, T, T> {
						const T h12 = power2(h1), h14 = power2(h12);
						return {(h12 - 1.) * (h14 - 1.) - h1 * (c_sqrt_e * (h14 + 1.) - cubic_coeff * h12),
								h1 * (6. * h14 - 4. * h12 - 2.) - c_sqrt_e * (5. * h14 + 1.) + 3. * cubic_coeff * h12,
								30. * h14 - (c_sqrt_e * 20. * h1 + 12.) * h12 + 6. * cubic_coeff * h1 - 2.};
					},
					// f(0) = 1, f(1) = -sqrt(e) * (d/e)^2 < 0
					T{0.5}, T{0.0}, T{1.0}, std::numeric_limits<T>::digits);
				const T h2 = 1. / h1, g1 = d_e / (h2 - h1);
				I_u_0 = sqrt_e_1 * EllipticIntegral4Complex(0, u_obs, u_obj, T{1.}, T{0.}, sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				I_u_plus_0 = sqrt_e_1 * -EllipticIntegral4Complex(-2, u_obs, u_obj, T{1.}, -u_plus_1, sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				I_u_minus_0 = sqrt_e_1 * -EllipticIntegral4Complex(-2, u_obs, u_obj, T{1.}, -u_minus_1, sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				I_u_2_0 = sqrt_e_1 * EllipticIntegral4Complex(-2, u_obs, u_obj, T{0.}, T{1.}, sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				I_u_4_0 = sqrt_e_1 * EllipticIntegral4Complex(-4, u_obs, u_obj, T{0.}, T{1.}, sqrt_e_1, g1, h1, sqrt_e_1, -g1, h2);
				if (std::isnan(I_u_0) || std::isnan(I_u_plus_0) || std::isnan(I_u_minus_0) || std::isnan(I_u_2_0) || std::isnan(I_u_4_0))
					return StatusCode::NUMERIC_ERROR;
				return StatusCode::SUCCESS;
			} else if (root_num_U == 2) {
				const T f = e_1 / (u_roots[0] * u_roots[1]), g = (1. / u_roots[0] + 1. / u_roots[1]) * f; // 1. / (−a^2*q^2*u_1*u_4)
				if (q2 < 0.) {
					// root[0] < root[1] < 0.
					I_u_0 = sqrt_e_1 * EllipticIntegral2Complex(0, u_obs, u_obj, T{1.}, T{0.}, f, g, T{1.}, -u_roots[0], T{1.}, -u_roots[1], T{1.});
					I_u_plus_0 = sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, T{1.}, -u_plus_1, f, g, T{1.}, -u_roots[0], T{1.}, -u_roots[1], T{1.});
					I_u_minus_0 = sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, T{1.}, -u_minus_1, f, g, T{1.}, -u_roots[0], T{1.}, -u_roots[1], T{1.});
					I_u_2_0 = sqrt_e_1 * EllipticIntegral2Complex(-2, u_obs, u_obj, T{0.}, T{1.}, f, g, T{1.}, -u_roots[0], T{1.}, -u_roots[1], T{1.});
					I_u_4_0 = sqrt_e_1 * EllipticIntegral2Complex(-4, u_obs, u_obj, T{0.}, T{1.}, f, g, T{1.}, -u_roots[0], T{1.}, -u_roots[1], T{1.});
					if (std::isnan(I_u_0) || std::isnan(I_u_plus_0) || std::isnan(I_u_minus_0) || std::isnan(I_u_2_0) || std::isnan(I_u_4_0))
						return StatusCode::NUMERIC_ERROR;
					return StatusCode::SUCCESS;
				}
				if (u_roots[1] < u_obj) {
					I_u_0 = u_roots[1] / u_obj;
					return StatusCode::NO_INTERSECTION_WITH_TARGET_SURFACE;
				}
				if (u_roots[1] >= u_plus) { // photon falls into the BH
					I_u_0 = sqrt_e_1 * EllipticIntegral2Complex(0, u_obs, u_obj, T{1.}, T{0.}, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.});
					I_u_plus_0 = sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, T{1.}, -u_plus_1, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.});
					I_u_minus_0 = sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_obj, T{1.}, -u_minus_1, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.});
					I_u_2_0 = sqrt_e_1 * EllipticIntegral2Complex(-2, u_obs, u_obj, T{0.}, T{1.}, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.});
					I_u_4_0 = sqrt_e_1 * EllipticIntegral2Complex(-4, u_obs, u_obj, T{0.}, T{1.}, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.});
					if (std::isnan(I_u_0) || std::isnan(I_u_plus_0) || std::isnan(I_u_minus_0) || std::isnan(I_u_2_0) || std::isnan(I_u_4_0))
						return StatusCode::NUMERIC_ERROR;
					return StatusCode::SUCCESS;
				}
				aMinusPlusB(sqrt_e_1 * EllipticIntegral2Complex(0, u_obs, u_roots[1], T{1.}, T{0.}, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}), sqrt_e_1 * EllipticIntegral2Complex(0, u_obj, u_roots[1], T{1.}, T{0.}, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}), I_u_0, I_u_1);
				aMinusPlusB(sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_roots[1], T{1.}, -u_plus_1, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}), sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obj, u_roots[1], T{1.}, -u_plus_1, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}), I_u_plus_0, I_u_plus_1);
				aMinusPlusB(sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obs, u_roots[1], T{1.}, -u_minus_1, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}), sqrt_e_1 * -EllipticIntegral2Complex(-2, u_obj, u_roots[1], T{1.}, -u_minus_1, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}), I_u_minus_0, I_u_minus_1);
				aMinusPlusB(sqrt_e_1 * EllipticIntegral2Complex(-2, u_obs, u_roots[1], T{0.}, T{1.}, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}), sqrt_e_1 * EllipticIntegral2Complex(-2, u_obj, u_roots[1], T{0.}, T{1.}, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}), I_u_2_0, I_u_2_1);
				aMinusPlusB(sqrt_e_1 * EllipticIntegral2Complex(-4, u_obs, u_roots[1], T{0.}, T{1.}, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}), sqrt_e_1 * EllipticIntegral2Complex(-4, u_obj, u_roots[1], T{0.}, T{1.}, f, g, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}), I_u_4_0, I_u_4_1);
				return StatusCode::SUCCESS;
			}
			// root_num_U == 4
			if (u_roots[1] < u_obj) {
				I_u_0 = u_roots[1] / u_obj;
				return StatusCode::NO_INTERSECTION_WITH_TARGET_SURFACE;
			}
			if (u_roots[1] >= 1.) { // photon falls into the BH
				I_u_0 = sqrt_e_1 * EllipticIntegral(0, u_obs, u_obj, T{1.}, T{0.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.});
				I_u_plus_0 = sqrt_e_1 * -EllipticIntegral(-2, u_obs, u_obj, T{1.}, -u_plus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.});
				I_u_minus_0 = sqrt_e_1 * -EllipticIntegral(-2, u_obs, u_obj, T{1.}, -u_minus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.});
				I_u_2_0 = sqrt_e_1 * EllipticIntegral(-2, u_obs, u_obj, T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.});
				I_u_4_0 = sqrt_e_1 * EllipticIntegral(-4, u_obs, u_obj, T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.});
				return StatusCode::SUCCESS;
			}
			aMinusPlusB(sqrt_e_1 * EllipticIntegral(0, u_obs, u_roots[1], T{1.}, T{0.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.}), sqrt_e_1 * EllipticIntegral(0, u_obj, u_roots[1], T{1.}, T{0.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.}), I_u_0, I_u_1);
			aMinusPlusB(sqrt_e_1 * -EllipticIntegral(-2, u_obs, u_roots[1], T{1.}, -u_plus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.}), sqrt_e_1 * -EllipticIntegral(-2, u_obj, u_roots[1], T{1.}, -u_plus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.}), I_u_plus_0, I_u_plus_1);
			aMinusPlusB(sqrt_e_1 * -EllipticIntegral(-2, u_obs, u_roots[1], T{1.}, -u_minus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.}), sqrt_e_1 * -EllipticIntegral(-2, u_obj, u_roots[1], T{1.}, -u_minus_1, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.}), I_u_minus_0, I_u_minus_1);
			aMinusPlusB(sqrt_e_1 * EllipticIntegral(-2, u_obs, u_roots[1], T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.}), sqrt_e_1 * EllipticIntegral(-2, u_obj, u_roots[1], T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.}), I_u_2_0, I_u_2_1);
			aMinusPlusB(sqrt_e_1 * EllipticIntegral(-4, u_obs, u_roots[1], T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.}), sqrt_e_1 * EllipticIntegral(-4, u_obj, u_roots[1], T{0.}, T{1.}, u_roots[1], T{-1.}, -u_roots[0], T{1.}, u_roots[2], T{-1.}, u_roots[3], T{-1.}), I_u_4_0, I_u_4_1);
			return StatusCode::SUCCESS;
		}
		static int Mu0Integral(T a, T q2, T M_plus, T M_minus, T delta_M, T delta_M_plus, int beta_sign, T mu_obs, T &mu_plus, T &mu_minus, T &A, T &k, T &n, int &alpha_1, T &I_mu_0, T &I_mu_full_turn, T &I_t_mu_0, T &I_t_mu_full_turn, T &I_phi_mu_0, T &I_phi_mu_full_turn) {
			if (std::abs(q2) < absolute_accuracy) { // M_minus == 0.
				mu_plus = std::copysign(squareRoot(M_plus), mu_obs);
				mu_minus = 0.;
				alpha_1 = beta_sign * sign(mu_obs);
				A = mu_plus * std::abs(a);
				k = 1.;
				const T mu0_mu_plus = std::min(T{1.}, std::abs(mu_obs) / mu_plus);
				I_mu_full_turn = std::numeric_limits<T>::infinity();
				I_t_mu_full_turn = std::numeric_limits<T>::max();
				const T phi = std::acos(mu0_mu_plus);
				I_t_mu_0 = squareRoot(1. - power2(mu0_mu_plus));
				I_mu_0 = std::log((1. + I_t_mu_0) / mu0_mu_plus) / A;
				if (M_plus == 1.)
					return StatusCode::SUCCESS;
				n = M_plus / delta_M_plus;
				I_phi_mu_full_turn = std::numeric_limits<T>::max();
				I_phi_mu_0 = boost::math::ellint_3(k, -n, phi); // TODO: Simplify
				return StatusCode::SUCCESS;
			}
			if (M_minus > 0.) {
				mu_plus = std::copysign(squareRoot(M_plus), mu_obs);
				mu_minus = std::copysign(squareRoot(M_minus), mu_obs);
				alpha_1 = beta_sign * sign(mu_obs);
				A = mu_plus * std::abs(a);
				k = squareRoot(delta_M / M_plus);
				const T x = std::min(T{1.}, squareRoot((M_plus - power2(mu_obs)) / delta_M));
				I_mu_full_turn = boost::math::ellint_1(k) / A;
				I_t_mu_full_turn = boost::math::ellint_2(k);
				const T phi = asin(x);
				I_mu_0 = boost::math::ellint_1(k, phi) / A;
				I_t_mu_0 = boost::math::ellint_2(k, phi);
				if (M_plus == 1.)
					return StatusCode::SUCCESS;
				n = delta_M / delta_M_plus;
				I_phi_mu_full_turn = boost::math::ellint_3(k, -n);
				I_phi_mu_0 = boost::math::ellint_3(k, -n, phi);
				return StatusCode::SUCCESS;
			}
			// M_minus < 0.
			mu_plus = squareRoot(M_plus);
			mu_minus = -mu_plus;
			alpha_1 = beta_sign;
			A = squareRoot(delta_M) * std::abs(a);
			k = squareRoot(M_plus / delta_M);
			// const T x =squareRoot(1. -power2(mu_obs / mu_plus));
			I_mu_full_turn = 2. * boost::math::ellint_1(k) / A;
			I_t_mu_full_turn = 2. * boost::math::ellint_2(k);
			const T phi = acos(std::min(T{1.}, std::abs(mu_obs) / mu_plus));
			if (mu_obs >= 0.) {
				I_mu_0 = boost::math::ellint_1(k, phi) / A;
				I_t_mu_0 = boost::math::ellint_2(k, phi);
				if (M_plus == 1.)
					return StatusCode::SUCCESS;
				n = M_plus / delta_M_plus;
				I_phi_mu_full_turn = 2. * boost::math::ellint_3(k, -n);
				I_phi_mu_0 = boost::math::ellint_3(k, -n, phi);
				return StatusCode::SUCCESS;
			}
			// the photon has to cross the equatorial plane.
			I_mu_0 = I_mu_full_turn - boost::math::ellint_1(k, phi) / A;
			I_t_mu_0 = I_t_mu_full_turn - boost::math::ellint_2(k, phi);
			if (M_plus == 1.)
				return StatusCode::SUCCESS;
			n = M_plus / delta_M_plus;
			I_phi_mu_full_turn = 2. * boost::math::ellint_3(k, -n);
			I_phi_mu_0 = I_phi_mu_full_turn - boost::math::ellint_3(k, -n, phi);
			return StatusCode::SUCCESS;
		}
		static int MuFIntegral(T a, T a2, T l, T q2, T M_plus, T M_minus, T delta_M, T delta_M_plus, T mu_plus, T mu_minus, T I_u, T I_mu_0, T I_mu_full_turn, T I_t_mu_0, T I_t_mu_full_turn, T I_phi_mu_0, T I_phi_mu_full_turn, T A, T k, T n, int alpha_1, int &alpha_2, T &mu_f, T &t_mu, T &phi_mu) {
			const int N = alpha_1 > 0 ? ceil((I_u - I_mu_0) / I_mu_full_turn) : floor((I_u + I_mu_0) / I_mu_full_turn); // s_mu = alpha_1
			int alpha_3;
			// Here we follow the calculation in the code GEOKERR.
			// alpha_2 = (N & 1) ? alpha_1 : -alpha_1;
			// alpha_3 = 2 * floor((2 * N + 1 - alpha_1) / 4);
			if (N & 1) {
				alpha_2 = alpha_1;
				if (alpha_1 == 1)
					alpha_3 = N - 1;
				else
					alpha_3 = N + 1;
			} else {
				alpha_2 = -alpha_1;
				alpha_3 = N;
			}
			T phi;
			if (std::abs(q2) < absolute_accuracy) {
				// assert alpha_3 == 0
				const T cos_phi = 1. / cosh(std::abs(a * mu_plus) * (I_u - alpha_1 * I_mu_0));
				const T sin_phi = tanh(std::abs(a * mu_plus) * (I_u - alpha_1 * I_mu_0));
				mu_f = mu_plus * cos_phi;
				if (cos_phi < sin_phi)
					phi = acos(cos_phi);
				else
					phi = asin(sin_phi);
				t_mu = A * (alpha_1 * I_t_mu_0 + alpha_2 * sin_phi + alpha_3 * I_t_mu_full_turn);
				if (M_plus == 1.) {
					if (alpha_1 != alpha_2)
						phi_mu = 0.;
					else
						phi_mu = std::copysign(boost::math::constants::pi<T>(), l);
				} else
					phi_mu = l * (-I_u + (alpha_1 * I_phi_mu_0 + alpha_2 * boost::math::ellint_3(k, -n, phi) + alpha_3 * I_phi_mu_full_turn) / (A * delta_M_plus)); // TODO: simplify;
				return StatusCode::SUCCESS;
			}
			const T I = (I_u - alpha_1 * I_mu_0 - alpha_3 * I_mu_full_turn) / alpha_2;
			if (M_minus > 0.) {
				mu_f = mu_plus * boost::math::jacobi_dn(k, mu_plus * std::abs(a) * I);
				if (T sin2_phi = (M_plus - power2(mu_f)) / delta_M; sin2_phi < 0.5)
					phi = asin(squareRoot(sin2_phi));
				else
					phi = acos(squareRoot((power2(mu_f) - M_minus) / delta_M));
				t_mu = A * (alpha_1 * I_t_mu_0 + alpha_2 * boost::math::ellint_2(k, phi) + alpha_3 * I_t_mu_full_turn);
				if (M_plus == 1.) {
					if (alpha_1 != alpha_2)
						phi_mu = 0.;
					else
						phi_mu = std::copysign(boost::math::constants::pi<T>(), l);
				} else
					phi_mu = l * (-I_u + (alpha_1 * I_phi_mu_0 + alpha_2 * boost::math::ellint_3(k, -n, phi) + alpha_3 * I_phi_mu_full_turn) / (A * delta_M_plus));
				return StatusCode::SUCCESS;
			}
			T J;
			T sn, cn, dn;
			// M_minus < 0.
			if (I <= 0.5 * I_mu_full_turn) {
				sn = boost::math::jacobi_elliptic(k, squareRoot(delta_M) * std::abs(a) * I, &cn, &dn);
				mu_f = mu_plus * cn;
				if (cn < sn)
					phi = std::acos(cn);
				else
					phi = std::asin(sn);
				t_mu = a2 * M_minus * I_u + A * (alpha_1 * I_t_mu_0 + alpha_2 * boost::math::ellint_2(k, phi) + alpha_3 * I_t_mu_full_turn);
				if (M_plus == 1.) {
					if (alpha_1 != alpha_2)
						phi_mu = 0.;
					else
						phi_mu = std::copysign(boost::math::constants::pi<T>(), l);
				} else
					phi_mu = l * (-I_u + (alpha_1 * I_phi_mu_0 + alpha_2 * boost::math::ellint_3(k, -n, phi) + alpha_3 * I_phi_mu_full_turn) / (A * delta_M_plus));
				return StatusCode::SUCCESS;
			}
			// I > 0.5 * I_mu_full_turn, the photon has to cross the equatorial plane.
			sn = boost::math::jacobi_elliptic(k, squareRoot(delta_M) * std::abs(a) * (I_mu_full_turn - I), &cn, &dn);
			mu_f = mu_minus * cn;
			if (cn < sn)
				phi = std::acos(cn);
			else
				phi = std::asin(sn);
			t_mu = a2 * M_minus * I_u + A * (alpha_1 * I_t_mu_0 - alpha_2 * boost::math::ellint_2(k, phi) + (alpha_2 + alpha_3) * I_t_mu_full_turn);
			if (M_plus == 1.) {
				if (alpha_1 != alpha_2)
					phi_mu = 0.;
				else
					phi_mu = std::copysign(boost::math::constants::pi<T>(), l);
			} else
				phi_mu = l * (-I_u + (alpha_1 * I_phi_mu_0 - alpha_2 * boost::math::ellint_3(k, -n, phi) + (alpha_2 + alpha_3) * I_phi_mu_full_turn) / (A * delta_M_plus));
			return StatusCode::SUCCESS;
		}
		int DerivePhotonInformationFromSolverRoot() {
		}

	  public:
		Kerr(T spin) : a_(spin), a2_(a_ * a_), a4_(a2_ * a2_), sqrt_delta_a2_(std::sqrt(1. - a2_)), u_plus_1(1. + sqrt_delta_a2_), u_plus(1. / u_plus_1), u_minus_1(a2_ * u_plus), u_minus(1. / u_minus_1), u_r(-0.5 / sqrt_delta_a2_) {}
		std::string Name() const override {
			return "Kerr";
		}
		int MetricTensor(const boost::numeric::ublas::bounded_vector<T, 8> &position, boost::numeric::ublas::bounded_matrix<T, 4, 4> &metric) const override {
			const T r = position[1], r2 = r * r, Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * power2(std::cos(position[2])), sin2_theta = power2(std::sin(position[2]));
			if (rho2 == 0. || Delta == 0.)
				return StatusCode::DIVISION_BY_ZERO;
			const T r_rho_2 = 2. * r / rho2;
			metric.clear();
			metric(0, 0) = -(1. - r_rho_2);
			metric(0, 3) = -r_rho_2 * a_ * sin2_theta;
			metric(1, 1) = rho2 / Delta;
			metric(2, 2) = rho2;
			metric(3, 0) = -r_rho_2 * a_ * sin2_theta;
			metric(3, 3) = (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * sin2_theta;
			return StatusCode::SUCCESS;
		}
		T DotProduct3D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x1, T x2, T x3, T y1, T y2, T y3) const override {
			const T r = position[1], r2 = power2(r), Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * power2(std::cos(position[2])), r_rho_2 = 2. * r / rho2, sin2_theta = power2(std::sin(position[2]));
			return (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * sin2_theta * sin2_theta) * x3 * y3;
		}
		T DotProduct4D(const boost::numeric::ublas::bounded_vector<T, 8> &position, T x0, T x1, T x2, T x3, T y0, T y1, T y2, T y3) const override {
			const T r = position[1], r2 = power2(r), Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * power2(std::cos(position[2])), r_rho_2 = 2. * r / rho2, sin2_theta = power2(std::sin(position[2]));
			return (r_rho_2 - 1.) * x0 * y0 - r_rho_2 * a_ * sin2_theta * (x0 * y3 + x3 * y0) + (rho2 / Delta) * x1 * y1 + rho2 * x2 * y2 + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * sin2_theta * sin2_theta) * x3 * y3;
		}
		T DistanceSquare(const boost::numeric::ublas::bounded_vector<T, 8> &x, const boost::numeric::ublas::bounded_vector<T, 8> &y, bool include_time) const override {
			const T r = x[1], r2 = r * r, d0 = x[0] - y[0], d3 = phiDifference(x[3] - y[3]);
			const T Delta = r2 - 2. * r + a2_, rho2 = r2 + a2_ * power2(std::cos(x[2])), r_rho_2 = 2. * r / rho2, sin2_theta = power2(std::sin(x[2]));
			if (include_time)
				return (r_rho_2 - 1.) * power2(d0) - 2. * r_rho_2 * a_ * sin2_theta * d0 * d3 + (rho2 / Delta) * power2(x[1] - y[1]) + rho2 * power2(x[2] - y[2]) + ((r2 + a2_) * sin2_theta + r_rho_2 * a2_ * power2(sin2_theta)) * power2(d3);
			return (rho2 / Delta) * power2(x[1] - y[1]) + rho2 * power2(x[2] - y[2]) + ((r2 + a2_) * sin2_theta + 2. * r * a2_ * power2(sin2_theta) / rho2) * power2(d3);
		}
		int LagrangianToHamiltonian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
			const T dt_dtau = 1. / y[4], r2 = power2(y[1]), sin2_theta = power2(std::sin(y[2]));
			const T Delta = r2 - 2. * y[1] + a2_, rho2 = r2 + a2_ * power2(std::cos(y[2]));
			if (rho2 == 0. || Delta == 0.)
				return StatusCode::DIVISION_BY_ZERO;
			const T r_rho_2 = 2. * y[1] / rho2;
			y[4] = (y[4] - 1. + r_rho_2 * (1. - a_ * sin2_theta * y[7])) * dt_dtau; // 1 + p_t
			y[5] *= rho2 / Delta * dt_dtau;
			y[6] *= rho2 * dt_dtau;
			y[7] = (-r_rho_2 * a_ + (r2 + a2_ * (1. + r_rho_2 * sin2_theta)) * y[7]) * sin2_theta * dt_dtau;
			return StatusCode::SUCCESS;
		}
		int HamiltonianToLagrangian(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
			const T pt = y[4] - 1., r2 = power2(y[1]), sin2_theta = power2(std::sin(y[2]));
			const T Delta = r2 - 2. * y[1] + a2_, rho2 = r2 + a2_ * power2(std::cos(y[2]));
			if (rho2 == 0. || Delta == 0.)
				return StatusCode::DIVISION_BY_ZERO;
			const T rho_2 = 1. / rho2, r_rho_2 = 2. * y[1] * rho_2;
			y[4] = -Delta / ((Delta + r_rho_2 * (a2_ + r2)) * pt + r_rho_2 * a_ * y[7]);
			y[5] *= Delta * y[4] * rho_2;
			y[6] *= y[4] * rho_2;
			if (y[7] == 0.)
				y[7] = -r_rho_2 * a_ * pt / Delta * y[4];
			else if (sin2_theta == 0.)
				return StatusCode::DIVISION_BY_ZERO;
			else
				y[7] = (-r_rho_2 * a_ * pt + (1. - r_rho_2) / sin2_theta * y[7]) / Delta * y[4];
			return StatusCode::SUCCESS;
		}
		int FastTrace(const T r_observer, const T theta_observer, const T sin_theta_observer, const T cos_theta_observer, const T r_object, const T theta_object, const T phi_object, T &alpha, T &beta, boost::numeric::ublas::bounded_vector<T, 8> &photon, T &photon_time) const override {
			namespace ublas = boost::numeric::ublas;
			const T sin_theta_object = sinTheta(theta_object), cos_theta_object = cosTheta(theta_object), sin_phi_object = std::sin(phi_object), cos_phi_object = std::cos(phi_object);
			const T sin2_theta_object = power2(sin_theta_object), cos2_theta_object = power2(cos_theta_object);
			const T cos_observer_object = sin_theta_observer * sin_theta_object * cos_phi_object + cos_theta_observer * cos_theta_object;
			const T theta_observer_object = acos(cos_observer_object), sin_observer_object = squareRoot(1. - power2(cos_observer_object));
			const T r2_object = power2(r_object), r2_observer = power2(r_observer), u_observer = 1. / r_observer, u_object = 1. / r_object;
			const T Delta_object = r_object * (r_object - 2.) + a2_, rho2_object = r2_object + a2_ * cos2_theta_object;
			const T rho_2_object = 1. / rho2_object;
			const T r_rho_2_object = 2. * r_object * rho_2_object;
			const T g00_object = -(1. - r_rho_2_object), g03_object = -r_rho_2_object * a_ * sin2_theta_object;
			const T g11_object = rho2_object / Delta_object, g22_object = rho2_object, g33_object = (r2_object + a2_) * sin2_theta_object - g03_object * a_ * sin2_theta_object;
			T iota, sin_iota, cos_iota, phi1;
			int status;
			rotationAngleOfPerpendicularCoordiante(sin_theta_observer, cos_theta_observer, sin_theta_object, cos_theta_object, sin_phi_object, cos_phi_object, iota, sin_iota, cos_iota, phi1);
			KerrFastTraceParameters<T> fast_trace_parameters(this, r_observer, r2_observer, u_observer, u_object, cos_theta_observer, cos_theta_object, theta_observer, sin_theta_observer, iota, sin_iota, cos_iota, theta_observer_object, false);
			ublas::bounded_vector<T, 2> alpha_beta_initial_value;
			const T effective_radius = r_object * sin_observer_object + power3(theta_observer_object / boost::math::constants::half_pi<T>());
			// const T alpha_coefficient = sin_theta_object * sin_phi_object, beta_coefficient = cos_theta_object * sin_theta_observer - sin_theta_object * cos_phi_object * cos_theta_observer;
			alpha_beta_initial_value(0) = effective_radius * cos_iota;
			alpha_beta_initial_value(1) = effective_radius * sin_iota;
			GaussNewtonSolver<T, 2, KerrFastTraceParameters<T>> alpha_beta_solver(DeltaUMuPhi, fast_trace_parameters);
			if ((status = alpha_beta_solver.Set(alpha_beta_initial_value)) != StatusCode::SUCCESS || (status = alpha_beta_solver.Solve(boost::math::tools::root_epsilon<double>(), 512)) != StatusCode::SUCCESS) {
				fast_trace_parameters.line_of_sight_axis = true;
				alpha_beta_initial_value(0) = effective_radius;
				alpha_beta_initial_value(1) = iota;
				GaussNewtonSolver<T, 2, KerrFastTraceParameters<T>> alpha_beta_line_of_sight_solver(DeltaUMuPhi, fast_trace_parameters);
				if ((status = alpha_beta_line_of_sight_solver.Set(alpha_beta_initial_value)) != StatusCode::SUCCESS || (status = alpha_beta_line_of_sight_solver.Solve(boost::math::tools::root_epsilon<double>(), 512)) != StatusCode::SUCCESS)
					return status;
				// DNewtonRotationMultiFunctionSolver<T, 2, KerrFastTraceParameters<T>> alpha_beta_rotation_solver(DeltaUMuPhi, fast_trace_parameters);
				// if (alpha_beta_rotation_solver.Set(alpha_beta_solver.Root(), theta_observer, sin_theta_observer, cos_theta_observer, r_object, sin_theta_object, cos_theta_object, phi_object, sin_phi_object, cos_phi_object, SPIN_CENTER_OFFSET * a_ * sin_theta_object / r_object, 0., false, alpha_beta_solver.TrustRadius()) != StatusCode::SUCCESS || alpha_beta_rotation_solver.Solve(boost::math::tools::root_epsilon<double>()) != StatusCode::SUCCESS) {
				// 	return StatusCode::FAILURE;
				// 	// 	// PrintlnWarning("Kerr FastTrace() ROTATION failed with status = {}", status);
				// 	// 	DNewtonTranslationMultiFunctionSolver<T, 2, KerrFastTraceParameters<T>> alpha_beta_translation_solver(DeltaUMuPhi, fast_trace_parameters);
				// 	// 	if (alpha_beta_translation_solver.Set(alpha_beta_rotation_solver.Root(), theta_observer, sin_theta_observer, cos_theta_observer, r_object, sin_theta_object, cos_theta_object, phi_object, sin_phi_object, cos_phi_object, false) != StatusCode::SUCCESS || alpha_beta_translation_solver.Solve(boost::math::tools::root_epsilon<double>()) != StatusCode::SUCCESS) {
				// 	// 		DirectionMultiFunctionSolver<T, 2, KerrFastTraceParameters<T>> alpha_beta_direction_solver(DeltaUMuPhi, fast_trace_parameters);
				// 	// 		if (status = alpha_beta_direction_solver.Set(alpha_beta_translation_solver.Root()); status != StatusCode::SUCCESS)
				// 	// 			return status;
				// 	// 		// PrintlnWarning("Kerr FastTrace() set DIRECTION failed with status = {}", status);
				// 	// 		if (status = alpha_beta_direction_solver.Solve(boost::math::tools::root_epsilon<double>(), 2048); status != StatusCode::SUCCESS)
				// 	// 			return status;
				// 	// 		alpha = alpha_beta_direction_solver.Root()(0);
				// 	// 		beta = alpha_beta_direction_solver.Root()(1);
				// 	// 	} else {
				// 	// 		alpha = alpha_beta_translation_solver.Root()(0);
				// 	// 		beta = alpha_beta_translation_solver.Root()(1);
				// 	// 	}
				// } else {
				// 	alpha = alpha_beta_rotation_solver.Root()(0) * std::cos(alpha_beta_rotation_solver.Root()(1));
				// 	beta = alpha_beta_rotation_solver.Root()(0) * std::sin(alpha_beta_rotation_solver.Root()(1));
				// }
				else {
					alpha = alpha_beta_line_of_sight_solver.Root()(0) * std::cos(alpha_beta_line_of_sight_solver.Root()(1));
					beta = alpha_beta_line_of_sight_solver.Root()(0) * std::sin(alpha_beta_line_of_sight_solver.Root()(1));
				}
			} else {
				alpha = alpha_beta_solver.Root()(0);
				beta = alpha_beta_solver.Root()(1);
			}
			photon(0) = fast_trace_parameters.tau; // not important
			photon(1) = r_object;
			photon(2) = theta_object;
			photon(3) = phi_object;
			photon(4) = (power2(g03_object) - g00_object * g33_object) / (fast_trace_parameters.E * g33_object + fast_trace_parameters.L * g03_object);
			photon(7) = -(g00_object * fast_trace_parameters.L + g03_object * fast_trace_parameters.E) / (g03_object * fast_trace_parameters.L + g33_object * fast_trace_parameters.E);
			photon(6) = -fast_trace_parameters.mu_dir * std::sqrt(fast_trace_parameters.Q * power2(photon(4)) - cos2_theta_object * (power2(-r_rho_2_object * a_ + (a2_ + r2_object + r_rho_2_object * a2_ * sin2_theta_object) * photon(7)) * sin2_theta_object - a2_ * power2(r_rho_2_object * (1. - a_ * sin2_theta_object * photon(7)) - 1.))) * rho_2_object;
			photon(5) = -fast_trace_parameters.u_dir * std::sqrt(-(g00_object + g22_object * power2(photon(6)) + 2. * g03_object * photon(7) + g33_object * power2(photon(7))) / g11_object); // solved by normalization
			photon_time = fast_trace_parameters.t;
			return StatusCode::SUCCESS;
		}
		int CalcThetaPhi(const T r_observer, const T theta_observer, const T sin_theta_observer, const T cos_theta_observer, T alpha, T beta, const boost::numeric::ublas::vector<T> &u, T theta_0[], T theta_1[], T phi_0[], T phi_1[]) {
			namespace ublas = boost::numeric::ublas;
			// The photon in the observer's frame has the tetrad velocity: [1, r / R, beta / R, -alpha / R], where R = std::sqrt(r^2 + alpha^2 + beta^2).
			std::size_t N = u.size();
			ublas::bounded_vector<T, 8> photon;
			T photon_time;
			if (int status = this->InitializePhoton(photon, photon_time, alpha, beta, r_observer, power2(r_observer), theta_observer, sin_theta_observer); status != StatusCode::SUCCESS)
				return status;
			const T e_1 = 1. / Energy(photon, TimeSystem::T, DynamicalSystem::LAGRANGIAN);
			const T l = AngularMomentum(photon, TimeSystem::T, DynamicalSystem::LAGRANGIAN) * e_1, l2 = power2(l);
			const T q2 = CarterConstant(photon, 0., TimeSystem::T, DynamicalSystem::LAGRANGIAN) * power2(e_1);
			const T u_obs = 1. / r_observer;
			const T mu_obs = cos_theta_observer;
			// U=1+[a^2−q^2−l^2]u^2+2[(a−l)^2+q^2]u^3−a^2q^2u^4
			ublas::vector<T> I_u_0(N), I_u_1(N);
			I_u_1 = ublas::zero_vector<T>(N);
			ublas::vector<T> int_u_plus_0(N), int_u_plus_1(N);
			ublas::vector<T> int_u_minus_0(N), int_u_minus_1(N);
			ublas::vector<T> int_u_2_0(N), int_u_2_1(N);
			ublas::vector<T> int_u_4_0(N), int_u_4_1(N);
			for (int i = 0; i < N; ++i) {
				if (int status = UIntegral(a_, a2_, u_plus_1, u_minus_1, u_plus, u_minus, l, l2, q2, r_observer, u_obs, u(i), I_u_0(i), I_u_1(i), int_u_plus_0(i), int_u_plus_1(i), int_u_minus_0(i), int_u_minus_1(i), int_u_2_0(i), int_u_2_1(i), int_u_4_0(i), int_u_4_1(i)); status != StatusCode::SUCCESS)
					break;
			}
			// M = q^2 + (a^2-q^2-l^2)*mu^2-a^2*mu^4
			T M_minus_plus[2];
			T delta_M, delta_M_plus;
			if (int root_num_M = polySolveQuadratic(-a2_, a2_ - l2 - q2, q2, M_minus_plus); root_num_M == 0)
				return StatusCode::FAILURE;
			delta_M = std::sqrt(power2(a2_ - l2 - q2) + 4. * a2_ * q2) / a2_;
			if (M_minus_plus[1] > 1.) {
				M_minus_plus[1] = 1.;
				delta_M_plus = 0.;
			} else if (M_minus_plus[1] > 0.9)
				delta_M_plus = 2. * l2 / (a2_ + l2 + q2 + std::sqrt(power2(a2_ + l2 + q2) - 4. * a2_ * l2));
			else
				delta_M_plus = 1. - M_minus_plus[1];
			T mu_plus, mu_minus;
			T A, k, n;
			T I_mu_0, I_mu_full_turn;
			T I_t_mu_0, I_t_mu_full_turn;
			T I_phi_mu_0, I_phi_mu_full_turn;
			// Here the first turning point is determined by `beta` of the observer.
			int alpha_1, alpha_2;
			Mu0Integral(a_, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, sign(beta), mu_obs, mu_plus, mu_minus, A, k, n, alpha_1, I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn);
			T mu_f_0, mu_f_1, t_mu, phi_mu_0, phi_mu_1;
			for (int i = 0; i < N; ++i) {
				if (int status = MuFIntegral(a_, a2_, l, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, mu_plus, mu_minus, I_u_0[i], I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn, A, k, n, alpha_1, alpha_2, mu_f_0, t_mu, phi_mu_0); status != StatusCode::SUCCESS)
					return status;
				theta_0[i] = std::acos(mu_f_0);
				const T phi_u_0 = u_r * ((l * u_plus_1 + 2. * (a_ - l)) * int_u_plus_0[i] - (l * u_minus_1 + 2. * (a_ - l)) * int_u_minus_0[i]);
				phi_0[i] = phi_u_0 + phi_mu_0;
				if (std::isnan(phi_0[i]))
					return StatusCode::FAILURE;
			}
			for (int i = 0; i < N; ++i) {
				if (I_u_1[i] == 0.)
					continue;
				if (int status = MuFIntegral(a_, a2_, l, q2, M_minus_plus[1], M_minus_plus[0], delta_M, delta_M_plus, mu_plus, mu_minus, I_u_1[i], I_mu_0, I_mu_full_turn, I_t_mu_0, I_t_mu_full_turn, I_phi_mu_0, I_phi_mu_full_turn, A, k, n, alpha_1, alpha_2, mu_f_1, t_mu, phi_mu_1); status != StatusCode::SUCCESS)
					return status;
				theta_1[i] = std::acos(mu_f_1);
				const T phi_u_1 = u_r * ((l * u_plus_1 + 2. * (a_ - l)) * int_u_plus_1[i] - (l * u_minus_1 + 2. * (a_ - l)) * int_u_minus_1[i]);
				phi_1[i] = phi_u_1 + phi_mu_1;
				if (std::isnan(phi_1[i]))
					return StatusCode::FAILURE;
			}
			return StatusCode::SUCCESS;
		}
		int FastShadow(const T r_observer, const T theta_observer, const T sin_theta_observer, const T cos_theta_observer, const T alpha, const T beta, const T r_min) const override {
			boost::numeric::ublas::bounded_vector<T, 8> photon;
			T photon_time;
			if (int status = this->InitializePhoton(photon, photon_time, alpha, beta, r_observer, power2(r_observer), theta_observer, sin_theta_observer); status != StatusCode::SUCCESS)
				return status;
			const T E = Energy(photon, TimeSystem::T, DynamicalSystem::LAGRANGIAN), E_1 = 1. / E;
			const T L = AngularMomentum(photon, TimeSystem::T, DynamicalSystem::LAGRANGIAN), l = L * E_1, l2 = power2(l);
			const T Q = CarterConstant(photon, 0., TimeSystem::T, DynamicalSystem::LAGRANGIAN), q2 = Q * power2(E_1);
			const T c = a2_ - l2 - q2, d = 2. * (power2(a_ - l) + q2), e = -a2_ * q2;
			T u_roots[4];
			if (std::abs(e) < absolute_accuracy) {
				// U=1+[a^2−l^2]u^2+2(a−l)^2u^3
				// This means the observer and the object are on the equatorial plane.
				if (a_ == l) // U=1.
					return 0;
				const T d_1 = 1. / d, sqrt_d_1 = boost::math::constants::half_root_two<T>() / std::abs(a_ - l);
				if (int root_num_U = polySolveCubic(c * d_1, T{0.}, d_1, u_roots); root_num_U == 1) {
					return 0;
				}
				if (u_roots[1] >= 1.) { // photon falls into the BH
					return 0;
				}
				return 1;
			}
			const T e_1 = 1. / e, sqrt_e = std::sqrt(std::abs(e)), sqrt_e_1 = 1. / sqrt_e;
			const T c_sqrt_e = c * sqrt_e_1, d_e = d * e_1;
			if (int root_num_U = polySolveQuartic(d_e, c * e_1, T{0.}, e_1, u_roots); root_num_U == 0) { // q2 < 0., e > 0.
				return 0;
			} else if (root_num_U == 2) {
				const T f = e_1 / (u_roots[0] * u_roots[1]), g = (1. / u_roots[0] + 1. / u_roots[1]) * f; // 1. / (−a^2*q^2*u_1*u_4)
				if (q2 < 0.) {
					// root[0] < root[1] < 0.
					return 0;
				}
				if (u_roots[1] >= 1.) { // photon falls into the BH
					return 0;
				}
				return 1;
			}
			// root_num_U == 4
			if (u_roots[1] >= 1.) { // photon falls into the BH
				return 0;
			}
			return 1;
		}
		T Energy(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time, DynamicalSystem dynamics) const override {
			if (dynamics == DynamicalSystem::LAGRANGIAN) {
				if (time == TimeSystem::T)
					return (1. - 2. * y[1] / (power2(y[1]) + a2_ * power2(std::cos(y[2]))) * (1. - a_ * power2(std::sin(y[2])) * y[7])) / y[4];
				// time == TimeSystem::TAU
				return y[4] - 2. * y[1] / (power2(y[1]) + a2_ * power2(std::cos(y[2]))) * (y[4] - a_ * power2(std::sin(y[2])) * y[7]);
			} // dynamics == DynamicalSystem::HAMILTONIAN
			return 1. - y[4];
		}
		T AngularMomentum(const boost::numeric::ublas::bounded_vector<T, 8> &y, TimeSystem time, DynamicalSystem dynamics) const override {
			if (dynamics == DynamicalSystem::LAGRANGIAN) {
				const T r2 = power2(y[1]), sin2_theta = power2(std::sin(y[2]));
				const T r_rho_2 = 2. * y[1] / (r2 + a2_ * power2(std::cos(y[2])));
				if (time == TimeSystem::T)
					return (-r_rho_2 * a_ + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta / y[4];
				// time == TimeSystem::TAU
				return (-r_rho_2 * a_ * y[4] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta;
			} // dynamics == DynamicalSystem::HAMILTONIAN
			return y[7];
		}
		T CarterConstant(const boost::numeric::ublas::bounded_vector<T, 8> &y, const T mu2, TimeSystem time, DynamicalSystem dynamics) const override {
			if (dynamics == DynamicalSystem::LAGRANGIAN) {
				const T r2 = power2(y[1]), sin2_theta = power2(std::sin(y[2])), cos2_theta = power2(std::cos(y[2]));
				const T rho2 = r2 + a2_ * cos2_theta;
				const T r_rho_2 = 2. * y[1] / rho2;
				if (time == TimeSystem::T)
					return mu2 * cos2_theta * a2_ + (power2(rho2 * y[6]) + cos2_theta * (power2(-r_rho_2 * a_ + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta - a2_ * power2(r_rho_2 * (1. - a_ * sin2_theta * y[7]) - 1.))) / power2(y[4]);
				// time == TimeSystem::TAU
				return mu2 * cos2_theta * a2_ + power2(rho2 * y[6]) + cos2_theta * (power2(-r_rho_2 * a_ * y[4] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * y[7]) * sin2_theta - a2_ * power2(r_rho_2 * (y[4] - a_ * sin2_theta * y[7]) - y[4]));
			} // dynamics == DynamicalSystem::HAMILTONIAN
			return power2(y[6]) + power2(std::cos(y[2])) * (a2_ * (mu2 - power2(1. - y[4])) + power2(y[7] / std::sin(y[2])));
		}
		int NormalizeTimelikeGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y) const override {
			const T r = y[1], r2 = power2(r), a2_r2 = a2_ + r2;
			const T sin2_theta = power2(std::sin(y[2])), sin4_theta = power2(sin2_theta);
			const T rho2 = r2 + a2_ * power2(std::cos(y[2]));
			const T r_rho_2 = 2. * r / rho2;
			// y[7] += 2. * a_ * r / (power2(a2_r2) - a2_ * Delta * sin2_theta);
			y[4] = std::sqrt(1. - r_rho_2 + 2. * r_rho_2 * a_ * sin2_theta * y[7] - (rho2 / (r2 - 2. * r + a2_) * power2(y[5]) + rho2 * power2(y[6]) + (a2_r2 * sin2_theta + r_rho_2 * a2_ * sin4_theta) * power2(y[7])));
			return std::isnan(y[4]) ? StatusCode::DOMAIN_ERROR : StatusCode::SUCCESS;
		}
		int NormalizeNullGeodesic(boost::numeric::ublas::bounded_vector<T, 8> &y, T frequency = 1.) const override {
			const T r = y[1], r2 = r * r;
			const T sin2_theta = power2(std::sin(y[2]));
			const T rho2 = r2 + a2_ * power2(std::cos(y[2]));
			const T r_rho_2 = 2. * r / rho2;
			const T a = rho2 / (r2 - 2. * r + a2_) * y[5] * y[5] + rho2 * y[6] * y[6] + (a2_ + r2 + r_rho_2 * a2_ * sin2_theta) * sin2_theta * y[7] * y[7];
			const T b = -2. * r_rho_2 * a_ * sin2_theta * y[7];
			const T c = r_rho_2 - 1.;
			const T coefficient = std::copysign(0.5 / a, frequency) * (-b + std::sqrt(b * b - 4. * a * c));
			if (std::isnan(coefficient))
				return StatusCode::DOMAIN_ERROR;
			y[4] = frequency;
			y[5] *= coefficient;
			y[6] *= coefficient;
			y[7] *= coefficient;
			return StatusCode::SUCCESS;
		}
		int GetIntegrationSystem(std::function<void(const boost::numeric::ublas::bounded_vector<T, 8> &, boost::numeric::ublas::bounded_vector<T, 8> &, T)> &system, TimeSystem time, DynamicalSystem dynamics, MotionMode motion = MotionMode::GEODESIC) const override {
			if (time == TimeSystem::T) {
				if (dynamics == DynamicalSystem::LAGRANGIAN) {
					if (motion == MotionMode::GEODESIC) {
						system = [a = this->a_, a2 = this->a2_, a4 = this->a4_](const boost::numeric::ublas::bounded_vector<T, 8> &y, boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = y[6]; // d\theta/dt
							dydt[3] = y[7]; // d\phi/dt
							const T r = y[1], r2 = power2(r), r4 = power2(r2), a2_r2 = a2 + r2;
							const T sin_theta = sinTheta(y[2]), sin2_theta = power2(sin_theta), sin4_theta = power2(sin2_theta), cos_theta = cosTheta(y[2]), cos2_theta = power2(cos_theta), sin_theta_cos_theta = sin_theta * cos_theta, cot_theta = cos_theta / sin_theta;
							const T Delta = r2 - 2. * r + a2, Delta_1 = 1. / Delta;
							const T rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho4 = power2(rho2), rho_4 = power2(rho_2), rho_6 = rho_2 * rho_4, r2_a2_cos2_theta = r2 - a2 * cos2_theta;
							const T dydt4 = 2. * rho_4 * (Delta_1 * a2_r2 * r2_a2_cos2_theta * y[5] - 2. * a2 * r * sin_theta_cos_theta * y[6] * (1. - a * sin2_theta * y[7]) - Delta_1 * a * (2. * r4 + r2 * rho2 + a2 * r2_a2_cos2_theta) * sin2_theta * y[5] * y[7]);
							// d^2\tau/dt^2=-(d\tau/dt)^3*(d^2t/d\tau^2)
							dydt[4] = dydt4 * y[4];
							// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[5] = (-Delta * r2_a2_cos2_theta * power2(1. - a * sin2_theta * y[7]) - (r * (a2 * sin2_theta - r) + a2 * cos2_theta) * Delta_1 * rho4 * power2(y[5]) + 2. * a2 * sin_theta_cos_theta * rho4 * y[5] * y[6] + r * Delta * rho4 * power2(y[6]) + Delta * sin2_theta * r * rho4 * power2(y[7])) * rho_6 + dydt4 * y[5];
							// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[6] = (2. * a2 * r * sin_theta_cos_theta - a2 * sin_theta_cos_theta * rho4 * (Delta_1 * power2(y[5]) - power2(y[6])) - 2. * r * rho4 * y[5] * y[6] - 4. * a * r * sin_theta_cos_theta * a2_r2 * y[7] + sin_theta_cos_theta * (2. * a4 * r * sin4_theta + 4. * a2 * r * sin2_theta * rho2 + a2_r2 * rho4) * power2(y[7])) * rho_6 + dydt4 * y[6];
							// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
							dydt[7] = (-2. * a * r2_a2_cos2_theta * Delta_1 * y[5] + 4. * a * r * cot_theta * y[6] - 2. * Delta_1 * (r * rho4 - 2. * r2 * rho2 - r2_a2_cos2_theta * a2 * sin2_theta) * y[5] * y[7] - 2. * cot_theta * (rho4 + 2. * a2 * r * sin2_theta) * y[6] * y[7]) * rho_4 + dydt4 * y[7];
						};
						return StatusCode::SUCCESS;
					} else if (motion == MotionMode::CIRCULAR) {
						system = [](const boost::numeric::ublas::bounded_vector<T, 8> &y, boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = 0.;	// dr/dt
							dydt[2] = 0.;	// d\theta/dt
							dydt[3] = y[7]; // d\phi/dt
							dydt[4] = 0.;
							dydt[5] = 0.;
							dydt[6] = 0.;
							dydt[7] = 0.;
						};
						return StatusCode::SUCCESS;
					} else if (motion == MotionMode::CONICAL) {
						system = [a = this->a_, a2 = this->a2_](const boost::numeric::ublas::bounded_vector<T, 8> &y, boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
							dydt[0] = y[4]; // d\tau/dt
							dydt[1] = y[5]; // dr/dt
							dydt[2] = 0.;	// d\theta/dt = 0.
							dydt[3] = y[7]; // d\phi/dt
							const T r = y[1], r2 = power2(r);
							const T sin2_theta = power2(std::sin(y[2]));
							const T Delta_1 = 1. / (r2 - 2. * r + a2);
							const T rho2 = r2 + a2 * power2(std::cos(y[2])), rho_2 = 1. / rho2;
							const T r_rho_2 = 2. * y[1] * rho_2;
							const T g00 = -(1. - r_rho_2), g11 = rho2 * Delta_1, g03 = -r_rho_2 * a * sin2_theta, g33 = (r2 + a2 - g03 * a) * sin2_theta;
							const T L_1 = y[4] / (g03 + g33 * y[7]), L_2 = power2(L_1);
							const T A = g33 * (g33 * L_2 + 1.), A_1 = 1. / A, B = 2. * g03 * (g33 * L_2 + 1.), C = g00 + g11 * power2(y[5]) + power2(g03) * L_2;
							const T dg00_dr = (2. - 4. * r2 * rho_2) * rho_2, dg11_dr = 2. * (r - g11 * (r - 1.)) * Delta_1, dg03_dr = -dg00_dr * a * sin2_theta, dg33_dr = 2. * r * sin2_theta - dg03_dr * a * sin2_theta;
							const T dA_dr = dg33_dr * (2. * g33 * L_2 + 1.), dB_dr = 2. * (dg03_dr * (g33 * L_2 + 1.) + g03 * dg33_dr * L_2), dC_dr = dg00_dr + dg11_dr * power2(y[5]) + 2. * g03 * dg03_dr * L_2;
							// d^2r/dt^2 = 0.
							dydt[5] = 0.;
							// d^2\theta/dt^2 = 0.
							dydt[6] = 0.;
							// d^2\phi/dt^2 = d(d\phi/dt)/dr * dr/dt
							dydt[7] = y[5] * A_1 * (-y[7] * dA_dr + 0.5 * (-dB_dr + (B * dB_dr - 2. * (A * dC_dr + C * dA_dr)) / (2. * A * y[7] + B)));
							// d^2\tau/dt^2 = d^2\tau/dtdr * dr/dt = d((g03 + g33 * d\phi/dt) / L)/dr * dr/dt
							dydt[4] = (y[5] * (dg03_dr + dg33_dr * y[7]) + g33 * dydt[7]) * L_1;
						};
						return StatusCode::SUCCESS;
					}
				} else if (dynamics == DynamicalSystem::HAMILTONIAN && motion == MotionMode::GEODESIC) {
					system = [a = this->a_, a2 = this->a2_](const boost::numeric::ublas::bounded_vector<T, 8> &y, boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
						const T r = y[1], r2 = y[1] * y[1], a2_r2 = a2 + r2, pr2 = y[5] * y[5], ptheta2 = y[6] * y[6];
						const T E = 1. - y[4], E2 = E * E, delta_E2 = (2. - y[4]) * y[4], L2 = y[7] * y[7];
						const T sin_theta = sinTheta(y[2]), sin2_theta = sin_theta * sin_theta, sin_2_theta = 1. / sin2_theta, sin_4_theta = sin_2_theta * sin_2_theta, cos_theta = cosTheta(y[2]), cos2_theta = cos_theta * cos_theta;
						const T Delta = a2_r2 - 2. * r, Delta_1 = 1. / Delta, Delta_2 = Delta_1 * Delta_1;
						const T rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho_4 = rho_2 * rho_2;
						const T Q = ptheta2 + cos2_theta * (a2 * delta_E2 + L2 * sin_2_theta);
						const T R = -a2_r2 * r2 * delta_E2 + E2 * 2. * r * a2 + 2. * r * r2 - Delta * Q - (r2 - 2. * r) * L2 - 4. * r * a * E * y[7]; // R =power2(E * a2_r2 - a * y[7]) - Delta * (r2 +power2(y[7] - a * E) + Q);
						//[\tau,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
						dydt[0] = rho2 * Delta / (a2_r2 * (E * a2_r2 - a * y[7]) + a * Delta * (y[7] - a * E * sin2_theta));			// d\tau/dt
						dydt[1] = Delta * rho_2 * y[5] * dydt[0];																		// dr/dt
						dydt[2] = rho_2 * y[6] * dydt[0];																				// d\theta/dt
						dydt[3] = (2. * E * a * r - y[7] * (a2 - Delta * (1. + cos2_theta * sin_2_theta))) * Delta_1 * rho_2 * dydt[0]; // d\phi/dt
						dydt[4] = 0.;
						dydt[5] = ((a2 * cos2_theta + r * a2 * sin2_theta - r2) * pr2 + ((-2. * delta_E2 * r2 - a2 * delta_E2 + 3. * r - L2 - Q) * r + a2 * E2 + L2 - 2. * a * E * y[7] + Q) * rho2 * Delta_1 - ((r - 1.) * rho2 + Delta * r) * R * Delta_2) * rho_4 * dydt[0];
						dydt[6] = (-(Delta * pr2 - R * Delta_1) * a2 * rho_2 + a2 * delta_E2 + L2 * sin_2_theta + cos2_theta * sin_4_theta * L2) * sin_theta * cos_theta * rho_2 * dydt[0];
						dydt[7] = 0.;
					};
					return StatusCode::SUCCESS;
				}
			} else if (time == TimeSystem::TAU) {
				if (dynamics == DynamicalSystem::LAGRANGIAN) {
					system = [a = this->a_, a2 = this->a2_, a4 = this->a4_](const boost::numeric::ublas::bounded_vector<T, 8> &y, boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
						dydt[0] = y[4]; // dt/d\tau
						dydt[1] = y[5]; // dr/d\tau
						dydt[2] = y[6]; // d\theta/d\tau
						dydt[3] = y[7]; // d\phi/d\tau
						const T r = y[1], r2 = power2(r), a2_r2 = a2 + r2;
						const T sin_theta = sinTheta(y[2]), sin2_theta = power2(sin_theta), sin4_theta = power2(sin2_theta), cos_theta = cosTheta(y[2]), cos2_theta = power2(cos_theta), sin_theta_cos_theta = sin_theta * cos_theta, cot_theta = cos_theta / sin_theta;
						const T Delta = r2 - 2. * r + a2, Delta_1 = 1. / Delta;
						const T rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho4 = power2(rho2), rho_4 = power2(rho_2), rho_6 = rho_2 * rho_4, r2_a2_cos2_theta = r2 - a2 * cos2_theta;
						dydt[4] = -2. * rho_4 * (Delta_1 * a2_r2 * r2_a2_cos2_theta * y[5] * (y[4] - a * sin2_theta * y[7]) - 2. * a2 * r * sin_theta_cos_theta * y[6] * (y[4] - a * sin2_theta * y[7]) - 2. * Delta_1 * r2 * rho2 * a * sin2_theta * y[5] * y[7]);
						// d^2r/dt^2=(d^2r/d\tau^2)*(d\tau/dt)^2+(dr/dt)*(d^2\tau/dt^2)*(dt/d\tau)
						dydt[5] = (-Delta * r2_a2_cos2_theta * power2(y[4] - a * sin2_theta * y[7]) - (r * (a2 * sin2_theta - r) + a2 * cos2_theta) * Delta_1 * rho4 * power2(y[5]) + 2. * a2 * sin_theta_cos_theta * rho4 * y[5] * y[6] + r * Delta * rho4 * power2(y[6]) + Delta * sin2_theta * r * rho4 * power2(y[7])) * rho_6;
						// d^2\theta/dt^2=(d^2\theta/d\tau^2)*(d\tau/dt)^2+(d\theta/dt)*(d^2\tau/dt^2)*(dt/d\tau)
						dydt[6] = (2. * a2 * r * sin_theta_cos_theta * power2(y[4]) - a2 * sin_theta_cos_theta * rho4 * (Delta_1 * power2(y[5]) - power2(y[6])) - 2. * r * rho4 * y[5] * y[6] - 4. * a * r * sin_theta_cos_theta * a2_r2 * y[4] * y[7] + sin_theta_cos_theta * (2. * a4 * r * sin4_theta + 4. * a2 * r * sin2_theta * rho2 + a2_r2 * rho4) * power2(y[7])) * rho_6;
						// d^2\phi/dt^2=(d^2\phi/d\tau^2)*(d\tau/dt)^2+(d\phi/dt)*(d^2\tau/dt^2)*(dt/d\tau)
						dydt[7] = (-2. * a * r2_a2_cos2_theta * Delta_1 * y[4] * y[5] + 4. * a * r * cot_theta * y[4] * y[6] - 2. * Delta_1 * (r * rho4 - 2. * r2 * rho2 - r2_a2_cos2_theta * a2 * sin2_theta) * y[5] * y[7] - 2. * cot_theta * (rho4 + 2. * a2 * r * sin2_theta) * y[6] * y[7]) * rho_4;
					};
					return StatusCode::SUCCESS;
				} else if (dynamics == DynamicalSystem::HAMILTONIAN) {
					system = [a = this->a_, a2 = this->a2_](const boost::numeric::ublas::bounded_vector<T, 8> &y, boost::numeric::ublas::bounded_vector<T, 8> &dydt, T t) {
						const T r = y[1], r2 = power2(y[1]), a2_r2 = a2 + r2, pr2 = power2(y[5]), ptheta2 = power2(y[6]);
						const T E = 1. - y[4], E2 = power2(E), delta_E2 = (2. - y[4]) * y[4], L2 = power2(y[7]);
						const T sin_theta = sinTheta(y[2]), sin2_theta = power2(sin_theta), sin_2_theta = 1. / sin2_theta, sin_4_theta = power2(sin_2_theta), cos_theta = cosTheta(y[2]), cos2_theta = power2(cos_theta);
						const T Delta = a2_r2 - 2. * r, Delta_1 = 1. / Delta, Delta_2 = power2(Delta_1);
						const T rho2 = r2 + a2 * cos2_theta, rho_2 = 1. / rho2, rho_4 = power2(rho_2);
						const T Q = ptheta2 + cos2_theta * (a2 * delta_E2 + L2 * sin_2_theta);
						const T R = -a2_r2 * r2 * delta_E2 + E2 * 2. * r * a2 + 2. * r * r2 - Delta * Q - (r2 - 2. * r) * L2 - 4. * r * a * E * y[7]; // R =power2(E * a2_r2 - a * y[7]) - Delta * (r2 +power2(y[7] - a * E) + Q);
						//[t,r,\theta>\pi/2?\theta-\pi:\theta,\phi,1+p_t,p_r,p_\theta,p_\phi]
						dydt[0] = (a2_r2 * (E * a2_r2 - a * y[7]) * Delta_1 + a * (y[7] - a * E * sin2_theta)) * rho_2;		  // dt/d\tau
						dydt[1] = Delta * rho_2 * y[5];																		  // dr/d\tau
						dydt[2] = rho_2 * y[6];																				  // d\theta/d\tau
						dydt[3] = (2. * E * a * r - y[7] * (a2 - Delta * (1. + cos2_theta * sin_2_theta))) * Delta_1 * rho_2; // d\phi/d\tau
						dydt[4] = 0.;
						dydt[5] = ((a2 * cos2_theta + r * a2 * sin2_theta - r2) * pr2 + ((-2. * delta_E2 * r2 - a2 * delta_E2 + 3. * r - L2 - Q) * r + a2 * E2 + L2 - 2. * a * E * y[7] + Q) * rho2 * Delta_1 - ((r - 1.) * rho2 + Delta * r) * R * Delta_2) * rho_4;
						dydt[6] = (-(Delta * pr2 - R * Delta_1) * a2 * rho_2 + a2 * delta_E2 + L2 * sin_2_theta + cos2_theta * sin_4_theta * L2) * sin_theta * cos_theta * rho_2;
						dydt[7] = 0.;
					};
					return StatusCode::SUCCESS;
				}
			}
			throw std::logic_error("Kerr::GetIntegrationSystem() invaild");
		}
		T innermostStableCircularOrbitRadius(bool prograde) const override {
			const T Z1 = 1. + std::cbrt(1. - a2_) * (std::cbrt(1. + a_) + std::cbrt(1. - a_));
			const T Z2 = std::sqrt(3. * a2_ + power2(Z1));
			if (prograde)
				return 3. + Z2 - std::sqrt((3. - Z1) * (3. + Z1 + 2. * Z2));
			return 3. + Z2 + std::sqrt((3. - Z1) * (3. + Z1 + 2. * Z2));
		}
	};

	template <typename T>
	struct KerrFastTraceParameters {
		const Kerr<T> *const kerr;
		const T r, r2;
		const T u_obs, u_obj;
		const T mu_obs, mu_obj;
		const T theta_obs, sin_theta_obs;
		const T iota, sin_iota, cos_iota;
		const T theta_obs_obj;
		bool line_of_sight_axis;
		T E, L, Q;
		T t;
		T tau;
		int u_dir, mu_dir;
		KerrFastTraceParameters(const Kerr<T> *kerr, T r, T r2, T u_obs, T u_obj, T mu_obs, T mu_obj, T theta_obs, T sin_theta_obs, T iota, T sin_iota, T cos_iota, T theta_obs_obj, bool line_of_sight_axis) : kerr(kerr), r(r), r2(r2), u_obs(u_obs), u_obj(u_obj), mu_obs(mu_obs), mu_obj(mu_obj), theta_obs(theta_obs), sin_theta_obs(sin_theta_obs), iota(iota), sin_iota(sin_iota), cos_iota(cos_iota), theta_obs_obj(theta_obs_obj), line_of_sight_axis(line_of_sight_axis) {}
	};
} // namespace sbody

#endif
