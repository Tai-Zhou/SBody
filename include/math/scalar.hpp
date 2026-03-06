/**
 * @file scalar.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-03
 *
 * @copyright Copyright (c) 2026
 */

#ifndef SBODY_SCALAR_HPP
#define SBODY_SCALAR_HPP

namespace sbody {
	template <typename T>
	T sign(T x) {
		if (x >= 0.)
			return 1.;
		return -1.;
	}

	template <typename T>
	T power2(T x) {
		return x * x;
	}

	template <typename T>
	T power3(T x) {
		return x * x * x;
	}

	template <typename T>
	T power4(T x) {
		x *= x;
		return x * x;
	}

	/**
	 * @brief Square root of `x`. Return `0` if `x` is negative.
	 *
	 * @param x number
	 * @return result
	 */
	template <typename T>
	T squareRoot(T x) {
		if (x <= 0.)
			return 0.;
		return std::sqrt(x);
	}

	template <typename T>
	T hypotSquare(T x, T y) {
		return x * x + y * y;
	}
	template <typename T>
	T hypotSquare(T x, T y, T z) {
		return x * x + y * y + z * z;
	}

	/**
	 * @brief Return `1` if `x`, `y` have opposite signs, else `0`.
	 *
	 * @param x number
	 * @param y number
	 * @return result
	 */
	template <typename T>
	bool oppositeSign(T x, T y) {
		return (x > 0. && y < 0.) || (x < 0. && y > 0.);
	}

	template <typename T>
	void aMinusPlusB(T a, T b, T &a_minus_b, T &a_plus_b) {
		a_minus_b = a - b;
		a_plus_b = a + b;
	}
} // namespace sbody

#endif
