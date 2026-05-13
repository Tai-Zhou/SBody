/**
 * @file Utility.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <cstddef>

#include <boost/array.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fmt/core.h>

#include "math/geometry.hpp"

namespace sbody {
/**
 * @brief Linear interpolation of points (`x0`, `y0`) and (`x1`, `y1`) at `x`. \f[y=\frac{y_0(x_1-x)+y_1(x-x_0)}{x_1-x_0}\f]
 *
 * @param x position to evaluate \f$y\f$
 * @param x0 \f$x_0\f$
 * @param x1 \f$x_1\f$
 * @param y0 \f$y_0\f$
 * @param y1 \f$y_1\f$
 * @return result
 */
template <typename T>
T linearInterpolation(T x, T x0, T x1, T y0, T y1) {
    if (x0 == x1)
        return 0.5 * (y0 + y1);
    return (y0 * (x1 - x) + y1 * (x - x0)) / (x1 - x0);
}

/**
 * @brief Linear interpolation of vectors (`x0`, `y0`) and (`x1`, `y1`) at `x`, stored in `y`.
 * \f[y=\frac{y_0(x_1-x)+y_1(x-x_0)}{x_1-x_0}\f]
 *
 * @param x position to evaluate \f$y\f$
 * @param x0 \f$x_0\f$
 * @param x1 \f$x_1\f$
 * @param y0 \f$y_0\f$
 * @param y1 \f$y_1\f$
 * @param y \f$y\f$
 * @param size size of the vectors
 * @return status
 */
template <typename T, std::size_t N>
void linearInterpolation(T x, T x0, T x1, const boost::numeric::ublas::bounded_vector<T, N> &y0,
                         const boost::numeric::ublas::bounded_vector<T, N> &y1,
                         boost::numeric::ublas::bounded_vector<T, N> &y) {
    if (x0 == x1)
        y = 0.5 * (y0 + y1);
    else
        y = ((x - x0) * y1 + (x1 - x) * y0) / (x1 - x0);
}
template <typename T, std::size_t N>
void linearInterpolation(T x, T x0, T x1, const boost::numeric::ublas::bounded_vector<T, N> &y0,
                         const boost::numeric::ublas::bounded_vector<T, N> &y1, double y[]) {
    if (x0 == x1) {
        for (std::size_t i = 0; i < N; ++i)
            y[i] = 0.5 * (y0(i) + y1(i));
        return;
    }
    const T x01_1 = 1. / (x1 - x0);
    for (std::size_t i = 0; i < N; ++i)
        y[i] = (x - x0) * x01_1 * y1(i) + (x1 - x) * x01_1 * y0(i);
}

/**
 * @brief Linear interpolation of two spherical positions `y0` and `y1` at `t`, stored in `y`.
 * \f[y=\frac{\text{SphericalToCartesian}(y_0)(t_1-t)+\text{SphericalToCartesian}(y_1)(t-t_0)}{t_1-t_0}\f]
 *
 * @param t time to evaluate \f$y\f$
 * @param t0 \f$t_0\f$
 * @param t1 \f$t_1\f$
 * @param y0 8 dimensional vector, spherical, \f$y_0\f$
 * @param y1 8 dimensional vector, spherical, \f$y_1\f$
 * @param y 8 dimensional vector, cartesian, \f$y\f$
 * @return status
 */
template <typename T>
void interpolateSphericalPositionToCartesian(T t, T t0, T t1, const boost::numeric::ublas::bounded_vector<T, 8> &y0,
                                             const boost::numeric::ublas::bounded_vector<T, 8> &y1,
                                             boost::numeric::ublas::bounded_vector<T, 8> &y) {
    boost::numeric::ublas::bounded_vector<T, 8> c0, c1;
    sphericalToCartesian(y0, c0);
    sphericalToCartesian(y1, c1);
    linearInterpolation(t, t0, t1, c0, c1, y);
}
} // namespace sbody
