/**
 * @file special_functions.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2026-03-04
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>

#include "math/scalar.hpp"

namespace sbody {
template <typename T>
T carlsonRC(T x, T y) {
    if (y < 0.)
        return std::sqrt(x / (x - y)) * boost::math::ellint_rc(x - y, -y);
    return boost::math::ellint_rc(x, y);
}

template <typename T>
T carlsonRJ(T x, T y, T z, T p) {
    if (p <= 0.) {
        if ((x >= y && z >= x) || (x <= y && z <= x))
            std::swap(x, y);
        else if ((z >= y && x >= z) || (z <= y && x <= z))
            std::swap(y, z);
        const T y_1 = 1. / y, y_p_1 = 1. / (y - p), gamma = y + (z - y) * (y - x) * y_p_1; // make sure gamma > 0.
        return ((gamma - y) * boost::math::ellint_rj(x, y, z, gamma) -
                3. * (boost::math::ellint_rf(x, y, z) - carlsonRC(x * z * y_1, p * gamma * y_1))) *
               y_p_1;
    }
    return boost::math::ellint_rj(x, y, z, p);
}
template <typename T>
T ellipticIntegral(int p5, T y, T x, T a5, T b5, T a1, T b1, T a2, T b2, T a3, T b3, T a4 = 1., T b4 = 0.) {
    if (x == y)
        return 0.;
    if (p5 != 0 && p5 != -2 && p5 != -4)
        return std::numeric_limits<T>::quiet_NaN();
    const T d12 = a1 * b2 - a2 * b1, d13 = a1 * b3 - a3 * b1, d14 = a1 * b4 - a4 * b1;
    const T d23 = a2 * b3 - a3 * b2, d24 = a2 * b4 - a4 * b2, d34 = a3 * b4 - a4 * b3;
    const T X1 = std::sqrt(a1 + b1 * x), X2 = std::sqrt(a2 + b2 * x), X3 = std::sqrt(a3 + b3 * x), X4 = std::sqrt(a4 + b4 * x),
            X52 = a5 + b5 * x;
    const T Y1 = std::sqrt(a1 + b1 * y), Y2 = std::sqrt(a2 + b2 * y), Y3 = std::sqrt(a3 + b3 * y), Y4 = std::sqrt(a4 + b4 * y),
            Y52 = a5 + b5 * y;
    const T U2_12 = power2((X1 * X2 * Y3 * Y4 + Y1 * Y2 * X3 * X4) / (x - y));
    const T U2_13 = U2_12 - d14 * d23;
    const T U2_14 = U2_12 - d13 * d24;
    const T I1 = 2. * boost::math::ellint_rf(U2_12, U2_13, U2_14);
    if (p5 == 0)
        return I1;
    const T b52 = b5 * b5;
    const T d15 = a1 * b5 - a5 * b1, d15_1 = 1. / d15, d25 = a2 * b5 - a5 * b2, d35 = a3 * b5 - a5 * b3,
            d45 = a4 * b5 - a5 * b4;
    const T W2 = U2_12 - d13 * d14 * d25 * d15_1;
    const T Q2 = X52 * Y52 * W2 / power2(X1 * Y1);
    const T P2 = Q2 + d25 * d35 * d45 * d15_1;
    const T RC_P2_Q2 = X1 * Y1 == 0. ? 0. : carlsonRC(P2, Q2);
    const T I3 = 2. * (d12 * d13 * d14 * d15_1 / 3. * carlsonRJ(U2_12, U2_13, U2_14, W2) + RC_P2_Q2);
    if (p5 == -2)
        return (b5 * I3 - b1 * I1) * d15_1;
    const T I2 = 2. * (d12 * d13 * boost::math::ellint_rd(U2_12, U2_13, U2_14) / 3. + X1 * Y1 / (X4 * Y4 * std::sqrt(U2_14)));
    return -0.5 * d15_1 * (b1 / d15 + b2 / d25 + b3 / d35 + b4 / d45) * b5 * I3 +
           b52 * d24 * d34 / (2. * d15 * d25 * d35 * d45) * I2 +
           power2(b1 * d15_1) * (1. - d12 * d13 * b52 / (2. * b1 * b1 * d25 * d35)) * I1 -
           b52 / (d15 * d25 * d35) * (X1 * X2 * X3 / (X4 * X52) - Y1 * Y2 * Y3 / (Y4 * Y52));
}

template <typename T>
T ellipticIntegral2Complex(int p5, T y, T x, T a5, T b5, T f, T g, T h, T a1, T b1, T a4 = 1., T b4 = 0.) {
    if (x == y)
        return 0.;
    if (p5 != 0 && p5 != -2 && p5 != -4)
        return std::numeric_limits<T>::quiet_NaN();
    const T b12 = b1 * b1, b52 = b5 * b5;
    const T X1 = std::sqrt(a1 + b1 * x), X4 = std::sqrt(a4 + b4 * x);
    const T Y1 = std::sqrt(a1 + b1 * y), Y4 = std::sqrt(a4 + b4 * y);
    const T xi = std::sqrt(f + (g + h * x) * x), eta = std::sqrt(f + (g + h * y) * y);
    const T M2 = power2(X1 * Y4 + X4 * Y1) * (power2((xi + eta) / (x - y)) - h);
    const T c2_11 = 2 * (f * b12 - g * a1 * b1 + h * a1 * a1), c2_44 = 2 * (f * b4 * b4 - g * a4 * b4 + h * a4 * a4);
    const T c2_14 = 2. * f * b1 * b4 - g * (a1 * b4 + a4 * b1) + 2. * h * a1 * a4,
            c2_15 = 2. * f * b1 * b5 - g * (a1 * b5 + a5 * b1) + 2. * h * a1 * a5;
    const T c11 = std::sqrt(c2_11), c44 = std::sqrt(c2_44), c11_c44 = c11 * c44;
    const T L2m = std::max(T{0.}, M2 + c2_14 - c11_c44), L2p = std::max(T{0.}, M2 + c2_14 + c11_c44);
    const T I1 = 4. * boost::math::ellint_rf(M2, L2m, L2p);
    if (p5 == 0)
        return I1;
    const T X52 = a5 + b5 * x, Y52 = a5 + b5 * y;
    const T c2_55 = 2 * (f * b52 - g * a5 * b5 + h * a5 * a5), c55 = std::sqrt(c2_55);
    const T d14 = a1 * b4 - a4 * b1, d15 = a1 * b5 - a5 * b1, d15_1 = 1. / d15, d45 = a4 * b5 - a5 * b4;
    const T W2p = M2 + d14 * (c2_15 + c11 * c55) * d15_1;
    const T U = (X1 * X4 * eta + Y1 * Y4 * xi) / (x - y), U2 = U * U;
    const T W2 = U2 - 0.5 * c2_11 * d45 * d15_1;
    const T Q2 = X52 * Y52 * W2 / power2(X1 * Y1);
    const T P2 = Q2 + 0.5 * c2_55 * d45 * d15_1;
    const T RC_P2_Q2 = X1 * Y1 == 0. ? 0. : carlsonRC(P2, Q2);
    const T I3 = 2. * (c11 / (3. * c55) * (4. * (W2p - M2) * carlsonRJ(M2, L2m, L2p, W2p) - 1.5 * I1 + 3. * carlsonRC(U2, W2)) +
                       RC_P2_Q2);
    if (p5 == -2)
        return (b5 * I3 - b1 * I1) * d15_1;
    const T I2 = 2. * (c11 / (3. * c44) * (4. * (c2_14 + c11_c44) * boost::math::ellint_rd(M2, L2m, L2p) - 1.5 * I1 + 3. / U) +
                       X1 * Y1 / (X4 * Y4 * U));
    return -0.5 * d15_1 * ((b1 * b5) / d15 + 2. * b5 * (g * b5 - 2. * h * a5) / c2_55 + (b4 * b5) / d45) * I3 +
           b52 * 0.5 * c2_44 / (d15 * d45 * c2_55) * I2 + power2(b1 * d15_1) * (1. - 0.5 * c2_11 * b52 / (b12 * c2_55)) * I1 -
           b52 / (0.5 * d15 * c2_55) * (X1 * xi / (X4 * X52) - Y1 * eta / (Y4 * Y52));
}
template <typename T>
T ellipticIntegral4Complex(int p5, T y, T x, T a5, T b5, T f1, T g1, T h1, T f2, T g2, T h2) {
    if (x == y)
        return 0.;
    if (p5 != 0 && p5 != -2 && p5 != -4)
        return std::numeric_limits<T>::quiet_NaN();
    const T x_y_1 = 1. / (x - y);
    const T xi12 = f1 + (g1 + h1 * x) * x, eta12 = f1 + (g1 + h1 * y) * y;
    const T xi22 = f2 + (g2 + h2 * x) * x, eta22 = f2 + (g2 + h2 * y) * y;
    const T xi1 = std::sqrt(xi12), eta1 = std::sqrt(eta12);
    const T xi2 = std::sqrt(xi22), eta2 = std::sqrt(eta22);
    const T xi1p = (g1 + 2. * h1 * x) / (2. * xi1), eta_1p = (g1 + 2. * h1 * y) / (2. * eta1);
    const T B = xi1p * xi2 - eta_1p * eta2;
    const T theta1 = xi12 + eta12 - h1 * power2(x - y);
    const T theta2 = xi22 + eta22 - h2 * power2(x - y);
    const T zeta1 = std::sqrt(2. * xi1 * eta1 + theta1), zeta2 = std::sqrt(2. * xi2 * eta2 + theta2);
    const T U = (xi1 * eta2 + xi2 * eta1) * x_y_1, U2 = U * U;
    const T M = zeta1 * zeta2 * x_y_1;
    const T M2 = M * M;
    const T delta11_2 = 4. * f1 * h1 - g1 * g1, delta12_2 = 2. * (f1 * h2 + f2 * h1) - g1 * g2,
            delta22_2 = 4. * f2 * h2 - g2 * g2;
    const T Delta = std::sqrt(delta12_2 * delta12_2 - delta11_2 * delta22_2);
    const T Delta_m = delta12_2 - Delta, Delta_p = delta12_2 + Delta;
    const T L2m = M2 + Delta_m, L2p = M2 + Delta_p;
    const T RF = boost::math::ellint_rf(M2, L2m, L2p);
    if (p5 == 0)
        return 4. * RF;
    const T G = 2. * Delta * Delta_p * boost::math::ellint_rd(M2, L2m, L2p) / 3. + Delta / (2. * U) +
                (delta12_2 * theta1 - delta11_2 * theta2) / (4. * xi1 * eta1 * U);
    const T Sigma = G - Delta_p * RF + B;
    const T alpha15 = 2. * f1 * b5 - g1 * a5, beta15 = g1 * b5 - 2. * h1 * a5;
    const T alpha25 = 2. * f2 * b5 - g2 * a5, beta25 = g2 * b5 - 2. * h2 * a5;
    const T gamma1 = 0.5 * (alpha15 * b5 - beta15 * a5), gamma2 = 0.5 * (alpha25 * b5 - beta25 * a5);
    const T gamma1_1 = 1. / gamma1, gamma2_1 = 1. / gamma2;
    const T Lambda = delta11_2 * gamma2 * gamma1_1;
    const T Omega2 = M2 + Lambda;
    const T psi = 0.5 * (alpha15 * beta25 - alpha25 * beta15);
    const T xi5 = a5 + b5 * x, eta5 = a5 + b5 * y;
    const T X = 0.5 * x_y_1 * (xi5 * (alpha15 + beta15 * y) * eta2 / eta1 + eta5 * (alpha15 + beta15 * x) * xi2 / xi1);
    const T S = 0.5 * (M2 + delta12_2) - U2, S2 = S * S;
    const T mu = gamma1 * xi5 * eta5 / (xi1 * eta1);
    const T T2 = power2(mu * S + 2. * gamma1 * gamma2);
    const T V2 = mu * mu * (S2 + Lambda * U2);
    const T a = S * Omega2 / U + 2. * Lambda * U, a2 = a * a;
    const T b2 = (S2 / U2 + Lambda) * Omega2 * Omega2;
    const T H = delta11_2 * psi * gamma1_1 * gamma1_1 * (carlsonRJ(M2, L2m, L2p, Omega2) / 3. + 0.5 * carlsonRC(a2, b2)) -
                X * carlsonRC(T2, V2);
    if (p5 == -2)
        return -2. * (b5 * H + beta15 * RF * gamma1_1);
    return b5 * (beta15 * gamma1_1 + beta25 * gamma2_1) * H + power2(beta15 * gamma1_1) * RF +
           b5 * b5 * (Sigma - b5 * (xi1 * xi2 / xi5 - eta1 * eta2 / eta5)) * gamma1_1 * gamma2_1;
}
} // namespace sbody
