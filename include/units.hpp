/**
 * @file Unit.hpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

#include <stdexcept>

namespace sbody {
/**
 * @brief SBody uses the Geometrized units, in which \f$c=G=1\f$. In addition, the central black hole mass is set to 1.
 *
 */
template <typename T>
class Unit {
  public:
    /// Mass of the Sun, \f$M_\odot\f$
    const T M_sun;
    /// Mass of Earth, \f$M_\oplus\f$
    const T M_earth;
    /// Mass of Jupiter, \f$M_\mathrm{jup}\f$
    const T M_jup;
    /// Second
    const T s;
    /// Minute
    const T min;
    /// Hour
    const T hr;
    /// Day
    const T day;
    /// Year
    const T yr;
    /// Centimeter
    const T cm;
    /// Meter
    const T m;
    /// Kilometer
    const T km;
    /// Radius of the Sun, \f$R_\odot\f$
    const T R_sun;
    /// The astronomical unit
    const T AU;
    /// Milliparsec
    const T mpc;
    /// Parsec
    const T pc;
    /// Kiloparsec
    const T kpc;
    /// Megaparsec
    const T Mpc;
    /// Gigaparsec
    const T Gpc;
    /// Erg
    const T erg;
    /// Joules
    const T J;
    /**
     * @brief Initialize the unit system
     *
     * @param mass mass of the central black hole, in unit of \f$M_\odot\f$
     * @return status
     */
    Unit(T mass)
        : M_sun(1. / mass), M_earth(3.0034893488507934e-06 * M_sun), M_jup(9.545942339693249e-4 * M_sun),
          s(2.0302544672808357e5 * M_sun), min(60. * s), hr(60. * min), day(86400. * s), yr(365.25 * day),
          cm(6.772199944005381e-6 * M_sun), m(100. * cm), km(1000. * m), R_sun(6.957e8 * m), AU(1.495978707e11 * m),
          mpc(3.085677581491367e13 * m), pc(1000. * mpc), kpc(1000. * pc), Mpc(1000. * kpc), Gpc(1000. * Mpc),
          erg(5.595677593689533e-55 * M_sun), J(1.e7 * erg) {
        if (mass <= 0)
            throw std::invalid_argument("Mass must be positive.");
    }
};
} // namespace sbody
