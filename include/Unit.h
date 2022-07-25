/**
 * @file Unit.h
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_UNIT_H
#define SBODY_UNIT_H

namespace SBody {
	class Unit {
	  public:
		/**
		 * @brief Initialize the unit system
		 *
		 * @param mass mass of the central black hole, in unit of \f$M_\odot\f$
		 * @return int
		 */
		static int init(double mass);
		/// Degree, in unit of rad
		static constexpr double deg = 1.7453292519943295e-2;
		/// Minute of arc, in unit of rad
		static constexpr double min = 2.908882086657216e-4;
		/// Second of arc, in unit of rad
		static constexpr double sec = 4.84813681109536e-6;
		/// Millisecond of arc, in unit of rad
		static constexpr double mas = 4.84813681109536e-9;
		/// Microsecond of arc, in unit of rad
		static constexpr double uas = 4.84813681109536e-12;
		/// Second
		static double s;
		/// Hour
		static double hr;
		/// Day
		static double day;
		/// Year
		static double yr;
		/// Centimeter
		static double cm;
		/// Meter
		static double m;
		/// Kilometer
		static double km;
		/// Radius of the Sun, \f$R_\odot\f$
		static double R_sun;
		/// The astronomical unit
		static double AU;
		/// Milliparsec
		static double mpc;
		/// Parsec
		static double pc;
		/// Kiloparsec
		static double kpc;
		/// Megaparsec
		static double Mpc;
		/// Gigaparsec
		static double Gpc;
		/// Mass of the Sun, \f$M_\odot\f$
		static double M_sun;
		/// Mass of Earth, \f$M_\oplus\f$
		static double M_earth;
		/// Mass of Jupiter, \f$M_{jup}\f$
		static double M_jup;
		/// Erg
		static double erg;
		/// Joules
		static double J;
	};
} // namespace SBody

#endif
