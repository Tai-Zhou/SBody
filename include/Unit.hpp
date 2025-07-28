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

#ifndef SBODY_UNIT_H
#define SBODY_UNIT_H

namespace SBody {
	/**
	 * @brief SBody uses the Geometrized units, in which \f$c=G=1\f$. In addition, the central black hole mass is set to 1.
	 *
	 */
	template <typename Type>
	class Unit {
	  public:
		/// Degree, in unit of rad
		static constexpr Type deg = 1.745329251994329576923690768488612713442871888541728e-2l;
		/// Minute of arc, in unit of rad
		static constexpr Type min = 2.908882086657215961539484614147687855738119814236213e-4l;
		/// Second of arc, in unit of rad
		static constexpr Type sec = 4.8481368110953599358991410235794797595635330237270218e-6l;
		/// Millisecond of arc, in unit of rad
		static constexpr Type mas = 4.8481368110953599358991410235794797595635330237270218e-9l;
		/// Microsecond of arc, in unit of rad
		static constexpr Type uas = 4.8481368110953599358991410235794797595635330237270218e-12l;
		/// Second
		Type s;
		/// Hour
		Type hr;
		/// Day
		Type day;
		/// Year
		Type yr;
		/// Centimeter
		Type cm;
		/// Meter
		Type m;
		/// Kilometer
		Type km;
		/// Radius of the Sun, \f$R_\odot\f$
		Type R_sun;
		/// The astronomical unit
		Type AU;
		/// Milliparsec
		Type mpc;
		/// Parsec
		Type pc;
		/// Kiloparsec
		Type kpc;
		/// Megaparsec
		Type Mpc;
		/// Gigaparsec
		Type Gpc;
		/// Mass of the Sun, \f$M_\odot\f$
		Type M_sun;
		/// Mass of Earth, \f$M_\oplus\f$
		Type M_earth;
		/// Mass of Jupiter, \f$M_\mathrm{jup}\f$
		Type M_jup;
		/// Erg
		Type erg;
		/// Joules
		Type J;
		/**
		 * @brief Initialize the unit system
		 *
		 * @param mass mass of the central black hole, in unit of \f$M_\odot\f$
		 * @return status
		 */
		Unit(Type mass) {
			if (mass <= 0)
				throw std::invalid_argument("Mass must be positive.");
			M_sun = 1. / mass;
			M_earth = 3.0034893488507934e-06 * M_sun;
			M_jup = 9.545942339693249e-4 * M_sun;
			s = 2.0302544672808357e5 * M_sun;
			hr = 3600. * s;
			day = 86400. * s;
			yr = 365.25 * day;
			cm = 6.772199944005381e-6 * M_sun;
			m = 100. * cm;
			km = 1000. * m;
			R_sun = 6.957e8 * m;
			AU = 1.495978707e11 * m;
			mpc = 3.085677581491367e13 * m;
			pc = 1000. * mpc;
			kpc = 1000. * pc;
			Mpc = 1000. * kpc;
			Gpc = 1000. * Mpc;
			erg = 5.595677593689533e-55 * M_sun;
			J = 1.e7 * erg;
		}
	};
} // namespace SBody

#endif
