#ifndef SBODY_CONSTANT_H
#define SBODY_CONSTANT_H

namespace SBody {
	// Math
	constexpr double M_2PI = 6.28318530717958647692528676655900576; // 2*pi
	constexpr double M_PI2 = 9.86960440108935861883449099987615111; // pi^2
	namespace Constant {
		// Physical
		constexpr double e = 1.60217662e-19;
		constexpr double h = 7.52764601e-76;
		// Time
		constexpr double s = 2.0302544672808357e5;
		constexpr double day = 1.754139859730642e10;
		constexpr double yr = 6.40699583766617e12;
		// Length
		constexpr double cm = 6.772199944005381e-06;
		constexpr double R_sun = 4.7114195010445436e5;
		constexpr double AU = 1.0131066915778643e8;
		constexpr double mpc = 2.0896825544594498e10;
		constexpr double pc = 2.0896825544594498e13;
		constexpr double kpc = 2.0896825544594498e16;
		constexpr double Mpc = 2.0896825544594498e19;
		constexpr double Gpc = 2.0896825544594498e22;
		// Mass
		constexpr double M_earth = 3.00348935e-6;
		constexpr double M_jup = 9.5459e-4;
		// Energy
		constexpr double erg = 5.59567759e-55;
		constexpr double J = 5.59567759e-48;
	} // namespace Constant
} // namespace SBody

#endif
