#ifndef SBODY_UNIT_H
#define SBODY_UNIT_H

namespace SBody {
	class Unit {
	  public:
		static int init(double mass);
		static constexpr double deg = 1.7453292519943295e-2;
		static constexpr double min = 2.908882086657216e-4;
		static constexpr double sec = 4.84813681109536e-6;
		static constexpr double mas = 4.84813681109536e-9;
		static constexpr double uas = 4.84813681109536e-12;
		// Time
		static double s;
		static double day;
		static double yr;
		// Length
		static double cm;
		static double m;
		static double R_sun;
		static double AU;
		static double mpc;
		static double pc;
		static double kpc;
		static double Mpc;
		static double Gpc;
		// Mass
		static double M_sun;
		static double M_earth;
		static double M_jup;
		// Energy
		static double erg;
		static double J;
	};
} // namespace SBody

#endif
