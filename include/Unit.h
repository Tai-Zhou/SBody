#ifndef SBODY_UNIT_H
#define SBODY_UNIT_H

#include "Constant.h"

namespace SBody {
	class Unit {
	  public:
		static int init(double mass);
		static double M_sun;
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
		static double M_earth;
		static double M_jup;
		// Energy
		static double erg;
		static double J;
	};
} // namespace SBody

#endif
