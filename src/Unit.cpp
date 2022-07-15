#include "Unit.h"

namespace SBody {
	double Unit::M_sun = 1.;
	double Unit::M_earth = 3.0034893488507934e-06;
	double Unit::M_jup = 9.545942339693249e-4;
	double Unit::s = 2.0302544672808357e5;
	double Unit::day = 86400. * s;
	double Unit::yr = 365.25 * day;
	double Unit::cm = 6.772199944005381e-6;
	double Unit::m = 6.772199944005381e-4;
	double Unit::R_sun = 6.957e8 * m;
	double Unit::AU = 1.495978707e11 * m;
	double Unit::mpc = 3.085677581491367e13 * m;
	double Unit::pc = 1000. * mpc;
	double Unit::kpc = 1000. * pc;
	double Unit::Mpc = 1000. * kpc;
	double Unit::Gpc = 1000. * Mpc;
	double Unit::erg = 5.595677593689533e-55;
	double Unit::J = 1.e7 * erg;
	void Unit::init(double mass) {
		M_sun = 1. / mass;
		M_earth = 3.0034893488507934e-06 * M_sun;
		M_jup = 9.545942339693249e-4 * M_sun;
		s = 2.0302544672808357e5 * mass;
		day = 86400. * s;
		yr = 365.25 * day;
		cm = 6.772199944005381e-6 * mass;
		m = 100. * cm;
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
} // namespace SBody
