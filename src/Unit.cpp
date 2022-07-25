/**
 * @file Unit.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include "Unit.h"

namespace SBody {
	double Unit::M_sun = 1.;
	double Unit::M_earth = 3.0034893488507934e-06;
	double Unit::M_jup = 9.545942339693249e-4;
	double Unit::s = 2.0302544672808357e5;
	double Unit::hr = 7.3089160822110085e8;
	double Unit::day = 1.754139859730642e10;
	double Unit::yr = 6.40699583766617e12;
	double Unit::cm = 6.772199944005381e-6;
	double Unit::m = 6.772199944005381e-4;
	double Unit::km = 6.772199944005381e-1;
	double Unit::R_sun = 4.7114195010445436e5;
	double Unit::AU = 1.0131066915778643e8;
	double Unit::mpc = 2.0896825544594498e10;
	double Unit::pc = 2.0896825544594498e13;
	double Unit::kpc = 2.0896825544594498e16;
	double Unit::Mpc = 2.0896825544594498e19;
	double Unit::Gpc = 2.0896825544594498e22;
	double Unit::erg = 5.595677593689533e-55;
	double Unit::J = 5.595677593689533e-48;
	int Unit::init(double mass) {
		if (mass <= 0)
			return 1;
		M_sun = 1. / mass;
		M_earth = 3.0034893488507934e-06 * M_sun;
		M_jup = 9.545942339693249e-4 * M_sun;
		s = 2.0302544672808357e5 * M_sun;
		hr = 3600. * s;
		day = 86400. * s;
		yr = 365.25 * day;
		cm = 6.772199944005381e-6 * mass;
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
		return 0;
	}
} // namespace SBody
