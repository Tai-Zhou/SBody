/**
 * @file Python.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <array>
#include <chrono>
#include <cmath>
#include <csignal>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include <gsl/gsl_math.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include "IO.h"
#include "Metric.h"
#include "Object.h"
#include "Unit.h"
#include "Utility.h"
#include "View.h"

namespace py = pybind11;
using namespace std;
using namespace SBody;

int Fibo(int n) {
	int i = 0, j = 1, k;
	while (n--) {
		k = j;
		j = (i + j) % 1000000007;
		i = k;
	}
	return i;
}

void saveFibo(int n) {
	SBody::NumPy record("Fibo", {1});
	int i = 0, j = 1, k;
	while (n--) {
		k = j;
		j = (i + j) % 1000000007;
		i = k;
		record.Save({sqrt(i)});
	}
}

py::array_t<double> MCMC(double mass, int metric, double PN, double R, double tp, double a = 5., double e = 0.884649, double inclination = 134.567, double ascendingNode = 228.171, double periapsis = 66.263) {
	double spin = 0.;
	vector<double> tList = {1992.224, 1994.314, 1995.534, 1996.253,
							1996.427, 1997.544, 1998.373, 1999.465,
							2000.472, 2000.487, 2000.523, 2001.502,
							2002.25, 2002.335, 2002.393, 2002.409,
							2002.412, 2002.414, 2002.41438188, 2002.418,
							2002.421, 2002.488, 2002.578, 2002.66,
							2002.66, 2003.214, 2003.271, 2003.34419768,
							2003.351, 2003.353, 2003.356, 2003.446,
							2003.446, 2003.451, 2003.452, 2003.454,
							2003.454, 2003.55, 2003.676, 2003.678,
							2003.761, 2004.24, 2004.325, 2004.347,
							2004.443, 2004.511, 2004.513, 2004.516,
							2004.516, 2004.535, 2004.537, 2004.574,
							2004.574, 2004.632, 2004.664, 2004.67,
							2004.73, 2004.73, 2005.15264675, 2005.158,
							2005.212, 2005.215, 2005.27, 2005.366,
							2005.371, 2005.374, 2005.455, 2005.461,
							2005.467, 2005.57, 2005.576, 2005.677,
							2005.769, 2006.204, 2006.305, 2006.324,
							2006.624, 2007.21700843, 2007.23, 2007.304,
							2007.545, 2007.55, 2007.55, 2007.686,
							2007.687, 2007.73, 2008.148, 2008.197,
							2008.262, 2008.268, 2008.431, 2008.456,
							2008.472, 2008.473, 2008.593, 2008.601,
							2008.708, 2009.185, 2009.273, 2009.3,
							2009.303, 2009.336, 2009.336, 2009.371,
							2009.38373515, 2009.385, 2009.502, 2009.505,
							2009.557, 2009.557, 2009.606, 2009.718,
							2009.776, 2010.234, 2010.239, 2010.239,
							2010.245, 2010.351, 2010.354, 2010.444,
							2010.455, 2010.455, 2010.455, 2010.46,
							2010.616, 2010.619, 2010.622, 2010.624,
							2010.627, 2010.676, 2010.679, 2011.238,
							2011.241, 2011.244, 2011.249, 2011.312,
							2011.313, 2011.315, 2011.317, 2011.31739855,
							2011.337, 2011.443, 2011.553, 2011.567,
							2011.613, 2011.689, 2011.695, 2011.695,
							2011.698, 2011.722, 2012.202, 2012.21,
							2012.339, 2012.342, 2012.494, 2012.497,
							2012.513, 2012.533, 2012.544, 2012.552,
							2012.552, 2012.604, 2012.7, 2012.70150083,
							2012.705, 2013.161, 2013.24, 2013.262,
							2013.317, 2013.366, 2013.42, 2013.437,
							2013.494, 2013.502, 2013.587, 2013.59,
							2013.617, 2013.655, 2013.726, 2014.185,
							2014.263, 2014.51917732, 2014.521, 2015.299,
							2015.432, 2015.517, 2015.706, 2015.706,
							2015.747, 2016.221, 2016.284, 2016.287,
							2016.325, 2016.369, 2016.39849, 2016.519,
							2016.525, 2016.53, 2016.72579643, 2017.2394057,
							2017.28599, 2017.50222741, 2017.5449, 2017.60102297,
							2017.71494, 2018.22132, 2018.24783858, 2018.32475,
							2018.33844779, 2018.40972, 2018.41130907, 2018.483251,
							2018.53875, 2018.56407834, 2018.62677, 2018.65422043,
							2019.27084392, 2019.32185, 2019.42827, 2019.50532956,
							2019.62025885, 2019.69998329};
	double inc = M_PI * 0. / 180., eps = 0.;
	Unit::Initialize(mass);
	a *= 4.84813681109536e-09 * R * Unit::pc;
	inclination *= M_PI / 180.;
	ascendingNode *= M_PI / 180.;
	periapsis *= M_PI / 180.;
	double trueAnomaly = M_PI;
	double t = tp - sqrt(gsl_pow_3(a / Unit::AU) / mass) * 0.5;
	unique_ptr<View> viewPtr;
	// metric::setMetric(metric, PN != 0 ? 1 : 0, mass, 0, 0, 0);
	Newton main_metric(PN, T);
	char strFormat[1024]; // TODO: waiting for C++20
	snprintf(strFormat, 1024, " (%.1f,%.1f,%.1f)[%f,%f]", spin, 0., 0., inc, eps);
	if (metric == 1)
		viewPtr = make_unique<View>(make_unique<Schwarzschild>(HAMILTONIAN), 8246.7 * Unit::pc, inc, string("view") + strFormat);
	double x[8], y[8];
	if (metric == 0) {
		y[0] = 0.;
		y[4] = 1.;
		double r = a * (1 - e * e) / (1 + e * cos(trueAnomaly));
		double tp1 = -r * cos(periapsis + trueAnomaly), tp2 = -r * sin(periapsis + trueAnomaly) * cos(inclination);
		double xp1 = tp1 * cos(ascendingNode) - tp2 * sin(ascendingNode), xp2 = tp2 * cos(ascendingNode) + tp1 * sin(ascendingNode), xp3 = -r * sin(periapsis + trueAnomaly) * sin(inclination);
		y[1] = (xp1 * cos(eps) + xp2 * sin(eps)) * cos(inc) + xp3 * sin(inc);
		y[2] = xp2 * cos(eps) - xp1 * sin(eps);
		y[3] = xp3 * cos(inc) - (xp1 * cos(eps) + xp2 * sin(eps)) * sin(inc);
		double vtheta = sqrt((1 - e * e) * mass * a) / r, vr = GSL_SIGN(M_PI - ModBy2Pi(trueAnomaly)) * sqrt(max(0., 2. * mass / r - mass / a - vtheta * vtheta));
		double tp5 = vtheta * sin(periapsis + trueAnomaly) - vr * cos(periapsis + trueAnomaly), tp6 = -(vtheta * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * cos(inclination);
		double xp5 = tp5 * cos(ascendingNode) - tp6 * sin(ascendingNode), xp6 = tp5 * sin(ascendingNode) + tp6 * cos(ascendingNode), xp7 = -(vtheta * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * sin(inclination);
		y[5] = (xp5 * cos(eps) + xp6 * sin(eps)) * cos(inc) + xp7 * sin(inc);
		y[6] = xp6 * cos(eps) - xp5 * sin(eps);
		y[7] = xp7 * cos(inc) - (xp5 * cos(eps) + xp6 * sin(eps)) * sin(inc);
	} else {
		x[0] = 0.;
		double r = a * (1 - e * e) / (1 + e * cos(trueAnomaly));
		double tp1 = -r * cos(periapsis + trueAnomaly), tp2 = -r * sin(periapsis + trueAnomaly) * cos(inclination);
		double xp1 = tp1 * cos(ascendingNode) - tp2 * sin(ascendingNode), xp2 = tp2 * cos(ascendingNode) + tp1 * sin(ascendingNode), xp3 = -r * sin(periapsis + trueAnomaly) * sin(inclination);
		x[1] = (xp1 * cos(eps) + xp2 * sin(eps)) * cos(inc) + xp3 * sin(inc);
		x[2] = xp2 * cos(eps) - xp1 * sin(eps);
		x[3] = xp3 * cos(inc) - (xp1 * cos(eps) + xp2 * sin(eps)) * sin(inc);
		double vtheta = sqrt((1 - e * e) * mass * a) / r, vr = GSL_SIGN(M_PI - ModBy2Pi(trueAnomaly)) * sqrt(max(0., 2. * mass / r - mass / a - vtheta * vtheta));
		double tp5 = vtheta * sin(periapsis + trueAnomaly) - vr * cos(periapsis + trueAnomaly), tp6 = -(vtheta * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * cos(inclination);
		double xp5 = tp5 * cos(ascendingNode) - tp6 * sin(ascendingNode), xp6 = tp5 * sin(ascendingNode) + tp6 * cos(ascendingNode), xp7 = -(vtheta * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * sin(inclination);
		x[5] = (xp5 * cos(eps) + xp6 * sin(eps)) * cos(inc) + xp7 * sin(inc);
		x[6] = xp6 * cos(eps) - xp5 * sin(eps);
		x[7] = xp7 * cos(inc) - (xp5 * cos(eps) + xp6 * sin(eps)) * sin(inc);
		metric::c2s(x, y);
		main_metric.NormalizeTimelikeGeodesic(y);
	}
	Integrator &&integrator = main_metric.GetIntegrator(PN);
	Star star_0(make_unique<Newton>(PN, T), Unit::R_sun, y, 0);
	auto result = py::array_t<double>(tList.size() * 8);
	double *result_ptr = (double *)result.request().ptr;
	int status = 0, rayNO = 0;
	double h = -1., tPoint = tList[0] * Unit::yr;
	while (t > tPoint)
		status = integrator.Apply(&t, tPoint, &h, star_0.pos);
	h = 1.;
	integrator.Reset();
	for (double tPoint : tList) {
		tPoint = tPoint * Unit::yr;
		while (status <= 0 && t < tPoint)
			status = integrator.Apply(&t, tPoint, &h, star_0.pos);
		if (status > 0)
			py::print("[!] status =", status);
		if (metric == 0)
			copy(star_0.pos, star_0.pos + 8, result_ptr + (rayNO++) * 8);
		else
			viewPtr->TraceStar(star_0, rayNO++);
	}
	return result.reshape({tList.size(), 8UL});
}

PYBIND11_MODULE(SBoPy, m) {
	m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: cmake_example
        .. autosummary::
           :toctree: _generate
           add
           subtract
    )pbdoc";

	m.def("Fibo", &Fibo, R"pbdoc(
        Add two numbers
        Some other explanation about the add function.
    )pbdoc");

	m.def("saveFibo", &saveFibo, R"pbdoc(
        Add two numbers
        Some other explanation about the add function.
    )pbdoc");

	m.def(
		"subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers
        Some other explanation about the subtract function.
    )pbdoc");

	m.def("MCMC", &MCMC, R"pbdoc(
		MCMC function
	)pbdoc");
	m.attr("s") = Unit::s;
	m.attr("day") = Unit::day;
	m.attr("yr") = Unit::yr;
	m.attr("__version__") = VERSION;
}
