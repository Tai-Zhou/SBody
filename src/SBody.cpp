// Main program

#include "SBody.h"

#include <array>
#include <chrono>
#include <cmath>
#include <csignal>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include <getopt.h>

#include "Constant.h"
#include "IO.h"
#include "Metric.h"
#include "Object.h"
#include "Utility.h"
#include "View.h"

using namespace std;
using namespace SBody;

void interruptHandler(int signum) {
	indicators::show_console_cursor(true);
	exit(signum);
}

void help() {
	cout << "SBody (" << VERSION << ")" << endl;
	cout << endl;
	cout << "Options:" << endl;
	cout << "  -m --mass   [f]: mass of the source [" << mass << "] (M_sun)" << endl;
	cout << "  -s --spin   [f]: spin of the source [" << spin << "]" << endl;
	cout << "  -l --NUT    [f]: NUT charge of the source [" << NUT << "]" << endl;
	cout << "  -A --sm     [f]: semimajor axis of the star (mpc)" << endl;
	cout << "  -E --ec     [f]: eccentricity of the star" << endl;
	cout << "  -I --in     [f]: inclination of the star (deg)" << endl;
	cout << "  -o --pe     [f]: position angle of periapsis of the star (deg)" << endl;
	cout << "  -O --an     [f]: ascending node of the star (deg)" << endl;
	cout << "  -t --time   [f]: time limit [" << tFinal << "] (s)" << endl;
	cout << "  -k --rec    [f]: record interval [" << tRec << "] (s)" << endl;
	cout << "  -c --tcal   [f]: time limit of calculation [" << TCal << "] (s)" << endl;
	cout << "  -n --metric [i]: Newton (0)/Schwarzschild (1)/Kerr (2)/KerrTaubNUT (3) [" << metric << "]" << endl;
	cout << "  -P --PN     [i]: PN1 (1) + PN2 (2) + PN2.5 (4) + PN3 (8) + PN3.5 (16) [" << PN << "]" << endl;
	cout << "  -R --ray    [i]: ray tracing, view (1) + shadow (4) + camera (2) [" << ray << "]" << endl;
	cout << "  -a --abs    [f]: absolute accuracy [" << absAcc << "]" << endl;
	cout << "  -r --rel    [f]: relative accuracy [" << relAcc << "]" << endl;
	cout << "  -i --inc    [f]: inclination of the BH" << endl;
	cout << "  -e --eps    [f]: epsilon of the BH" << endl;
	cout << "  -f --format [s]: storage format [" << storeFormat << "]" << endl; // TODO: use int instead of str
	cout << "  -h --hamilton  : use hamilonian instead of geodesic" << endl;
	cout << "  -L --light     : light instead of particle" << endl;
	cout << "  -b --bar       : show progress bar" << endl;
	cout << "  -h --help      : this help information" << endl;
	cout << endl;
	cout << "Support:" << endl;
	cout << "  github.com/Tai-Zhou/SBody" << endl;
	exit(1);
}

int main(int argc, char *argv[]) {
	auto TStart = chrono::steady_clock::now();
	signal(SIGINT, interruptHandler);
	mass = 4.261e6;
	spin = 0.;
	charge = 0.;
	NUT = 0.;
	// tFinal = 9.35062553 * 3.15576e7;
	tRec = 1e-4 * tFinal;
	TCal = 36000000;
	metric = 0;
	Hamiltonian = 0;
	PN = 1;
	ray = 1;
	int restMass = 1;
	IO::displayProgressBar = 1;
	storeFormat = "NumPy";
	double t = 0, tStep = 0;
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
	double a = 5. * Constant::mpc, e = 0.884649, inclination = M_PI * 134.567 / 180., periapsis = M_PI * 66.263 / 180., ascendingNode = M_PI * 228.171 / 180., trueAnomaly = M_PI; // phi=[2.73633242 3.92974873 3.32166381 3.2093593 3.67372211 5.18824159 | 3.19861806 2.63708292 3.05259405]
	unique_ptr<view> viewPtr;
	unique_ptr<thread> shadowPtr;
	unique_ptr<camera> cameraPtr;
	unique_ptr<thread> lensPtr;
	const char *optShort = "m:s:l:A:E:I:o:O:t:k:c:n:P:R:a:r:i:e:f:Lbh";
	const struct option optLong[] = {
		{"mass", required_argument, NULL, 'm'},
		{"spin", required_argument, NULL, 's'},
		{"NUT", required_argument, NULL, 'l'},
		{"sm", required_argument, NULL, 'A'},
		{"ec", required_argument, NULL, 'E'},
		{"in", required_argument, NULL, 'I'},
		{"pe", required_argument, NULL, 'o'},
		{"an", required_argument, NULL, 'O'},
		{"time", required_argument, NULL, 't'},
		{"rec", required_argument, NULL, 'k'},
		{"tcal", required_argument, NULL, 'c'},
		{"metric", required_argument, NULL, 'n'},
		{"PN", required_argument, NULL, 'P'},
		{"ray", required_argument, NULL, 'R'},
		{"abs", required_argument, NULL, 'a'},
		{"rel", required_argument, NULL, 'r'},
		{"inc", required_argument, NULL, 'i'},
		{"eps", required_argument, NULL, 'e'},
		{"format", required_argument, NULL, 'f'},
		{"light", no_argument, NULL, 'L'},
		{"bar", no_argument, NULL, 'b'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}};
	int opt;
	while ((opt = getopt_long(argc, argv, optShort, optLong, NULL)) != -1)
		switch (opt) {
		case 'm':
			mass = atof(optarg);
			break;
		case 's':
			spin = atof(optarg);
			break;
		case 'l':
			NUT = atof(optarg);
			break;
		case 'A':
			a = atof(optarg) * Constant::mpc;
			break;
		case 'E':
			e = atof(optarg);
			break;
		case 'I':
			inclination = atof(optarg) * M_PI / 180.;
			break;
		case 'o':
			periapsis = atof(optarg) * M_PI / 180.;
			break;
		case 'O':
			ascendingNode = atof(optarg) * M_PI / 180.;
			break;
		case 't':
			tFinal = atof(optarg);
			break;
		case 'k':
			tRec = atof(optarg);
			break;
		case 'c':
			TCal = atof(optarg);
			break;
		case 'n':
			metric = atoi(optarg);
			if (metric > 3)
				help();
			break;
		case 'P':
			PN = atoi(optarg);
			break;
		case 'R':
			ray = atoi(optarg);
			break;
		case 'a':
			absAcc = atof(optarg);
			break;
		case 'r':
			relAcc = atof(optarg);
			break;
		case 'i':
			inc = atof(optarg);
			break;
		case 'e':
			eps = atof(optarg);
			break;
		case 'f':
			storeFormat = atoi(optarg);
			break;
		case 'L':
			restMass = 0;
			break;
		case 'b':
			IO::displayProgressBar = 1;
			break;
		default:
			help();
		}
	tFinal *= Constant::s;
	tRec *= Constant::s;
	if (metric == 0)
		ray = 0;
	Metric::setMetric(metric, PN, mass, spin, charge, NUT);
	char strFormat[1024]; // TODO: waiting for C++20
	snprintf(strFormat, 1024, " (%.1f,%.1f,%.1f)[%f,%f]", spin, charge, NUT, inc, eps);
	if (IO::displayProgressBar) {
		indicators::show_console_cursor(false);
		IO::progressBars[0].set_option(indicators::option::ForegroundColor{indicators::Color(metric)});
		IO::progressBars[0].set_option(indicators::option::PrefixText{string("?") + strFormat});
	}
	if (ray & 5) {
		viewPtr = make_unique<view>(8246.7 * Constant::pc, inc, string("view") + strFormat);
		if (ray & 4)
			shadowPtr = make_unique<thread>(&view::shadow, viewPtr.get(), 100);
	}
	if (ray & 10) {
		cameraPtr = make_unique<camera>(1000, 5e-2, mass * 1.e3, inc, string("camera") + strFormat);
		if (ray & 8)
			lensPtr = make_unique<thread>(&camera::lens, cameraPtr.get());
	}
	double x[8], y[8];
	if (metric == 0) {
		y[0] = 0;
		y[4] = 1;
		double r = a * (1 - e * e) / (1 + e * cos(trueAnomaly));
		double tp1 = -r * cos(periapsis + trueAnomaly), tp2 = -r * sin(periapsis + trueAnomaly) * cos(inclination);
		double xp1 = tp1 * cos(ascendingNode) - tp2 * sin(ascendingNode), xp2 = tp2 * cos(ascendingNode) + tp1 * sin(ascendingNode), xp3 = -r * sin(periapsis + trueAnomaly) * sin(inclination);
		y[1] = (xp1 * cos(eps) + xp2 * sin(eps)) * cos(inc) + xp3 * sin(inc);
		y[2] = xp2 * cos(eps) - xp1 * sin(eps);
		y[3] = xp3 * cos(inc) - (xp1 * cos(eps) + xp2 * sin(eps)) * sin(inc);
		double vtheta = sqrt((1 - e * e) * mass * a) / r, vr = sign(M_PI - mod2Pi(trueAnomaly)) * sqrt(max(0., 2. * mass / r - mass / a - vtheta * vtheta));
		// vtheta = 0.3822615764261866;
		// vr = -0.16707659553531468;
		double tp5 = vtheta * sin(periapsis + trueAnomaly) - vr * cos(periapsis + trueAnomaly), tp6 = -(vtheta * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * cos(inclination);
		double xp5 = tp5 * cos(ascendingNode) - tp6 * sin(ascendingNode), xp6 = tp5 * sin(ascendingNode) + tp6 * cos(ascendingNode), xp7 = -(vtheta * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * sin(inclination);
		y[5] = (xp5 * cos(eps) + xp6 * sin(eps)) * cos(inc) + xp7 * sin(inc);
		y[6] = xp6 * cos(eps) - xp5 * sin(eps);
		y[7] = xp7 * cos(inc) - (xp5 * cos(eps) + xp6 * sin(eps)) * sin(inc);
	}
	else {
		x[0] = 0;
		if (restMass) {
			double r = a * (1 - e * e) / (1 + e * cos(trueAnomaly));
			double tp1 = -r * cos(periapsis + trueAnomaly), tp2 = -r * sin(periapsis + trueAnomaly) * cos(inclination);
			double xp1 = tp1 * cos(ascendingNode) - tp2 * sin(ascendingNode), xp2 = tp2 * cos(ascendingNode) + tp1 * sin(ascendingNode), xp3 = -r * sin(periapsis + trueAnomaly) * sin(inclination);
			x[1] = (xp1 * cos(eps) + xp2 * sin(eps)) * cos(inc) + xp3 * sin(inc);
			x[2] = xp2 * cos(eps) - xp1 * sin(eps);
			x[3] = xp3 * cos(inc) - (xp1 * cos(eps) + xp2 * sin(eps)) * sin(inc);
			double vtheta = sqrt((1 - e * e) * mass * a) / r, vr = sign(M_PI - mod2Pi(trueAnomaly)) * sqrt(max(0., 2. * mass / r - mass / a - vtheta * vtheta));
			// vtheta = 0.3822615764261866;
			// vr = -0.16707659553531468;
			double tp5 = vtheta * sin(periapsis + trueAnomaly) - vr * cos(periapsis + trueAnomaly), tp6 = -(vtheta * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * cos(inclination);
			double xp5 = tp5 * cos(ascendingNode) - tp6 * sin(ascendingNode), xp6 = tp5 * sin(ascendingNode) + tp6 * cos(ascendingNode), xp7 = -(vtheta * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * sin(inclination);
			x[5] = (xp5 * cos(eps) + xp6 * sin(eps)) * cos(inc) + xp7 * sin(inc);
			x[6] = xp6 * cos(eps) - xp5 * sin(eps);
			x[7] = xp7 * cos(inc) - (xp5 * cos(eps) + xp6 * sin(eps)) * sin(inc);
		}
		else {
			x[1] = mass * 1000;
			x[5] = -1.;
			x[6] = 0;
			x[7] = 0;
		}
		Metric::c2s(x, y);
		/*y[0] = 0;
		y[1] = 3.588 * Constant::mpc;
		y[2] = M_PI / 3;
		y[3] = M_PI;
		y[5] = 0;
		y[6] = 0;
		y[7] = sqrt(4e6 / (gsl_pow_3(y[1]) * 0.75)) * 1.000000578443; //[1.0000005784, 1.00000057845]*/
		if (restMass)
			Metric::particleNormalization(y);
		else
			Metric::lightNormalization(y, 1.);
		if (Hamiltonian)
			Metric::qdq2qp(y);
	}
	integrator integ(Metric::function, Metric::jacobian, metric != 0, -1.);
	Object::star star_0(Constant::R_sun, y, 0);
	Object::objectList.push_back(&star_0);
	IO::NumPy rec(Metric::name + strFormat, {12});
	vector<double> temp(12);
	int status = 0, rayNO = 0, TUse, TLastUse = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - TStart).count();
	double tFinal_1 = 1 / tFinal;
	double tPoint = tList[0] * Constant::yr;
	t = 2010.3566473381456 * Constant::yr;
	integ.resetHstart(-1.);
	while (t > tPoint)
		status = integ.apply(&t, tPoint, star_0.pos);
	integ.resetHstart(1.);
	for (double tPoint : tList) {
		// tStep = min(tStep + tRec, tFinal);
		tPoint = tPoint * Constant::yr;
		while (status <= 0 && t < tPoint)
			status = integ.apply(&t, tPoint, star_0.pos);
		if (status > 0)
			cerr << "[!] main status = " << status << endl;
		if (metric == 0)
			copy(star_0.pos, star_0.pos + 8, temp.begin());
		else {
			if (Hamiltonian) {
				copy(star_0.pos, star_0.pos + 8, temp.begin());
				Metric::qp2qdq(temp.data()); // TODO: need s2c()
			}
			else
				Metric::s2c(star_0.pos, temp.data());
		}
		temp[8] = t / Constant::s;
		temp[9] = Metric::energy(star_0.pos);
		temp[10] = Metric::angularMomentum(star_0.pos);
		temp[11] = Metric::carter(star_0.pos, 1.);
		if (ray & 1)
			viewPtr->traceBack(star_0, rayNO++);
		if (ray & 2)
			cameraPtr->traceBack();
		rec.save(temp);
		TUse = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - TStart).count();
		if (IO::displayProgressBar) {
			double percent = 100. * t * tFinal_1;
			if (floor(percent) > floor(IO::progressBars[0].current()) || TUse >= TLastUse + 500) {
				TLastUse = TUse;
				IO::progressBars[0].set_progress(percent);
			}
		}
		if (TUse >= TCal)
			break;
	}
	if (IO::displayProgressBar)
		IO::progressBarComplete(0, string("!") + strFormat);
	if (ray & 2)
		cameraPtr->save();
	if (ray & 4)
		shadowPtr->join();
	if (ray & 8)
		lensPtr->join();
	return 0;
}
