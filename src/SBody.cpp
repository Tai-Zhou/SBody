/**
 * @file SBody.cpp
 * @author Tai Zhou
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

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

#include "IO.h"
#include "Metric.h"
#include "Object.h"
#include "Unit.h"
#include "Utility.h"
#include "View.h"

using namespace std;
using namespace SBody;

void interruptHandler(int signum) {
	indicators::show_console_cursor(true);
	exit(signum);
}

void help(double mass, double spin, double NUT, double tFinal, size_t tStepNumber, double TCal, int metric, int PN, int ray, double absAcc, double relAcc, string storeFormat) {
	cout << "SBody (" << VERSION << ")" << endl;
	cout << endl;
	cout << "Options:" << endl;
	cout << "  -m --mass   [f]: mass of the source [" << mass << "] (M_sun)" << endl;
	cout << "  -s --spin   [f]: spin of the source [" << spin << "]" << endl;
	cout << "  -l --NUT    [f]: NUT charge of the source [" << NUT << "]" << endl;
	cout << "  -A --sm     [f]: semimajor axis of the star (r_g)" << endl;
	cout << "  -E --ec     [f]: eccentricity of the star" << endl;
	cout << "  -I --in     [f]: inclination of the star (deg)" << endl;
	cout << "  -o --pe     [f]: position angle of periapsis of the star (deg)" << endl;
	cout << "  -O --an     [f]: ascending node of the star (deg)" << endl;
	cout << "  -t --time   [f]: time limit [" << tFinal << "] (s)" << endl;
	cout << "  -k --rec    [i]: record number [" << tStepNumber << "]" << endl;
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
	// Application Entry
	auto TStart = chrono::steady_clock::now();
	signal(SIGINT, interruptHandler);
	double mass = 4.15e6, spin = 0., charge = 0., NUT = 0.;
	double tFinal = 2000.;
	size_t tStepNumber = 10000UL;
	double TCal = 36000000;
	size_t metric = 1;
	size_t Hamiltonian = 0;
	size_t PN = 1;
	size_t ray = 1;
	int restMass = 1;
	IO::displayProgressBar = 1;
	string storeFormat = "NumPy";
	double inc = M_PI * 0. / 180., eps = 0.;
	double a = 8.3, e = 0., inclination = M_PI * 162. / 180., periapsis = M_PI * 198.9 / 180., ascendingNode = M_PI * 25.1 / 180., trueAnomaly = M_PI_2; // phi=[2.73633242 3.92974873 3.32166381 3.2093593 3.67372211 5.18824159 | 3.19861806 2.63708292 3.05259405]
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
			a = atof(optarg);
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
			tStepNumber = atoi(optarg);
			break;
		case 'c':
			TCal = atof(optarg);
			break;
		case 'n':
			metric = atoi(optarg);
			if (metric > 3)
				help(mass, spin, NUT, tFinal, tStepNumber, TCal, metric, PN, ray, absAcc, relAcc, storeFormat);
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
			help(mass, spin, NUT, tFinal, tStepNumber, TCal, metric, PN, ray, absAcc, relAcc, storeFormat);
		}
	Unit::init(mass);
	tFinal *= Unit::s;
	double t = 0, tStep = 0, tRec = tFinal / tStepNumber;
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
		viewPtr = make_unique<view>(8180. * Unit::pc, inc, string("view") + strFormat);
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
		double vphi = sqrt((1 - e * e) * mass * a) / r, vr = sign(M_PI - mod2Pi(trueAnomaly)) * sqrt(max(0., 2. * mass / r - mass / a - vphi * vphi));
		// vphi = 0.3822615764261866;
		// vr = -0.16707659553531468;
		double tp5 = vphi * sin(periapsis + trueAnomaly) - vr * cos(periapsis + trueAnomaly), tp6 = -(vphi * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * cos(inclination);
		double xp5 = tp5 * cos(ascendingNode) - tp6 * sin(ascendingNode), xp6 = tp5 * sin(ascendingNode) + tp6 * cos(ascendingNode), xp7 = -(vphi * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * sin(inclination);
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
			double vphi = sqrt((1 - e * e) * mass * a) / r, vr = sign(M_PI - mod2Pi(trueAnomaly)) * sqrt(max(0., 2. * mass / r - mass / a - vphi * vphi));
			vphi = 0.3822615764261866; // 0.3840905791842067; // 0.3822615764261866;
			vr = -0.16707659553531468; // 0.16386176496281388; // 0.16707659553531468;
			double tp5 = vphi * sin(periapsis + trueAnomaly) - vr * cos(periapsis + trueAnomaly), tp6 = -(vphi * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * cos(inclination);
			double xp5 = tp5 * cos(ascendingNode) - tp6 * sin(ascendingNode), xp6 = tp5 * sin(ascendingNode) + tp6 * cos(ascendingNode), xp7 = -(vphi * cos(periapsis + trueAnomaly) + vr * sin(periapsis + trueAnomaly)) * sin(inclination);
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
		if (restMass)
			Metric::particleNormalization(y);
		else
			Metric::lightNormalization(y, 1.);
		if (Hamiltonian)
			Metric::qdq2qp(y);
	}
	integrator integ(Metric::function, Metric::jacobian, metric != 0);
	Object::star star_0(Unit::R_sun, y, 0);
	Object::objectList.push_back(&star_0);
	IO::NumPy rec(Metric::name + strFormat, {12});
	vector<double> temp(12);
	int status = 0, TUse, TLastUse = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - TStart).count();
	double h = 1., stepPercent = 100. / tStepNumber;
	for (int i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		while (status <= 0 && t < tStep)
			status = integ.apply(&t, tStep, &h, star_0.pos);
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
		temp[8] = t / Unit::s;
		temp[9] = Metric::energy(star_0.pos);
		temp[10] = Metric::angularMomentum(star_0.pos);
		temp[11] = Metric::carter(star_0.pos, 1.);
		if (ray & 1)
			viewPtr->traceStar(star_0, i);
		if (ray & 2)
			cameraPtr->traceStar();
		rec.save(temp);
		TUse = chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - TStart).count();
		if (IO::displayProgressBar) {
			double percent = i * stepPercent;
			if (percent > IO::progressBars[0].current() + 0.3 || TUse >= TLastUse + 300) {
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
