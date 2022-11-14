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
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include <getopt.h>

#include <fmt/core.h>
#include <gsl/gsl_math.h>

#include "IO.h"
#include "Metric.h"
#include "Object.h"
#include "Unit.h"
#include "Utility.h"
#include "View.h"

using namespace std;
using namespace SBody;

void InterruptHandler(int signum) {
	fmt::print("\r\033[KSBody exiting...\n");
	indicators::show_console_cursor(true);
	exit(signum);
}

void Help(double mass, double spin, double NUT, double tFinal, size_t tStepNumber, double TCal, int metric, int PN, int ray, double absAcc, double relAcc, string store_format) {
	fmt::print("\033[1mSBODY\033[0m ({})\n", VERSION);
	fmt::print("  A relativistic ray tracing program.\n");
	fmt::print("\n\033[1mOPTIONS\033[0m\n");
	fmt::print("  -m --mass   [f]: mass of the source [{}] (M_sun)\n", mass);
	fmt::print("  -s --spin   [f]: spin of the source [{}]\n", spin);
	fmt::print("  -l --NUT    [f]: NUT charge of the source [{}]\n", NUT);
	fmt::print("  -A --sm     [f]: semimajor axis of the star (r_g)\n");
	fmt::print("  -E --ec     [f]: eccentricity of the star\n");
	fmt::print("  -I --in     [f]: inclination of the star (deg)\n");
	fmt::print("  -o --pe     [f]: position angle of periapsis of the star (deg)\n");
	fmt::print("  -O --an     [f]: ascending node of the star (deg)\n");
	fmt::print("  -t --time   [f]: time limit [{}] (s)\n", tFinal);
	fmt::print("  -k --rec    [i]: record number [{}] (s)\n", tStepNumber);
	fmt::print("  -c --tcal   [f]: time limit of calculation [{}] (s)\n", TCal);
	fmt::print("  -n --metric [i]: Newton (0)/Schwarzschild (1)/Kerr (2)/KerrTaubNUT (3) [{}]\n", metric);
	fmt::print("  -P --PN     [i]: PN1 (1) + PN2 (2) + PN2.5 (4) + PN3 (8) + PN3.5 (16) [{}]\n", PN);
	fmt::print("  -R --ray    [i]: ray tracing, view (1) + shadow (4) + camera (2) [{}]\n", ray);
	fmt::print("  -a --abs    [f]: absolute accuracy [{}]\n", absAcc);
	fmt::print("  -r --rel    [f]: relative accuracy [{}]\n", relAcc);
	fmt::print("  -i --inc    [f]: inclination of the BH\n");
	fmt::print("  -e --eps    [f]: epsilon of the BH\n");
	fmt::print("  -f --format [s]: storage format [{}]\n", store_format);
	fmt::print("  -L --light     : light instead of particle\n");
	fmt::print("  -b --bar       : show progress bar\n");
	fmt::print("  -h --help      : this help information\n");
	fmt::print("\n\033[1mEXAMPLES\033[0m\n  $ SBody\n");
	fmt::print("\n\033[1mSUPPORT\033[0m\n  github.com/Tai-Zhou/SBody\n");
}

int main(int argc, char *argv[]) {
	// Application Entry
	auto TStart = chrono::steady_clock::now();
	double mass = 4.15e6, spin = 0., charge = 0., NUT = 0.;
	double tFinal = 2000.;
	size_t tStepNumber = 10000UL;
	double TCal = 36000000;
	size_t metric = 1;
	size_t Hamiltonian = 0;
	size_t PN = 1;
	size_t ray = 1;
	int rest_mass = 1;
	ProgressBar::display_ = true;
	string store_format = "NumPy";
	double inc = M_PI * 0. / 180., eps = 0.;
	double a = 8.3, e = 0., inclination = M_PI * 162. / 180., periapsis = M_PI * 198.9 / 180. + M_PI_2, ascending_node = M_PI * 25.1 / 180., true_anomaly = M_PI_2; // phi=[2.73633242 3.92974873 3.32166381 3.2093593 3.67372211 5.18824159 | 3.19861806 2.63708292 3.05259405]
	unique_ptr<View> viewPtr;
	unique_ptr<thread> shadowPtr;
	unique_ptr<Camera> cameraPtr;
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
			ascending_node = atof(optarg) * M_PI / 180.;
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
				Help(mass, spin, NUT, tFinal, tStepNumber, TCal, metric, PN, ray, absolute_accuracy, relative_accuracy, store_format);
			break;
		case 'P':
			PN = atoi(optarg);
			break;
		case 'R':
			ray = atoi(optarg);
			break;
		case 'a':
			absolute_accuracy = atof(optarg);
			break;
		case 'r':
			relative_accuracy = atof(optarg);
			break;
		case 'i':
			inc = atof(optarg);
			break;
		case 'e':
			eps = atof(optarg);
			break;
		case 'f':
			store_format = atoi(optarg);
			break;
		case 'L':
			rest_mass = 0;
			break;
		case 'b':
			ProgressBar::display_ = true;
			break;
		default:
			Help(mass, spin, NUT, tFinal, tStepNumber, TCal, metric, PN, ray, absolute_accuracy, relative_accuracy, store_format);
			exit(opt == 'h' ? 0 : 1);
		}
	Unit::Initialize(mass);
	tFinal *= Unit::s;
	double t = 0, tStep = 0, tRec = tFinal / tStepNumber;
	if (metric == 0)
		ray = 0;
	shared_ptr<Metric> main_metric = make_shared<Schwarzschild>(T);
	// metric::setMetric(metric, PN, mass, spin, charge, NUT);
	string strFormat = fmt::format(" ({:.1f},{:.1f},{:.1f})[{:f},{:f}]", spin, charge, NUT, inc, eps);
	if (ProgressBar::display_) {
		indicators::show_console_cursor(false);
		ProgressBar::bars_[0].set_option(indicators::option::ForegroundColor{indicators::Color(metric)});
		ProgressBar::bars_[0].set_option(indicators::option::PrefixText{string("?") + strFormat});
	}
	if (ray & 5) {
		viewPtr = make_unique<View>(make_unique<Schwarzschild>(HAMILTONIAN), 8180. * Unit::pc, inc, string("view") + strFormat);
		if (ray & 4)
			shadowPtr = make_unique<thread>(&View::Shadow, viewPtr.get());
	}
	if (ray & 10) {
		cameraPtr = make_unique<Camera>(make_unique<Schwarzschild>(HAMILTONIAN), 1000, 5e-2, 1.e3, inc, string("camera") + strFormat);
		if (ray & 8)
			lensPtr = make_unique<thread>(&Camera::Lens, cameraPtr.get());
	}
	// Integrator integ(metric::function, metric::jacobian, metric != 0);
	Star star_0(main_metric, Unit::R_sun, 0);
	// star_0.InitializeKeplerian(6.5, 0., inclination, periapsis, ascending_node, true_anomaly, inc, eps);
	star_0.InitializeGeodesic(a, inclination, periapsis, ascending_node, -0.16707659553531468, 0.3822615764261866, inc, eps);
	Integrator &&integrator = main_metric->GetIntegrator(metric != 0);
	NumPy rec(main_metric->Name() + strFormat, {12});
	vector<double> temp(12);
	int status = 0;
	double h = 1., stepPercent = 100. / tStepNumber;
	for (int i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		while (status <= 0 && t < tStep)
			status = star_0.IntegratorApply(&t, tStep, &h);
		if (status > 0)
			fmt::print(stderr, "\033[101m[ERR]\033[0m main status = {}\n", status);
		if (metric == 0)
			star_0.Position(temp.data());
		else {
			if (Hamiltonian) {
				star_0.Position(temp.data());
				main_metric->HamiltonianToLagrangian(temp.data()); // TODO: might need s2c()
			} else {
				star_0.Position(temp.data());
				SphericalToCartesian(temp.data());
			}
		}
		temp[8] = t / Unit::s;
		temp[9] = star_0.Energy();
		temp[10] = star_0.AngularMomentum();
		temp[11] = star_0.CarterConstant();
		if (ray & 1)
			viewPtr->TraceStar(star_0, i);
		if (ray & 2)
			cameraPtr->TraceStar();
		rec.Save(temp);
		if (ProgressBar::display_)
			ProgressBar::bars_[0].set_progress(i * stepPercent);
		if (chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - TStart).count() >= TCal)
			break;
	}
	if (ProgressBar::display_) {
		ProgressBar::SetComplete(0, string("!") + strFormat);
		indicators::show_console_cursor(true);
	}
	if (ray & 2)
		cameraPtr->Save();
	if (ray & 4)
		shadowPtr->join();
	if (ray & 8)
		lensPtr->join();
	return 0;
}
