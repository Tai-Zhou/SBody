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

#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <indicators/block_progress_bar.hpp>
#include <indicators/cursor_control.hpp>

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
	cout << "  -t --time   [f]: time limit [" << tFinal << "] (s)" << endl;
	cout << "  -o --rec    [f]: record interval [" << tRec << "] (s)" << endl;
	cout << "  -c --tcal   [f]: time limit of calculation [" << tCal << "] (s)" << endl;
	cout << "  -n --NSK    [i]: Newton (0)/Schwarzschild (1)/Kerr (2)/KerrH (3) [" << NSK << "]" << endl;
	cout << "  -P --PN     [i]: PN1 (1) + PN2 (2) + PN2.5 (4) + PN3 (8) + PN3.5 (16) [" << PN << "]" << endl;
	cout << "  -R --ray    [i]: ray tracing, view (1) + camera (2) [" << ray << "]" << endl;
	cout << "  -a --abs    [f]: absolute accuracy [" << absAcc << "]" << endl;
	cout << "  -r --rel    [f]: relative accuracy [" << relAcc << "]" << endl;
	cout << "  -f --format [s]: storage format [" << storeFormat << "]" << endl; //TODO: use int instead of str
	cout << "  -h --hamilton  : use hamilonian instead of geodesic" << endl;
	cout << "  -l --light     : light instead of particle" << endl;
	cout << "  -b --bar       : show progress bar" << endl;
	cout << "  -h --help      : this help information" << endl;
	cout << endl;
	cout << "Support:" << endl;
	cout << "  github.com/Tai-Zhou/SBody" << endl;
	exit(1);
}

int main(int argc, char *argv[]) {
	auto start = chrono::steady_clock::now();
	indicators::BlockProgressBar bar{
		indicators::option::BarWidth{80},
		indicators::option::Start{"["},
		indicators::option::End{"]"},
		indicators::option::PrefixText{"Calculating..."},
		indicators::option::ForegroundColor{indicators::Color::cyan},
		indicators::option::ShowElapsedTime{true},
		indicators::option::ShowRemainingTime{true},
		indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}};
	signal(SIGINT, interruptHandler);
	double h = 1e-3;
	mass = 4e6;
	spin = 0;
	tFinal = 1e4;
	tRec = 1e-4 * tFinal;
	tCal = 3600;
	NSK = 1;
	PN = 1;
	ray = 1;
	int PL = 1;
	int progressBar = 1;
	storeFormat = "NumPy";
	double t = 0, tStep = 0;
	const char *optShort = "m:s:t:o:c:n:P:R:a:r:f:lbh";
	const struct option optLong[] = {
		{"mass", required_argument, NULL, 'm'},
		{"spin", required_argument, NULL, 's'},
		{"time", required_argument, NULL, 't'},
		{"rec", required_argument, NULL, 'o'},
		{"tcal", required_argument, NULL, 'c'},
		{"NSK", required_argument, NULL, 'n'},
		{"PN", required_argument, NULL, 'P'},
		{"ray", required_argument, NULL, 'R'},
		{"abs", required_argument, NULL, 'a'},
		{"rel", required_argument, NULL, 'r'},
		{"format", required_argument, NULL, 'f'},
		{"light", no_argument, NULL, 'l'},
		{"ray", no_argument, NULL, 'R'},
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
		case 't':
			tFinal = atof(optarg);
			break;
		case 'o':
			tRec = atof(optarg);
			break;
		case 'c':
			tCal = atof(optarg);
			break;
		case 'n':
			NSK = atoi(optarg);
			if (NSK > 3)
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
		case 'f':
			storeFormat = atoi(optarg);
			break;
		case 'l':
			PL = 0;
			break;
		case 'b':
			progressBar = 1;
			break;
		default:
			help();
		}
	tFinal *= Constant::s;
	tRec *= Constant::s;
	Metric::setMetric(NSK, mass, spin);
	Metric::Newton::PN = PN;
	if (progressBar)
		indicators::show_console_cursor(false);
	double x[8], y[8];
	if (NSK == 0) {
		y[0] = Constant::AU / 100;
		y[1] = 0;
		y[2] = 0;
		y[3] = -sqrt(100 * mass / Constant::AU) * cos(M_PI / 10);
		y[4] = sqrt(100 * mass / Constant::AU) * sin(M_PI / 10);
		y[5] = 0;
	}
	else {
		x[0] = 0;
		x[2] = 0;
		x[3] = 0;
		if (PL) {
			x[1] = 10 * mass;
			x[5] = 0;
			x[6] = 1.05 * sqrt(0.1);
			x[7] = 0;
		}
		else {
			x[1] = mass * 442990;
			x[5] = 0;
			x[6] = 1;
			x[7] = 0;
		}
		Metric::c2s(x, y);
		if (PL)
			Metric::particleNormalization(y);
		else
			Metric::lightNormalization(y, 1.);
		if (NSK == 3)
			Metric::KerrH::qdq2qp(y);
	}
	integrator integ(Metric::function, Metric::jacobian, NSK == 0);
	int status = 0;
	Object::star star_0(Constant::R_sun, y, 0);
	Object::objectList.push_back(&star_0);
	view vie(mass * 1e3, M_PI_4, tFinal, 3000);
	camera cam(100, 4e-2, mass * 1000., M_PI_4, tFinal, 3000);
	if (ray & 2)
		cam.initialize();
	vector<vector<double>> rec;
	vector<double> temp(12);
	while (tStep < tFinal) {
		tStep = min(tStep + tRec, tFinal);
		while (status == 0 && t < tStep)
			status = integ.apply(&t, tStep, &h, star_0.pos);
		if (NSK == 0)
			temp = {star_0.pos[0], star_0.pos[1], star_0.pos[2], star_0.pos[3], star_0.pos[4], star_0.pos[5], star_0.pos[6], star_0.pos[7]};
		else {
			if (NSK == 3)
				Metric::KerrH::qp2qdq(y); //FIXME: must not be replaced!
			Metric::s2c(star_0.pos, &temp[0]);
		}
		temp[8] = t / Constant::s;
		temp[9] = Metric::energy(star_0.pos);
		temp[10] = Metric::angularMomentum(star_0.pos);
		temp[11] = Metric::carter(star_0.pos);
		if (ray & 1)
			vie.traceBack(star_0);
		if (ray & 2)
			cam.traceBack();
		rec.push_back(temp);
		if (progressBar)
			bar.set_progress(100 * t / tFinal);
		auto tpass = chrono::steady_clock::now() - start; // nano seconds
		if (tpass.count() >= tCal * 1000000000)
			break;
	}
	if (progressBar) {
		bar.set_option(indicators::option::PrefixText{"Complete!"});
		bar.mark_as_completed();
		indicators::show_console_cursor(true);
	}
	if (ray & 1)
		vie.save("view");
	if (ray & 2)
		cam.save("camera");
	IO::NumPy<double> output(Metric::name);
	output.save(rec);
	return 0;
}
