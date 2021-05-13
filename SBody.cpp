// Main program

#include "SBody.h"

#include <array>
#include <chrono>
#include <cmath>
#include <csignal>
#include <iostream>
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
#include "Unit.h"
#include "Utility.h"
#include "View.h"

using namespace std;
using namespace SBody;

void interruptHandler(int signum) {
	if (progressBar)
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
	cout << "  -p --PN     [i]: PN1 (1) + PN2 (2) + PN2.5 (4) + PN3 (8) + PN3.5 (16) [" << PN << "]" << endl;
	cout << "  -a --abs    [f]: absolute accuracy [" << absAcc << "]" << endl;
	cout << "  -r --rel    [f]: relative accuracy [" << relAcc << "]" << endl;
	cout << "  -f --format [s]: storage format [" << storeFormat << "]" << endl;
	cout << "  -h --hamilton  : use hamilonian instead of geodesic" << endl;
	cout << "  -l --light     : light instead of particle" << endl;
	cout << "  -b --proBar    : show progress bar" << endl;
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
	h = 1e-3;
	mass = 4e6;
	spin = 0;
	absAcc = 1e-15;
	relAcc = 1e-15;
	tFinal = 1e0;
	tRec = 1e-5 * tFinal;
	tCal = 3600;
	NSK = 3;
	PN = 1;
	PL = 1;
	progressBar = 0;
	storeFormat = "NumPy";
	double t = 0, tStep = 0;
	const char *optShort = "m:s:t:o:c:n:p:a:r:f:lbh";
	const struct option optLong[] = {
		{"mass", required_argument, NULL, 'm'},
		{"spin", required_argument, NULL, 's'},
		{"time", required_argument, NULL, 't'},
		{"rec", required_argument, NULL, 'o'},
		{"tcal", required_argument, NULL, 'c'},
		{"NSK", required_argument, NULL, 'n'},
		{"PN", required_argument, NULL, 'p'},
		{"abs", required_argument, NULL, 'a'},
		{"rel", required_argument, NULL, 'r'},
		{"format", required_argument, NULL, 'f'},
		{"light", no_argument, NULL, 'l'},
		{"proBar", no_argument, NULL, 'b'},
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
			if (NSK < 0 || NSK > 3)
				help();
			break;
		case 'p':
			PN = atoi(optarg);
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
	if (progressBar)
		indicators::show_console_cursor(false);
	source params(mass, spin * mass);
	Metric::Newton::PN = PN;
	const gsl_odeiv2_step_type *ode_type = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_step *ode_step;
	gsl_odeiv2_control *ode_control;
	gsl_odeiv2_evolve *ode_evolve;
	gsl_odeiv2_system ode_system;
	double x[8], y[8], z[8];
	ode_control = gsl_odeiv2_control_y_new(absAcc, relAcc);
	ode_system.params = &params;
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
			x[6] = sqrt(0.1) * sin(151. / 180. * M_PI);
			x[7] = sqrt(0.1) * cos(151. / 180. * M_PI);
		}
		else {
			x[1] = mass * 442990;
			x[5] = 0;
			x[6] = 1;
			x[7] = 0;
		}
	}
	if (NSK == 0) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, Metric::Newton::dimension);
		ode_evolve = gsl_odeiv2_evolve_alloc(Metric::Newton::dimension);
		ode_system.function = Metric::Newton::function;
		ode_system.jacobian = Metric::Newton::jacobian;
		ode_system.dimension = Metric::Newton::dimension;
		output = "Newton";
	}
	else if (NSK == 1) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, Metric::Schwarzschild::dimension);
		ode_evolve = gsl_odeiv2_evolve_alloc(Metric::Schwarzschild::dimension);
		ode_system.function = Metric::Schwarzschild::function;
		ode_system.jacobian = Metric::Schwarzschild::jacobian;
		ode_system.dimension = Metric::Schwarzschild::dimension;
		Metric::c2s(x, y, Metric::Schwarzschild::dimension);
		if (PL)
			Metric::Schwarzschild::particleNormalization(y, &params);
		else
			Metric::Schwarzschild::lightNormalization(y, &params);
		output = "Schwarzschild";
	}
	else if (NSK == 2) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, Metric::Kerr::dimension);
		ode_evolve = gsl_odeiv2_evolve_alloc(Metric::Kerr::dimension);
		ode_system.function = Metric::Kerr::function;
		ode_system.jacobian = Metric::Kerr::jacobian;
		ode_system.dimension = Metric::Kerr::dimension;
		Metric::c2s(x, y, Metric::Kerr::dimension);
		if (PL)
			Metric::Kerr::particleNormalization(y, &params);
		else
			Metric::Kerr::lightNormalization(y, &params);
		output = "Kerr";
	}
	else if (NSK == 3) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, Metric::KerrH::dimension);
		ode_evolve = gsl_odeiv2_evolve_alloc(Metric::KerrH::dimension);
		ode_system.function = Metric::KerrH::function;
		ode_system.jacobian = Metric::KerrH::jacobian;
		ode_system.dimension = Metric::KerrH::dimension;
		Metric::c2s(x, z, Metric::KerrH::dimension);
		if (PL)
			Metric::KerrH::particleNormalization(z, &params);
		else
			Metric::KerrH::lightNormalization(z, &params);
		Metric::KerrH::qdq2qp(z, y, &params);
		output = "KerrH";
	}
	int status = 0;
	vector<vector<double>> rec;
	vector<double> temp(ode_system.dimension + 3);
	while (tStep < tFinal) {
		tStep = min(tStep + tRec, tFinal);
		while (status == 0 && t != tStep)
			status = gsl_odeiv2_evolve_apply(ode_evolve, ode_control, ode_step, &ode_system, &t, tStep, &h, y);
		if (progressBar)
			bar.set_progress(100 * t / tFinal);
		if (NSK == 0) {
			for (int i = 0; i < 6; ++i)
				temp[i] = y[i];
			temp[6] = t / Constant::s;
			temp[7] = Metric::Newton::energy(y, &params);
			temp[8] = Metric::Newton::angularMomentum(y, &params);
		}
		else if (NSK == 1) {
			Metric::s2c(y, &temp[0], Metric::Schwarzschild::dimension);
			temp[8] = t / Constant::s;
			temp[9] = Metric::Schwarzschild::energy(y, &params);
			temp[10] = Metric::Schwarzschild::angularMomentum(y, &params);
		}
		else if (NSK == 2) {
			Metric::s2c(y, &temp[0], Metric::Kerr::dimension);
			temp[8] = t / Constant::s;
			temp[9] = Metric::Kerr::energy(y, &params);
			//temp[10] = Metric::Kerr::angularMomentum(y, &params);
			temp[10] = Metric::Kerr::carter(y, &params);
		}
		else if (NSK == 3) {
			Metric::KerrH::qp2qdq(y, x, &params);
			Metric::s2c(x, &temp[0], Metric::KerrH::dimension);
			temp[8] = t / Constant::s;
			temp[9] = Metric::KerrH::energy(y, &params);
			//temp[10] = Metric::KerrH::angularMomentum(y, &params);
			temp[10] = Metric::KerrH::carter(y, &params);
		}
		rec.push_back(temp);
		auto tpass = chrono::steady_clock::now() - start; // nano seconds
		if (tpass.count() >= tCal * 1000000000)
			return 0;
	}
	if (progressBar) {
		bar.set_option(indicators::option::PrefixText{"Complete!"});
		bar.mark_as_completed();
		indicators::show_console_cursor(true);
	}
	IO::NumPySave(rec, output);
	return 0;
}
