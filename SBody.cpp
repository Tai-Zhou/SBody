// Main program

#include "SBody.h"

#include <chrono>
#include <cmath>
#include <csignal>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include <getopt.h>

#include <gsl/gsl_odeiv2.h>
#include <indicators/block_progress_bar.hpp>
#include <indicators/cursor_control.hpp>

#include "Constant.h"
#include "IO.h"
#include "Kerr.h"
#include "Newton.h"
#include "PostNewtonian.h"
#include "Schwarzschild.h"
#include "Unit.h"
#include "Utility.h"

using namespace std;

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
	cout << "  -n --NPSK   [i]: Newton (0)/PN (1)/Schwarzschild (2)/Kerr (3) [" << NPSK << "]" << endl;
	cout << "  -p --PN     [i]: PN1 (1) + PN2 (2) + PN2.5 (4) + PN3 (8) + PN3.5 (16) [" << PN << "]" << endl;
	cout << "  -a --abs    [f]: absolute accuracy [" << absAcc << "]" << endl;
	cout << "  -r --rel    [f]: relative accuracy [" << relAcc << "]" << endl;
	cout << "  -f --format [s]: storage format [" << storeFormat << "]" << endl;
	cout << "  -l --light     : light instead of particle" << endl;
	cout << "  -b --proBar    : show progress bar" << endl;
	cout << "  -h --help      : this help information" << endl;
	cout << endl;
	cout << "Support:" << endl;
	cout << "  github.com/Tai-Zhou/SBody" << endl;
}

int main(int argc, char *argv[]) {
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
	mass = 1;
	spin = 0;
	absAcc = 1e-15;
	relAcc = 1e-15;
	tFinal = 1e8;
	tRec = tFinal / 100000;
	tCal = 3600;
	NPSK = 1;
	PN = 1;
	PL = 1;
	progressBar = 1;
	storeFormat = "NumPy";
	double t = 0, tStep = 0;
	const char *optShort = "m:s:t:o:c:n:p:a:r:f:lbh";
	static struct option optLong[] = {
		{"mass", required_argument, NULL, 'm'},
		{"spin", required_argument, NULL, 's'},
		{"time", required_argument, NULL, 't'},
		{"rec", required_argument, NULL, 'o'},
		{"tcal", required_argument, NULL, 'c'},
		{"NPSK", required_argument, NULL, 'n'},
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
			NPSK = atoi(optarg);
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
			return 1;
		}
	tFinal *= constant::s;
	tRec *= constant::s;
	if (progressBar)
		indicators::show_console_cursor(false);
	source params(mass, spin);
	postnewtonian::PN = PN;
	const gsl_odeiv2_step_type *ode_type = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_step *ode_step;
	gsl_odeiv2_control *ode_control;
	gsl_odeiv2_evolve *ode_evolve;
	gsl_odeiv2_system ode_system;
	double x[8], y[8], z[8];
	ode_control = gsl_odeiv2_control_y_new(absAcc, relAcc);
	ode_system.params = &params;
	if (NPSK < 2) {
		y[0] = constant::AU / 100;
		y[1] = 0;
		y[2] = 0;
		y[3] = -sqrt(100 * mass / constant::AU) * cos(constant::pi / 10);
		y[4] = sqrt(100 * mass / constant::AU) * sin(constant::pi / 10);
		y[5] = 0;
	}
	else {
		x[0] = 0;
		x[2] = 0;
		x[3] = 0;
		x[7] = 0;
		if (PL) {
			x[1] = constant::AU / 100;
			x[5] = -sqrt(100 * mass / constant::AU) * cos(constant::pi / 10);
			x[6] = sqrt(100 * mass / constant::AU) * sin(constant::pi / 10);
		}
		else {
			x[1] = mass * 442990;
			x[5] = 0;
			x[6] = 1;
		}
	}
	if (NPSK == 0) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, newton::dimension);
		ode_evolve = gsl_odeiv2_evolve_alloc(newton::dimension);
		ode_system.function = newton::function;
		ode_system.jacobian = newton::jacobian;
		ode_system.dimension = newton::dimension;
		output = "newton";
	}
	else if (NPSK == 1) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, postnewtonian::dimension);
		ode_evolve = gsl_odeiv2_evolve_alloc(postnewtonian::dimension);
		ode_system.function = postnewtonian::function;
		ode_system.jacobian = postnewtonian::jacobian;
		ode_system.dimension = postnewtonian::dimension;
		output = "postnewtonian";
	}
	else if (NPSK == 2) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, schwarzschild::dimension);
		ode_evolve = gsl_odeiv2_evolve_alloc(schwarzschild::dimension);
		ode_system.function = schwarzschild::function;
		ode_system.jacobian = schwarzschild::jacobian;
		ode_system.dimension = schwarzschild::dimension;
		if (PL) {
			schwarzschild::c2s(x, y);
			schwarzschild::particle::normalization(y, &params);
		}
		else {
			schwarzschild::c2s(x, y);
			schwarzschild::light::normalization(y, &params);
		}
		output = "schwarzschild";
	}
	else if (NPSK == 3) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, kerr::dimension);
		ode_evolve = gsl_odeiv2_evolve_alloc(kerr::dimension);
		ode_system.function = kerr::function;
		ode_system.jacobian = kerr::jacobian;
		ode_system.dimension = kerr::dimension;
		if (PL) {
			kerr::c2s(x, y);
			kerr::particle::normalization(y, &params);
		}
		else {
			kerr::c2s(x, y);
			kerr::light::normalization(y, &params);
		}
		output = "kerr";
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
		if (NPSK == 0) {
			for (int i = 0; i < 6; ++i)
				temp[i] = y[i];
			temp[6] = t / constant::s;
			temp[7] = newton::energy(y, &params);
			temp[8] = newton::angularMomentum(y, &params);
		}
		else if (NPSK == 1) {
			for (int i = 0; i < 6; ++i)
				temp[i] = y[i];
			temp[6] = t / constant::s;
			temp[7] = postnewtonian::energy(y, &params);
			temp[8] = postnewtonian::angularMomentum(y, &params);
		}
		else if (NPSK == 2) {
			schwarzschild::s2c(y, &temp[0]);
			temp[8] = t / constant::s;
			temp[9] = schwarzschild::energy(y, &params);
			temp[10] = schwarzschild::angularMomentum(y, &params);
		}
		else if (NPSK == 3) {
			kerr::s2c(y, &temp[0]);
			temp[8] = t / constant::s;
			temp[9] = kerr::energy(y, &params);
			temp[10] = kerr::angularMomentum(y, &params);
		}
		rec.push_back(temp);
	}
	if (progressBar) {
		bar.set_option(indicators::option::PrefixText{"Complete!"});
		bar.mark_as_completed();
		indicators::show_console_cursor(true);
	}
	IO::NumPySave(rec, output);
	return 0;
}
