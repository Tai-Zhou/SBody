// Main program

#include "SBody.h"

#include <chrono>
#include <cmath>
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
#include "Newton.h"
#include "PostNewtonian.h"
#include "Schwarzschild.h"
#include "Unit.h"
#include "Utility.h"

using namespace std;

void help() {
	cout << "SBody (" << VERSION << ")" << endl;
	cout << endl;
	cout << "Options:" << endl;
	cout << "  -m --mass   [f]: mass of the source [" << mass << "] (M_sun)" << endl;
	cout << "  -s --spin   [f]: spin of the source [" << spin << "]" << endl;
	cout << "  -t --time   [f]: time limit [" << tFinal << "] (s)" << endl;
	cout << "  -c --tcal   [f]: time limit of calculation [" << tCal << "] (s)" << endl;
	cout << "  -n --NPSK   [i]: Newton (0)/PN (1)/Schwarzschild (2)/Kerr (3) [" << NPSK << "]" << endl;
	cout << "  -p --PN     [i]: PN1 (1) + PN2 (2) + PN2.5 (4) + PN3 (8) + PN3.5 (16) [" << PN << "]" << endl;
	cout << "  -o --rec    [i]: [" << recRatio << "] integration steps per record" << endl;
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
	h = 1e-3;
	mass = 1;
	spin = 0;
	absAcc = 1e-15;
	relAcc = 1e-15;
	tFinal = 1e9;
	tCal = 3600;
	NPSK = 1;
	PN = 1;
	PL = 1;
	recRatio = 1;
	progressBar = 1;
	storeFormat = "NumPy";
	output = "output";
	double t = 0;
	const char *optShort = "m:s:t:c:n:p:o:a:r:f:lbh";
	static struct option optLong[] = {
		{"mass", required_argument, NULL, 'm'},
		{"spin", required_argument, NULL, 's'},
		{"time", required_argument, NULL, 't'},
		{"tcal", required_argument, NULL, 'c'},
		{"NPSK", required_argument, NULL, 'n'},
		{"PN", required_argument, NULL, 'p'},
		{"rec", required_argument, NULL, 'o'},
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
		case 'c':
			tCal = atof(optarg);
			break;
		case 'n':
			NPSK = atoi(optarg);
			break;
		case 'o':
			recRatio = atoi(optarg);
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
	if (progressBar)
		indicators::show_console_cursor(false);
	source params(mass, spin);
	postnewtonian::PN = 1;
	const gsl_odeiv2_step_type *ode_type = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_step *ode_step;
	gsl_odeiv2_control *ode_control;
	gsl_odeiv2_evolve *ode_evolve;
	gsl_odeiv2_system ode_system;
	double x[8], y[8];
	if (NPSK == 0) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, newton::dimension);
		ode_control = gsl_odeiv2_control_y_new(absAcc, relAcc);
		ode_evolve = gsl_odeiv2_evolve_alloc(newton::dimension);
		ode_system.function = newton::function;
		ode_system.jacobian = newton::jacobian;
		ode_system.dimension = newton::dimension;
		ode_system.params = &params;
		y[0] = constant::AU;
		y[1] = 0;
		y[2] = 0;
		y[3] = -sqrt(mass * constant::G * constant::M_sun / constant::AU) * cos(constant::pi / 200);
		y[4] = sqrt(mass * constant::G * constant::M_sun / constant::AU) * sin(constant::pi / 200);
		y[5] = 0;
	}
	else if (NPSK == 1) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, postnewtonian::dimension);
		ode_control = gsl_odeiv2_control_y_new(absAcc, relAcc);
		ode_evolve = gsl_odeiv2_evolve_alloc(postnewtonian::dimension);
		ode_system.function = postnewtonian::function;
		ode_system.jacobian = postnewtonian::jacobian;
		ode_system.dimension = postnewtonian::dimension;
		ode_system.params = &params;
		y[0] = constant::AU;
		y[1] = 0;
		y[2] = 0;
		y[3] = -sqrt(mass * constant::G * constant::M_sun / constant::AU) * cos(constant::pi / 200);
		y[4] = sqrt(mass * constant::G * constant::M_sun / constant::AU) * sin(constant::pi / 200);
		y[5] = 0;
	}
	else if (NPSK == 2) {
		ode_step = gsl_odeiv2_step_alloc(ode_type, schwarzschild::dimension);
		ode_control = gsl_odeiv2_control_y_new(absAcc, relAcc);
		ode_evolve = gsl_odeiv2_evolve_alloc(schwarzschild::dimension);
		ode_system.function = schwarzschild::function;
		ode_system.jacobian = schwarzschild::jacobian;
		ode_system.dimension = schwarzschild::dimension;
		ode_system.params = &params;
		if (PL) {
			x[0] = 0;
			x[1] = constant::AU;
			x[2] = 0;
			x[3] = 0;
			x[5] = -sqrt(mass * constant::G * constant::M_sun / constant::AU) * cos(constant::pi / 200);
			x[6] = sqrt(mass * constant::G * constant::M_sun / constant::AU) * sin(constant::pi / 200);
			x[7] = 0;
			schwarzschild::c2s(x, y);
			schwarzschild::particle::normalization(y, &params);
		}
		else {
			x[0] = 0;
			x[1] = mass * 442990;
			x[2] = 0;
			x[3] = 0;
			x[5] = 0;
			x[6] = 1;
			x[7] = 0;
			schwarzschild::c2s(x, y);
			schwarzschild::light::normalization(y, &params);
		}
	}
	else {
	}
	int status = 0, count = 0;
	vector<vector<double>> rec;
	vector<double> temp(ode_system.dimension + 2);
	while (status == 0 && t != tFinal) {
		status = gsl_odeiv2_evolve_apply(ode_evolve, ode_control, ode_step, &ode_system, &t, tFinal, &h, y);
		if (progressBar)
			bar.set_progress(100 * t / tFinal);
		if (count++ % recRatio == 0) {
			if (NPSK == 0) {
				for (int i = 0; i < 6; ++i)
					temp[i] = y[i];
				temp[6] = newton::energy(y, &params);
				temp[7] = newton::angularMomentum(y, &params);
			}
			else if (NPSK == 1) {
				for (int i = 0; i < 6; ++i)
					temp[i] = y[i];
				temp[6] = postnewtonian::energy(y, &params);
				temp[7] = postnewtonian::angularMomentum(y, &params);
			}
			else if (NPSK == 2) {
				schwarzschild::s2c(y, &temp[0]);
				temp[8] = schwarzschild::energy(y, &params);
				temp[9] = schwarzschild::angularMomentum(y, &params);
			}
			rec.push_back(temp);
		}
	}
	if (progressBar)
		indicators::show_console_cursor(true);
	IO::NumPySave(rec, output);
	return 0;
}
