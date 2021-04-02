// Main program

#include "SBody.h"

#include <chrono>
#include <cmath>
#include <iostream>
#include <thread>
#include <vector>

#include <unistd.h>

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
	cout << "Usage:" << endl;
	cout << "  SBody -h" << endl;
	cout << endl;
	cout << "Options:" << endl;
	cout << "  -h --help: the help information" << endl;
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
	double h = 1e-3;
	double t = 0;
	double mass = 1, spin = 0;
	double absAcc = 1e-15, relAcc = 1e-15; //Absolute & relative accuracy
	double tFinal = 1e10;				   //Time limit
	double tCal;						   //Calculation time limit
	int NPSK;							   //Newton/PN/Schwarzschild/Kerr
	int PL;								   //Particle/Light
	int recRatio;						   //N integration steps per record
	int progressBar;					   //Switch of the progress bar
	int storeFormat;					   //The way store the output
	string output = "output";			   //File name of the output
	//getopt
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
	switch (NPSK) {
	case 0:
		ode_step = gsl_odeiv2_step_alloc(ode_type, newton::dimension);
		ode_control = gsl_odeiv2_control_y_new(absAcc, relAcc);
		ode_evolve = gsl_odeiv2_evolve_alloc(newton::dimension);
		ode_system.function = newton::function;
		ode_system.jacobian = newton::jacobian;
		ode_system.dimension = newton::dimension;
		ode_system.params = &params;
		break;
	case 1:
		ode_step = gsl_odeiv2_step_alloc(ode_type, postnewtonian::dimension);
		ode_control = gsl_odeiv2_control_y_new(absAcc, relAcc);
		ode_evolve = gsl_odeiv2_evolve_alloc(postnewtonian::dimension);
		ode_system.function = postnewtonian::function;
		ode_system.jacobian = postnewtonian::jacobian;
		ode_system.dimension = postnewtonian::dimension;
		ode_system.params = &params;
		break;
	default:
	case 2:
		ode_step = gsl_odeiv2_step_alloc(ode_type, schwarzschild::dimension);
		ode_control = gsl_odeiv2_control_y_new(absAcc, relAcc);
		ode_evolve = gsl_odeiv2_evolve_alloc(schwarzschild::dimension);
		ode_system.function = schwarzschild::function;
		ode_system.jacobian = schwarzschild::jacobian;
		ode_system.dimension = schwarzschild::dimension;
		ode_system.params = &params;
		break;
	case 3:
		break;
	}
	switch (PL) {
	case 1:
		x[0] = 0;
		x[1] = mass * 442990;
		x[2] = 0;
		x[3] = 0;
		x[5] = 0;
		x[6] = 1;
		x[7] = 0;
		schwarzschild::c2s(x, y);
		schwarzschild::light::normalization(y, &params);
		break;
	default:
		x[0] = 0;
		x[1] = constant::AU;
		x[2] = 0;
		x[3] = 0;
		x[5] = -sqrt(mass * constant::G * constant::M_sun / constant::AU) * cos(constant::pi / 200);
		x[6] = sqrt(mass * constant::G * constant::M_sun / constant::AU) * sin(constant::pi / 200);
		x[7] = 0;
		schwarzschild::c2s(x, y);
		schwarzschild::particle::normalization(y, &params);
		break;
	}
	int status = 0;
	vector<vector<double>> rec;
	vector<double> temp(ode_system.dimension + 2);
	while (status == 0 && t != tFinal) {
		status = gsl_odeiv2_evolve_apply(ode_evolve, ode_control, ode_step, &ode_system, &t, tFinal, &h, y);
		schwarzschild::s2c(y, &temp[0]);
		if (progressBar)
			bar.set_progress(100 * t / tFinal);
		temp[ode_system.dimension] = schwarzschild::energy(y, &params);
		temp[ode_system.dimension + 1] = schwarzschild::angularMomentum(y, &params);
		rec.push_back(temp);
	}
	if (progressBar)
		indicators::show_console_cursor(true);
	IO::NumPySave(rec, output);
	return 0;
}
