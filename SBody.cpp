// Main program

#include "SBody.h"

#include <cmath>
#include <iostream>
#include <vector>

#include <unistd.h>

#include <gsl/gsl_odeiv2.h>
#include <indicators/progress_bar.hpp>

#include "Constant.h"
#include "IO.h"
#include "Newton.h"
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
	double absAcc = 1e-15, relAcc = 1e-15;
	double h = 1e-3;
	double t = 0, tfinal = 1e9;
	int NSK = 1;
	int PL = 0;
	const gsl_odeiv2_step_type *ode_type = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_step *ode_step;
	gsl_odeiv2_control *ode_control;
	gsl_odeiv2_evolve *ode_evolve;
	gsl_odeiv2_system ode_system;
	double x[8], y[8];
	switch (NSK) {
	case 0:
		ode_step = gsl_odeiv2_step_alloc(ode_type, newton::dimension);
		ode_control = gsl_odeiv2_control_y_new(absAcc, relAcc);
		ode_evolve = gsl_odeiv2_evolve_alloc(newton::dimension);
		ode_system.function = newton::function;
		ode_system.jacobian = newton::jacobian;
		ode_system.dimension = newton::dimension;
		break;
	case 1:
		ode_step = gsl_odeiv2_step_alloc(ode_type, schwarzschild::dimension);
		ode_control = gsl_odeiv2_control_y_new(absAcc, relAcc);
		ode_evolve = gsl_odeiv2_evolve_alloc(schwarzschild::dimension);
		ode_system.function = schwarzschild::function;
		ode_system.jacobian = schwarzschild::jacobian;
		ode_system.dimension = schwarzschild::dimension;
		break;
	default:
		break;
	}
	switch (PL) {
	case 1:
		x[0] = 0;
		x[1] = 442990;
		x[2] = 0;
		x[3] = 0;
		x[5] = 0;
		x[6] = 1;
		x[7] = 0;
		schwarzschild::c2s(x, y);
		schwarzschild::light::normalization(y);
		break;
	default:
		x[0] = 0;
		x[1] = constant::AU;
		x[2] = 0;
		x[3] = 0;
		x[5] = -sqrt(constant::G * constant::M_sun / constant::AU) * cos(constant::pi / 200);
		x[6] = sqrt(constant::G * constant::M_sun / constant::AU) * sin(constant::pi / 200);
		x[7] = 0;
		schwarzschild::c2s(x, y);
		schwarzschild::particle::normalization(y);
		break;
	}
	int status = 0;
	vector<vector<double>> rec;
	vector<double> temp(ode_system.dimension + 2);
	while (status == 0 && t != tfinal) {
		status = gsl_odeiv2_evolve_apply(ode_evolve, ode_control, ode_step, &ode_system, &t, tfinal, &h, y);
		schwarzschild::s2c(y, &temp[0]);
		temp[ode_system.dimension] = schwarzschild::energy(y);
		temp[ode_system.dimension + 1] = schwarzschild::angularMomentum(y);
		rec.push_back(temp);
	}
	IO::NumPySave(rec, "output");
	return 0;
}
