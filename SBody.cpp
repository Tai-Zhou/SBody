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
	signal(SIGINT, interruptHandler);
	double h = 1e-3;
	mass = 4e6;
	spin = 0.5;
	l = 0.5;
	tFinal = 6.5e-4 * 3.15576e7;
	tRec = 1e-4 * tFinal;
	tCal = 3600;
	NSK = 3;
	Hamiltonian = 0;
	PN = 1;
	ray = 2;
	int PL = 0;
	int displayProgressBar = 1;
	storeFormat = "NumPy";
	double t = 0, tStep = 0;
	const char *optShort = "m:s:t:o:c:n:P:R:a:r:i:e:f:lbh";
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
		{"inc", required_argument, NULL, 'i'},
		{"eps", required_argument, NULL, 'e'},
		{"format", required_argument, NULL, 'f'},
		{"light", no_argument, NULL, 'l'},
		{"ray", no_argument, NULL, 'R'},
		{"bar", no_argument, NULL, 'b'},
		{"help", no_argument, NULL, 'h'},
		{NULL, 0, NULL, 0}};
	int opt;
	double inc = M_PI_2, eps = 0.;
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
		case 'i':
			inc = atof(optarg);
			break;
		case 'e':
			eps = atof(optarg);
			break;
		case 'f':
			storeFormat = atoi(optarg);
			break;
		case 'l':
			PL = 0;
			break;
		case 'b':
			displayProgressBar = 1;
			break;
		default:
			help();
		}
	tFinal *= Constant::s;
	tRec *= Constant::s;
	Metric::setMetric(NSK, Hamiltonian, mass, spin, l);
	Metric::Newton::PN = PN;
	if (displayProgressBar)
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
			double a = 80. * Constant::AU, e = 0.88, inclination = M_PI * 45. / 180., ascendingNode = M_PI * 0. / 180., periapsis = M_PI * 0. / 180., trueAnomaly = 3.67372211; //phi=[2.73633242 3.92974873 3.32166381 3.2093593 3.67372211 5.18824159 3.25134498]
			double r = a * (1 - e * e) / (1 + e * cos(trueAnomaly));
			double tp1 = -r * cos(periapsis + trueAnomaly), tp2 = -r * sin(periapsis + trueAnomaly) * cos(inclination);
			double xp1 = tp1 * cos(ascendingNode) - tp2 * sin(ascendingNode), xp2 = tp2 * cos(ascendingNode) + tp1 * sin(ascendingNode), xp3 = -r * sin(periapsis + trueAnomaly) * sin(inclination);
			x[1] = (xp1 * cos(eps) + xp2 * sin(eps)) * cos(inc) + xp3 * sin(inc);
			x[2] = xp2 * cos(eps) - xp1 * sin(eps);
			x[3] = xp3 * cos(inc) - (xp1 * cos(eps) + xp2 * sin(eps)) * sin(inc);
			double vtheta = sqrt((1 - e * e) * mass * a) / r, vr = sign(M_PI - mod2Pi(trueAnomaly)) * sqrt(max(0., 2. * mass / r - mass / a - vtheta * vtheta));
			//vtheta = 0.3822615764261866;
			//vr = -0.16707659553531468;
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
		if (PL)
			Metric::particleNormalization(y);
		else
			Metric::lightNormalization(y, 1.);
		if (Hamiltonian)
			Metric::qdq2qp(y);
	}
	integrator integ(Metric::function, Metric::jacobian, NSK != 0);
	int status = 0;
	Object::star star_0(Constant::R_sun, y, 0);
	Object::objectList.push_back(&star_0);
	view vie(8. * Constant::kpc, inc);
	if (ray & 2)
		return vie.shadow(100);
	camera cam(100, 4e-2, mass * 1000., inc);
	if (ray & 4)
		cam.initialize();
	vector<vector<double>> rec;
	vector<double> temp(12);
	int rayNO = 0;
	while (tStep < tFinal) {
		tStep = min(tStep + tRec, tFinal);
		while (status == 0 && t < tStep)
			status = integ.apply(&t, tStep, &h, star_0.pos);
		if (status > 0)
			cerr << "[!] main status = " << status << endl;
		if (NSK == 0)
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
			vie.traceBack(star_0, rayNO++);
		if (ray & 4)
			cam.traceBack();
		rec.push_back(temp);
		if (displayProgressBar)
			IO::progressBar.set_progress(100 * t / tFinal);
		auto tpass = chrono::steady_clock::now() - start; // nano seconds
		if (tpass.count() >= tCal * 1000000000)
			break;
	}
	if (displayProgressBar) {
		IO::progressBar.set_option(indicators::option::PrefixText{string("Complete(") + to_string(inc) + "," + to_string(eps) + ")"});
		IO::progressBar.mark_as_completed();
		indicators::show_console_cursor(true);
	}
	if (ray & 1) {
		if (abs(spin) < epsilon)
			vie.save("view e" + to_string(eps));
		else
			vie.save(string("view i") + to_string(inc) + "e" + to_string(eps));
	}
	if (ray & 4)
		cam.save(string("camera i") + to_string(inc) + "e" + to_string(eps));
	IO::NumPy<double> output(Metric::name + (abs(spin) < epsilon ? " 0i" : " i") + to_string(inc) + 'e' + to_string(eps));
	output.save(rec);
	return 0;
}
