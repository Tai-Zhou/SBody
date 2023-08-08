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

void Help() {
	PrintlnBold("\n{}SBODY{} ({})", VERSION);
	fmt::println("  A relativistic ray tracing program.");
	PrintlnBold("\n{}OPTIONS{}");
	fmt::println("  -m --mass    [f]: mass of the source [{}] (M_sun)", MASS);
	fmt::println("  -s --spin    [f]: spin of the source [{}]", SPIN);
	fmt::println("  -l --NUT     [f]: NUT charge of the source [{}]", NUT);
	fmt::println("  -K --kepler [6f]: 6 doubles, separated by \",\", in order of:");
	fmt::println("                    semimajor axis (r_g)");
	fmt::println("                    eccentricity");
	fmt::println("                    inclination (deg)");
	fmt::println("                    position angle of periapsis (deg)");
	fmt::println("                    ascending node (deg)");
	fmt::println("                    true anomaly (deg)");
	fmt::println("  -E --ec      [f]: eccentricity of the star");
	fmt::println("  -I --in      [f]: inclination of the star (deg)");
	fmt::println("  -o --pe      [f]: position angle of periapsis of the star (deg)");
	fmt::println("  -O --an      [f]: ascending node of the star (deg)");
	fmt::println("  -t --time    [f]: time limit [{}] (s)", T_FINAL);
	fmt::println("  -k --rec     [i]: record number [{}] (s)", T_STEP_NUMBER);
	fmt::println("  -c --tcal    [i]: time limit of calculation [{}] (ms)", T_CAL);
	fmt::println("  -n --metric  [i]: Newton (0)/Schwarzschild (1)/Kerr (2)/KerrTaubNUT (3) [{}]", METRIC);
	fmt::println("  -P --PN      [i]: PN1 (1) + PN2 (2) + PN2.5 (4) + PN3 (8) + PN3.5 (16) [{}]", PN);
	fmt::println("  -R --ray     [i]: ray tracing, view (1) + shadow (4) + camera (2) [{}]", RAY);
	fmt::println("  -a --abs     [f]: absolute accuracy [{}]", ABS_ACC);
	fmt::println("  -r --rel     [f]: relative accuracy [{}]", REL_ACC);
	fmt::println("  -i --inc     [f]: inclination of the BH");
	fmt::println("  -e --eps     [f]: epsilon of the BH");
	fmt::println("  -f --format  [s]: storage format [{}]", FORMAT);
	fmt::println("  -L --light      : light instead of particle");
	fmt::println("  -b --bar        : show progress bar");
	fmt::println("  -h --help       : this help information");
	PrintlnBold("\n{}EXAMPLE{}\n  $ sbody"); // TODO: add example
	PrintlnBold("\n{}SUPPORT{}\n  github.com/Tai-Zhou/SBody");
}

int main(int argc, char *argv[]) {
	// Application Entry
	auto TStart = chrono::steady_clock::now();
	double mass = 4.261e6, spin = 0., charge = 0., NUT = 0.;
	double tFinal = 3600.; // 128.43325526 for 4.15e6 M_sun
	size_t tStepNumber = T_STEP_NUMBER;
	int TCal = T_CAL;
	size_t metric = 1;
	size_t Hamiltonian = 0;
	size_t PN = 1;
	size_t ray = 1;
	int rest_mass = 1;
	ProgressBar::display_ = true;
	string store_format = "NumPy";
	double inc = M_PI * 0. / 180., eps = 0.;
	// double a = 8.3, e = 0., inclination = M_PI * 162. / 180., periapsis = M_PI * 198.9 / 180. + M_PI_2, ascending_node = M_PI * 25.1 / 180., true_anomaly = M_PI_2; // phi=[2.73633242 3.92974873 3.32166381 3.2093593 3.67372211 5.18824159 | 3.19861806 2.63708292 3.05259405]
	double a = 4.999960135538237, e = 0.884649, inclination = M_PI * 134.567 / 180., periapsis = M_PI * 66.263 / 180., ascending_node = M_PI * 228.171 / 180., true_anomaly = M_PI;
	double input_parameter[10];
	unique_ptr<View> viewPtr;
	unique_ptr<thread> shadowPtr;
	unique_ptr<Camera> cameraPtr;
	unique_ptr<thread> lensPtr;
	const char *optShort = "m:s:l:K:E:I:o:O:t:k:c:n:P:R:a:r:i:e:f:Lbh";
	const char *separator = ",";
	char *parameter;
	const struct option optLong[] = {
		{"mass", required_argument, NULL, 'm'},
		{"spin", required_argument, NULL, 's'},
		{"NUT", required_argument, NULL, 'l'},
		{"kepler", required_argument, NULL, 'K'},
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
		case 'K':
			parameter = strtok(optarg, separator);
			for (int i = 0; i < 10 && parameter != nullptr; parameter = strtok(nullptr, separator))
				input_parameter[i++] = atof(parameter);
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
			ascending_node = atof(optarg);
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
				Help();
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
			Help();
			exit(opt == 'h' ? 0 : 1);
		}
	Unit::Initialize(mass);
	tFinal *= Unit::s;
	tFinal = 18 * Unit::yr;
	double t = 8.3564958220072 * Unit::yr, tStep = 0., tRec = tFinal / tStepNumber;
	if (metric == 0)
		ray = 0;
	shared_ptr<Metric> main_metric = make_shared<Schwarzschild>();
	// metric::setMetric(metric, PN, mass, spin, charge, NUT);
	string strFormat = fmt::format(" ({:f},{:f},{:f})[{:f},{:f}]", spin, charge, NUT, inc, eps);
	if (ProgressBar::display_) {
		indicators::show_console_cursor(false);
		ProgressBar::bars_[0].set_option(indicators::option::ForegroundColor{indicators::Color(metric)});
		ProgressBar::bars_[0].set_option(indicators::option::PrefixText{string("?") + strFormat});
	}
	if (ray & 5) {
		// viewPtr = make_unique<View>(make_unique<Schwarzschild>(HAMILTONIAN), 8180. * Unit::pc, inc, fmt::format("view ({:.1f},{:.1f},{:.1f})[{:f},{:f}]", spin, charge, NUT, inc, eps));
		// viewPtr = make_unique<View>(make_unique<Schwarzschild>(HAMILTONIAN), 8180. * Unit::pc, inc, fmt::format("view helical({:.2f},{:.6f},{:.6f},{:.6f},{:.6f})[{:f},{:f}]", 10.6, M_PI * 0.75, M_PI * 7. / 9., 0.01, 0.45 / 10.6, inc, eps));
		// viewPtr = make_unique<View>(make_unique<Schwarzschild>(HAMILTONIAN), 8180. * Unit::pc, inc, fmt::format("view helical({:f},{:f},{:f},{:f},{:f})[{:f},{:f}]", a, inclination, periapsis, e, ascending_node, inc, eps));
		viewPtr = make_unique<View>(make_unique<Schwarzschild>(), 8180. * Unit::pc, inc, fmt::format("view S2H"));
		if (ray & 4)
			shadowPtr = make_unique<thread>(&View::Shadow, viewPtr.get());
	}
	if (ray & 10) {
		cameraPtr = make_unique<Camera>(make_unique<Schwarzschild>(), 1000, 5e-2, 1.e3, inc, fmt::format("camera ({:.1f},{:.1f},{:.1f})[{:f},{:f}]", spin, charge, NUT, inc, eps));
		if (ray & 8)
			lensPtr = make_unique<thread>(&Camera::Lens, cameraPtr.get());
	}
	// Integrator integ(metric::function, metric::jacobian, metric != 0);
	Star star_0(main_metric, BASE, GEODESIC, Unit::R_sun, 0);
	star_0.InitializeSchwarzschildKeplerianApocenter(a * Unit::mpc, e, inclination, periapsis, ascending_node, inc, eps);
	// star_0.InitializeGeodesic(a, inclination, periapsis, ascending_node, -0.16707659553531468, 0.3822615764261866, inc, eps);
	// star_0.InitializeSchwarzschildKeplerianApocenterHarmonic(a * Unit::mpc, e, inclination, periapsis, ascending_node, inc, eps);
	// star_0.InitializeHelical(10.6, M_PI * 0.75, M_PI * 7. / 9., 0.01, 0.45 / 10.6);
	// Integrator &&integrator = main_metric->GetIntegrator(metric != 0);
	// NumPy rec(main_metric->Name() + strFormat, {12});
	// NumPy rec(fmt::format("HotSpot a={:.1f} e={:.2f} i={:.6f} o={:.2f}", a, e, inclination, periapsis), {12});
	NumPy rec(fmt::format("mcmc test", a, inclination, periapsis, e, ascending_node), {12});
	vector<double> temp(12);
	int status = 0;
	double h = 1., stepPercent = 100. / tStepNumber;
	h = -1;
	if (status = star_0.IntegratorApply(&t, 0., &h); status != 0)
		PrintlnWarning("star_0.IntegratorApply() = {}", status);
	h = 1;
	if (status = star_0.IntegratorReset(); status != 0)
		PrintlnWarning("star_0.IntegratorReset() = {}", status);
	for (int i = 0; i < tStepNumber; ++i) {
		tStep += tRec;
		status = star_0.IntegratorApply(&t, tStep, &h);
		if (status > 0)
			PrintlnError("main status = {}", status);
		star_0.Position(temp.data());
		if (Hamiltonian)
			main_metric->HamiltonianToBase(temp.data());
		SphericalToCartesian(temp.data());
		temp[8] = t / Unit::s;
		temp[9] = star_0.Energy();
		temp[10] = star_0.AngularMomentum();
		temp[11] = star_0.CarterConstant();
		if (ray & 1)
			viewPtr->TraceStar(star_0, nullptr, i, nullptr);
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
