#ifndef _UTILITY_H
#define _UTILITY_H

struct source {
	const double mass;
	const double spin;
	source(double _m, double _s = 0) : mass(_m), spin(_s) {}
};

double dot(const double x[], const double y[] = nullptr, int dimension = 3);
double norm(const double x[], int dimension = 3);
void cross(const double x[], const double y[], double z[]);
double _0x(const double x);
double _0x1(const double x);

#endif
