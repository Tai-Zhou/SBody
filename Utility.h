#ifndef _UTILITY_H
#define _UTILITY_H

double sqr(const double x);
double cub(const double x);
double quad(const double x);
double dot(const double x[], const double y[] = nullptr, int dimension = 3);
double norm(const double x[], int dimension = 3);
void cross(const double x[], const double y[], double z[]);
int c2s(const double x[], const double v[], double r[], double w[]);
int s2c(const double r[], const double w[], double x[], double v[]);
double _0x(const double x);
double _0x1(const double x);

#endif
