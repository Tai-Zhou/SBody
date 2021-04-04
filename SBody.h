// Main program

#ifndef _SBODY_H
#define _SBODY_H

#define VERSION "0.0.210331"

#include <string>

void help();
double h;
double mass, spin;
double absAcc, relAcc;
double tFinal;
double tCal;
int NPSK;
int PN;
int PL;
int recRatio;
int progressBar;
std::string storeFormat;
std::string output;

#endif
