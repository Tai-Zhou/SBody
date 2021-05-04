// Main program

#ifndef _SBODY_H
#define _SBODY_H

#include <string>

void interruptHandler(int signum);
void help();
double h;
double mass, spin;
double absAcc, relAcc;
double tFinal;
double tRec;
double tCal;
int NPSK;
int PN;
int PL;
int progressBar;
std::string storeFormat;
std::string output;

#endif
