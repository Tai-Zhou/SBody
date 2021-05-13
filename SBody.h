#ifndef _SBODY_H
#define _SBODY_H

#include <string>

void interruptHandler(int signum);
void help();
double mass, spin;
double tFinal;
double tRec;
double tCal;
size_t NSK;
size_t PN;
int PL;
int progressBar;
std::string storeFormat;
std::string output;

#endif
