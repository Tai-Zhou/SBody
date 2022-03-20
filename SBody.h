#ifndef _SBODY_H
#define _SBODY_H

#include <string>

void interruptHandler(int signum);
void help();
double mass, spin, charge, NUT;
double tFinal;
double tRec;
double TCal;
size_t NSK;
size_t Hamiltonian;
size_t PN;
size_t ray;
std::string storeFormat;

#endif
