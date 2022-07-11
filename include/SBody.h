#ifndef SBODY_H
#define SBODY_H

#include <string>

void interruptHandler(int signum);
void help(double mass, double spin, double NUT, double tFinal, int tStepNumber, double TCal, int metric, int PN, int ray, double absAcc, double relAcc, std::string storeFormat);

#endif
