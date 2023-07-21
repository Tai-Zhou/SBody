/**
 * @file SBody.h
 * @author Tai Zhou (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2022-07-25
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef SBODY_H
#define SBODY_H

#include <string>

constexpr double MASS = 4.15e6;
constexpr double SPIN = 0.0;
constexpr double NUT = 0.0;
constexpr double T_FINAL = 3600.;
constexpr int T_STEP_NUMBER = 10000;
constexpr int T_CAL = 3600000;
constexpr int METRIC = 1;
constexpr int PN = 1;
constexpr int RAY = 1;
constexpr double ABS_ACC = 1e-15;
constexpr double REL_ACC = 1e-15;
constexpr char FORMAT[6] = "NumPy";
constexpr double KEPLER[6] = {1, 2, 3, 4, 5, 6};

void InterruptHandler(int signum);
void Help();

#endif
