#ifndef __CONFIG_H
#define __CONFIG_H

#include <cppad/cppad.hpp>
#include <string>

const size_t kN = 10;  // the number of timesteps in the horizon
const CppAD::AD<double> kDt = 0.1; // how much time elapses between actuations

const double kMinDelta = -0.7;
const double kMaxDelta = 0.7;
const double kMinThrottle = -0.6;
const double kMaxThrottle = 1.0;

const double kMaxSpeed = 70;

const std::string kMaxCpuTime = "20.5";

const double kWcte = 1000.0;
const double kWepsi = 1000.0;
const double kWv = 1.0;
const double kWdelta = 10.0;
const double kWa = 10.0;
const double kWddelta = 100.0;
const double kWda = 15.0;

#endif
