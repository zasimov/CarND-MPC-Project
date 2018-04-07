#ifndef __CONFIG_H
#define __CONFIG_H

#include <cppad/cppad.hpp>

const size_t kN = 25;  // the number of timesteps in the horizon
const CppAD::AD<double> kDt = 0.05; // how much time elapses between actuations

#endif
