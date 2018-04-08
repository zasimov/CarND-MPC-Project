#ifndef __KINEMATIC_MODEL_H
#define __KINEMATIC_MODEL_H

/*
 * This module contains vehicle model equations.
 */

#include "Eigen-3.3/Eigen/Dense"
#include <cppad/cppad.hpp>

/*
 * A dimension of state space: x, y, psi, v, cte, epsi
 */
const int kStateDim = 6;
/*
 * A number of actuators: delta and throttle (acceleration)
 */
const int kActuatorCount = 2;


// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;


enum {
  kStateX = 0,
  kStateY = 1,
  kStatePsi = 2,
  kStateV = 3,
  kStateCte = 4,
  kStateEpsi = 5
};

/*
 * Calculate state(t+1) using state(t);
 */
template <typename D>
inline void calc_state_t_plus_1(
				// input
				const D x0,
				const D y0,
				const D psi0,
				const D v0,
				const D cte0,
				const D epsi0,
				// control
				const D delta,
				const D a,
				const D dt,
				// output
				D &x1,
				D &y1,
				D &psi1,
				D &v1,
				D &cte1,
				D &epsi1) {

  const auto vlf = v0 / Lf;
  const auto vdt = v0 * dt;
  const auto ddt = (-delta) * dt;

  x1 = x0 + CppAD::cos(psi0) * vdt;
  y1 = y0 + CppAD::sin(psi0) * vdt;
  psi1 = psi0 + vlf * ddt;
  v1 = v0 + a * dt;
  cte1 = cte0 + CppAD::sin(epsi0) * vdt;
  epsi1 = epsi0 + vlf * ddt;
}

/*
 * dfpoly3 calculates a derivate of polynom of degree 3 in a point x
 */
template <typename D>
inline D dfpoly3(const Eigen::VectorXd &coeffs, D x) {
  return 3.0 * coeffs[3] * x * x + 2.0 * coeffs[2] * x + coeffs[1];
}

/*
 * Evaluate a poly.
 */
template <typename D>
inline D poly3eval(const Eigen::VectorXd &coeffs, D x) {
  assert(coeffs.size() == 4);
  const D x2 = x * x;
  const D x3 = x2 * x;
  return coeffs[3] * x3 + coeffs[2] * x2 + coeffs[1] * x + coeffs[0];
}

#endif
