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

template <typename D>
struct State {
  D x_;
  D y_;
  D psi_;
  D v_;
  D cte_;
  D epsi_;

  State() {
    // I don't care about unitialized variables =)
  }

  State(D x, D y, D psi, D v, D cte, D epsi)
    : x_(x), y_(y), psi_(psi), v_(v), cte_(cte), epsi_(epsi) {
  }
};

/*
 * Calculate state(t+1) using state(t);
 */
template <typename D>
inline void calc_state_t_plus_1(const State<D> &t,
			 const D delta,
			 const D a,
			 const D dt,
			 State<D> &t_plus_1) {
  const auto x = t.x_;
  const auto y = t.y_;
  const auto psi = t.psi_;
  const auto v = t.v_;
  const auto cte = t.cte_;
  const auto epsi = t.epsi_;

  const auto vlf = v / Lf;
  const auto vdt = v * dt;
  const auto ddt = (-delta) * dt;

  t_plus_1.x_ = x + CppAD::cos(psi) * vdt;
  t_plus_1.y_ = y + CppAD::sin(psi) * vdt;
  t_plus_1.psi_ = psi + vlf * ddt;
  t_plus_1.v_ = v + a * dt;
  t_plus_1.cte_ = cte + CppAD::sin(epsi) * vdt;
  t_plus_1.epsi_ = epsi + vlf * ddt;
}


void calc_state_t_plus_1(Eigen::VectorXd &state_t,
			 const double delta,
			 const double a,
			 const double dt);

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
