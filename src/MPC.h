#ifndef MPC_H
#define MPC_H

#include <vector>
#include <cppad/cppad.hpp>
#include "Eigen-3.3/Eigen/Core"

using namespace std;
using CppAD::AD;

class MPC {
 public:
  MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  void Solve(const Eigen::VectorXd &state, const Eigen::VectorXd &coeffs);

 public:
  double steering_angle_;
  double throttle_;

  vector<double> xs_;
  vector<double> ys_;

 private:
  typedef CPPAD_TESTVECTOR(double) Dvector;
  Dvector vars_;
  Dvector vars_lowerbound_;
  Dvector vars_upperbound_;
  Dvector constraints_lowerbound_;
  Dvector constraints_upperbound_;
};

#endif /* MPC_H */
