#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  void Solve(const Eigen::VectorXd &state, const Eigen::VectorXd &coeffs);

 public:
  double steering_angle_;
  double throttle_;

  vector<double> xs_;
  vector<double> ys_;
};

#endif /* MPC_H */
