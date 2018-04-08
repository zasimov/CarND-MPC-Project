#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "kinematic_model.h"
#include "config.h"

using CppAD::AD;

/*
 * A number of variables
 */
const int kNumberOfVars = kStateDim * kN + kActuatorCount * (kN - 1);
/*
 * A number of constraints.
 */
const int kNumberOfConstraints = kStateDim * kN;


// helpers to work with vars

enum {
  kRowX = 0,    // a number of row for X
  kRowY = 1,    // a number of row for Y
  kRowPsi = 2,  // ...
  kRowV = 3,
  kRowCte = 4,
  kRowEpsi = 5,
  kRowDelta = 6,
  kRowA = 7
};

#define X_VAR(i) (0 + i)
#define Y_VAR(i) (X_VAR(kN) + i)
#define PSI_VAR(i) (Y_VAR(kN) + i)
#define V_VAR(i) (PSI_VAR(kN) + i)
#define CTE_VAR(i) (V_VAR(kN) + i)
#define EPSI_VAR(i) (CTE_VAR(kN) + i)
#define DELTA_VAR(i) (EPSI_VAR(kN) + i)
#define A_VAR(i) (DELTA_VAR(kN - 1) + i)


class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs_;
  FG_eval(const Eigen::VectorXd &coeffs):
    coeffs_(coeffs) {

    // nothing to do
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  inline AD<double> cost(const ADvector &vars) {
    AD<double> cost = 0;

    for (int i = 0; i < kN; ++i) {
      const auto cte = vars[CTE_VAR(i)];
      const auto epsi = vars[EPSI_VAR(i)];
      const auto v = vars[V_VAR(i)] - kMaxSpeed;
      cost += kWcte * cte * cte + kWepsi * epsi * epsi + kWv * v * v;
    }

    for (int i = 0; i < kN - 1; ++i) {
      const auto delta = vars[DELTA_VAR(i)];
      const auto a = vars[A_VAR(i)];
      cost += kWdelta * delta * delta + kWa * a * a;
    }

    for (int i = 0; i < kN - 2; ++i) {
      const auto ddelta = vars[DELTA_VAR(i + 1)] - vars[DELTA_VAR(i)];
      const auto da = vars[A_VAR(i + 1)] - vars[A_VAR(i)];
      cost += kWddelta * ddelta * ddelta + kWda * da * da;
    }

    return cost;
  }

  // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
  void operator()(ADvector& fg, const ADvector& vars) {
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = cost(vars);

    // Fill current state
    fg[1 + X_VAR(0)] = vars[X_VAR(0)];
    fg[1 + Y_VAR(0)] = vars[Y_VAR(0)];
    fg[1 + PSI_VAR(0)] = vars[PSI_VAR(0)];
    fg[1 + V_VAR(0)] = vars[V_VAR(0)];
    fg[1 + CTE_VAR(0)] = vars[CTE_VAR(0)];
    fg[1 + EPSI_VAR(0)] = vars[EPSI_VAR(0)];

    for (int t = 0; t < kN - 1; t++) {
      const auto x_t = vars[X_VAR(t)];
      const auto y_t = vars[Y_VAR(t)];
      const auto psi_t = vars[PSI_VAR(t)];
      const auto v_t = vars[V_VAR(t)];
      const auto cte_t = vars[CTE_VAR(t)];
      const auto epsi_t = vars[EPSI_VAR(t)];
      const auto delta = vars[DELTA_VAR(t)];
      const auto a = vars[A_VAR(t)];

      const auto x_t1 = vars[X_VAR(t + 1)];
      const auto y_t1 = vars[Y_VAR(t + 1)];
      const auto psi_t1 = vars[PSI_VAR(t + 1)];
      const auto v_t1 = vars[V_VAR(t + 1)];
      const auto cte_t1 = vars[CTE_VAR(t + 1)];
      const auto epsi_t1 = vars[EPSI_VAR(t + 1)];

      // calculate y desired and psi desired
      const auto y_des = poly3eval(coeffs_, x_t);
      const auto psi_des = CppAD::atan(dfpoly3(coeffs_, x_t));

      const auto cte = y_des - y_t;
      const auto epsi = psi_t - psi_des;

      AD<double> e_x_t1, e_y_t1, e_psi_t1, e_v_t1, e_cte_t1, e_epsi_t1;
      calc_state_t_plus_1(x_t, y_t, psi_t, v_t, cte, epsi,
			  delta, a, kDt,
			  e_x_t1, e_y_t1, e_psi_t1, e_v_t1, e_cte_t1, e_epsi_t1);

      fg[1 + X_VAR(t + 1)] = (x_t1 - e_x_t1);
      fg[1 + Y_VAR(t + 1)] = (y_t1 - e_y_t1);
      fg[1 + PSI_VAR(t + 1)] = (psi_t1 - e_psi_t1);
      fg[1 + V_VAR(t + 1)] = (v_t1 - e_v_t1);
      fg[1 + CTE_VAR(t + 1)] = (cte_t1 - e_cte_t1);
      fg[1 + EPSI_VAR(t + 1)] = (epsi_t1 - e_epsi_t1);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC()
  : vars_(kNumberOfVars),
    vars_lowerbound_(kNumberOfVars),
    vars_upperbound_(kNumberOfVars),
    constraints_lowerbound_(kNumberOfConstraints),
    constraints_upperbound_(kNumberOfConstraints) {


  // Set lower and upper limits for variables.
  for (int i = 0; i < DELTA_VAR(0); i++) {
    vars_lowerbound_[i] = -1.0e10;
    vars_upperbound_[i] = 1.0e10;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = DELTA_VAR(0); i < A_VAR(0); i++) {
    vars_lowerbound_[i] = kMinDelta;
    vars_upperbound_[i] = kMaxDelta;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = A_VAR(0); i < vars_.size(); i++) {
    vars_lowerbound_[i] = kMinThrottle;
    vars_upperbound_[i] = kMaxThrottle;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  for (int i = 0; i < constraints_lowerbound_.size(); i++) {
    constraints_lowerbound_[i] = 0;
    constraints_upperbound_[i] = 0;
  }
}

/*
 * state is the "current" state
 * coeffs is the coefficients of the fitting polynomial
 */
void MPC::Solve(const Eigen::VectorXd &state, const Eigen::VectorXd &coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  const double x0 = state[kStateX];
  const double y0 = state[kStateY];
  const double psi0 = state[kStatePsi];
  const double v0 = state[kStateV];
  const double cte0 = state[kStateCte];
  const double epsi0 = state[kStateEpsi];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  for (int i = 0; i < vars_.size(); i++) {
    vars_[i] = 0;
  }

  // Initial state
  vars_[X_VAR(0)] = x0;
  vars_[Y_VAR(0)] = y0;
  vars_[PSI_VAR(0)] = psi0;
  vars_[V_VAR(0)] = v0;
  vars_[CTE_VAR(0)] = cte0;
  vars_[EPSI_VAR(0)] = epsi0;

  constraints_lowerbound_[X_VAR(0)] = x0;
  constraints_lowerbound_[Y_VAR(0)] = y0;
  constraints_lowerbound_[PSI_VAR(0)] = psi0;
  constraints_lowerbound_[V_VAR(0)] = v0;
  constraints_lowerbound_[CTE_VAR(0)] = cte0;
  constraints_lowerbound_[EPSI_VAR(0)] = epsi0;

  constraints_upperbound_[X_VAR(0)] = x0;
  constraints_upperbound_[Y_VAR(0)] = y0;
  constraints_upperbound_[PSI_VAR(0)] = psi0;
  constraints_upperbound_[V_VAR(0)] = v0;
  constraints_upperbound_[CTE_VAR(0)] = cte0;
  constraints_upperbound_[EPSI_VAR(0)] = epsi0;


  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          " + kMaxCpuTime + "\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options,
      vars_,
      vars_lowerbound_, vars_upperbound_,
      constraints_lowerbound_, constraints_upperbound_,
      fg_eval,
      solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  if (! ok) {
    std::cerr << "unable to find solution (terminated with code 42)" << std::endl;
    exit(42);
  }

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  steering_angle_ = solution.x[DELTA_VAR(0)];
  throttle_ = solution.x[A_VAR(0)];

  xs_.clear();
  ys_.clear();

  for (int t = 0; t < kN - 1; t++) {
    xs_.push_back(solution.x[X_VAR(t + 1)]);
    ys_.push_back(solution.x[Y_VAR(t + 1)]);
  }
}
