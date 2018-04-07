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


const double ref_v = 70;


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

    // The part of the cost based on the reference state.
    for (int t = 0; t < kN; t++) {
      cost += CppAD::pow(vars[CTE_VAR(t)], 2);
      cost += CppAD::pow(vars[EPSI_VAR(t)], 2);
      cost += CppAD::pow(vars[V_VAR(t)] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (int t = 0; t < kN - 1; t++) {
      cost += CppAD::pow(vars[DELTA_VAR(t)], 2);
      cost += CppAD::pow(vars[A_VAR(t)], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (int t = 0; t < kN - 2; t++) {
      cost += CppAD::pow(vars[DELTA_VAR(t + 1)] - vars[DELTA_VAR(t)], 2);
      cost += CppAD::pow(vars[A_VAR(t + 1)] - vars[A_VAR(t)], 2);
    }

    return cost;
  }

  // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
  void operator()(ADvector& fg, const ADvector& vars) {
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = cost(vars);

    // Define constraints

    // Fill current state
    fg[1 + X_VAR(0)] = vars[X_VAR(0)];
    fg[1 + Y_VAR(0)] = vars[Y_VAR(0)];
    fg[1 + PSI_VAR(0)] = vars[PSI_VAR(0)];
    fg[1 + V_VAR(0)] = vars[V_VAR(0)];
    fg[1 + CTE_VAR(0)] = vars[CTE_VAR(0)];
    fg[1 + EPSI_VAR(0)] = vars[EPSI_VAR(0)];


    for (int t = 0; t < kN - 1; t++) {
      State<AD<double>> state_t(vars[X_VAR(t)],
				vars[Y_VAR(t)],
				vars[PSI_VAR(t)],
				vars[V_VAR(t)],
				vars[CTE_VAR(t)],
				vars[EPSI_VAR(t)]);
      const AD<double> delta = vars[DELTA_VAR(t)];
      const AD<double> a = vars[A_VAR(t)];

      State<AD<double>> state_t_plus_1(vars[X_VAR(t + 1)],
				       vars[Y_VAR(t + 1)],
				       vars[PSI_VAR(t + 1)],
				       vars[V_VAR(t + 1)],
				       vars[CTE_VAR(t + 1)],
				       vars[EPSI_VAR(t + 1)]);

      State<AD<double>> eval_t_plus_1;
      calc_state_t_plus_1(state_t, delta, a, kDt, eval_t_plus_1);

      fg[1 + X_VAR(t + 1)] = (state_t_plus_1.x_ - eval_t_plus_1.x_);
      fg[1 + Y_VAR(t + 1)] = (state_t_plus_1.y_ - eval_t_plus_1.y_);
      fg[1 + V_VAR(t + 1)] = (state_t_plus_1.v_ - eval_t_plus_1.v_);
      fg[1 + PSI_VAR(t + 1)] = (state_t_plus_1.psi_ - eval_t_plus_1.psi_);
      fg[1 + CTE_VAR(t + 1)] = (state_t_plus_1.cte_ - eval_t_plus_1.cte_);
      fg[1 + EPSI_VAR(t + 1)] = (state_t_plus_1.epsi_ - eval_t_plus_1.epsi_);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

/*
 * state is the "current" state
 * coeffs is the coefficients of the fitting polynomial
 */
void MPC::Solve(const Eigen::VectorXd &state, const Eigen::VectorXd &coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  const double x = state[kStateX];
  const double y = state[kStateY];
  const double psi = state[kStatePsi];
  const double v = state[kStateV];
  const double cte = state[kStateCte];
  const double epsi = state[kStateEpsi];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(kNumberOfVars);
  for (int i = 0; i < vars.size(); i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(kNumberOfVars);
  Dvector vars_upperbound(kNumberOfVars);
  // TODO: Set lower and upper limits for variables.

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).
  // NOTE: Feel free to change this to something else.
  for (int i = DELTA_VAR(0); i < A_VAR(0); i++) {
    vars_lowerbound[i] = -0.75;
    vars_upperbound[i] = 0.75;
  }

  // Acceleration/decceleration upper and lower limits.
  // NOTE: Feel free to change this to something else.
  for (int i = A_VAR(0); i < vars.size(); i++) {
    vars_lowerbound[i] = -0.5;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(kNumberOfConstraints);
  Dvector constraints_upperbound(kNumberOfConstraints);
  assert(constraints_lowerbound.size() == constraints_upperbound.size());
  for (int i = 0; i < constraints_lowerbound.size(); i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[X_VAR(0)] = x;
  constraints_lowerbound[Y_VAR(0)] = y;
  constraints_lowerbound[PSI_VAR(0)] = psi;
  constraints_lowerbound[V_VAR(0)] = v;
  constraints_lowerbound[CTE_VAR(0)] = cte;
  constraints_lowerbound[EPSI_VAR(0)] = epsi;

  constraints_upperbound[X_VAR(0)] = x;
  constraints_upperbound[Y_VAR(0)] = y;
  constraints_upperbound[PSI_VAR(0)] = psi;
  constraints_upperbound[V_VAR(0)] = v;
  constraints_upperbound[CTE_VAR(0)] = cte;
  constraints_upperbound[EPSI_VAR(0)] = epsi;


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
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
  if (! ok) {
    std::cerr << "unable to find solution (terminated with code 42)" << std::endl;
    exit(42);
  }

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  std::cout << "resuuulting" << std::endl;

  steering_angle_ = solution.x[DELTA_VAR(0)];
  throttle_ = solution.x[A_VAR(0)];

  xs_.clear();
  ys_.clear();

  for (int t = 0; t < kN - 1; t++) {
    xs_.push_back(solution.x[X_VAR(t + 1)]);
    ys_.push_back(solution.x[Y_VAR(t + 1)]);
  }

  std::cout << "finished" << std::endl;
}
