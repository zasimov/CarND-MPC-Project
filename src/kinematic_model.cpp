#include <cmath>
#include <cppad/cppad.hpp>
#include "kinematic_model.h"

void calc_state_t_plus_1(Eigen::VectorXd &state_t,
			 const double delta,
			 const double a,
			 const double dt) {
  assert(state_t.size() == kStateDim);

  const double x = state_t[kStateX];
  const double y = state_t[kStateY];
  const double psi = state_t[kStatePsi];
  const double v = state_t[kStateV];
  const double cte = state_t[kStateCte];
  const double epsi = state_t[kStateEpsi];

  State<double> t(x, y, psi, v, cte, epsi);
  State<double> t_plus_1;

  calc_state_t_plus_1(t, delta, a, dt, t_plus_1);

  state_t[kStateX] = t_plus_1.x_;
  state_t[kStateY] = t_plus_1.y_;
  state_t[kStatePsi] = t_plus_1.psi_;
  state_t[kStateV] = t_plus_1.v_;
  state_t[kStateCte] = t_plus_1.cte_;
  state_t[kStateEpsi] = t_plus_1.epsi_;
}
