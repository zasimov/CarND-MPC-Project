# The Model

The model includes a state vector, actuator outputs and equations to
predict state based on current state and current actuator outputs.

State vector is

  * x - the vehicle's x coordinate
  * y - the vehicle's y coordinate
  * psi - the vehicle's orientation angle
  * velocity - the vehicle's velocity
  * cte - the cross-track error
  * epsi - the psi error

Actuators are

  * delta - the steering angle
  * a - the acceleration (throttle)

Equations:

    x(t+1) = x(t) + v(t) * cos(psi(t)) * dt
    y(t+1) = y(t) + v(t) * sin(ps(t0i) * dt
    psi(t+1) = psi(t) + v(t) / Lf * (-delta) * dt
    v(t+1) = v(t) + a(t) * dt
	cte(t+1) = cte(t) - v(t) * sin(epsi(t)) * dt
    epsi(t+1) = epsi(t) +  v(t) / Lf * (-delta) * dt

    Lf - this is the length from front of vehicle to its center-of-gravity

`cte` and `epsi`

    cte = y_des - y
	y_des = f(x(t))
	epsi = psi - psi_des
	psi_des = atan(f'(x))

	f is a fitted poly

The equations are implemented in `kinematic_model.h`


# Timestep Length and Elapsed Duration

`N` chosen as `10` and `dt` chosen as `0.1`. These means that the
controller uses the horizon in `1` second to build control
trajectory. I think 1 second is enough. `dt` was chosen the same as
the latency of actuators.

`N` defines the depth of horizon (how many states we want to predict)
`dt` defines how much time is between predicted states.

I tried

  * N=5, dt=0.1 - the car gets pulled away from a track
  * N=20, dt=0.1 - there is no improvements over N=10, but much more computations

You can change `N` and `dt` in in `config.h`.


# Polynomial fitting

The waypoints were transformed to vehicle coordinate system (see
`main.cpp`). It is useful because the vehicle can be used as the
origin (0, 0).


# Latency

Controller uses "delayed" state to find control trajectory. See `main.cpp` line 136.
