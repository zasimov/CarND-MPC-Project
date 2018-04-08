# The Model

The model includes a state vector, actuator outputs and the equations to
predict a state based on the "current" state and the "current" actuator outputs.

The state vector is

  * `x` - the vehicle's x coordinate
  * `y` - the vehicle's y coordinate
  * `psi` - the vehicle's orientation angle
  * `velocity` - the vehicle's velocity
  * `cte` - the cross-track error
  * `epsi` - the psi error

The actuators are

  * `delta` - the steering angle
  * `a` - the acceleration (throttle)

The update equations:

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

The update equations are implemented in `kinematic_model.h`


# Timestep Length and Elapsed Duration

The chosen value of `N` is `10` and the chosen value of `dt` is `0.1`.
The effect is that the controller uses `1` second horizon to build control
trajectory. I think `1` second is enough. `dt` was chosen the same as
the latency of actuators.

I tried

  * N=5, dt=0.1 - the car gets pulled away from the track
  * N=20, dt=0.1 - there are no improvements over N=10, just much more computations

You can change `N` and `dt` in `config.h`.


# Polynomial fitting

The waypoints were transformed to the vehicle coordinate system (see
`main.cpp`, line 180). It is useful because the position of the
vehicle can be used as the origin (0, 0).


# Latency

Controller uses the "delayed" state to find the control trajectory. See `main.cpp`, line 136.
