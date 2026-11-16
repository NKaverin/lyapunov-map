#pragma once

#include "systems.hpp"

// 4th order Runge-Kutta
template <typename System, typename State>
State rk4_step(System system, State s, double dt) {
	const State k1 = system.rate(s);
	const State k2 = system.rate(s + (dt / 2.0) * k1);
	const State k3 = system.rate(s + (dt / 2.0) * k2);
	const State k4 = system.rate(s + dt * k3);
	return s + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

template <typename System, typename State>
State flow(System system, State s, double t, double dt) {
	const int steps = static_cast<int>(t / dt + 0.5);
	for (int n = 0; n < steps; ++n) {
		s = rk4_step(system, s, dt);
	}
	return s;
}