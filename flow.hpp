#pragma once

#include <cmath>

#include "dual.hpp"
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

inline double dot(Vec3<double> a, Vec3<double> b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

//direction component as derivative
inline Vec3<Dual> with_direction(Vec3<double> p, Vec3<double> v) {
	return Vec3<Dual>{ Dual(p.x, v.x), Dual(p.y, v.y), Dual(p.z, v.z) };
}

inline Vec3<double> point_of(Vec3<Dual> q) {
	return Vec3<double> { q.x.value, q.y.value, q.z.value };
}

inline Vec3<double> direction_of(Vec3<Dual> q) {
	return Vec3<double>{ q.x.derivative, q.y.derivative, q.z.derivative };
}

// Largest Lyapunov exponent of a flow
template <typename System>
double flow_exponent(System system, double dt, int steps) {
	Vec3<double> p{ 1.0, 1.0, 1.0 };
	for (int n = 0; n < 10000; ++n) {
		p = rk4_step(system, p, dt);
	}

	Vec3<double> v{ 1.0, 0.0, 0.0 };
	double sum = 0.0;
	for (int n = 0; n < steps; ++n) {
		const Vec3<Dual> q = rk4_step(system, with_direction(p, v), dt);
		const Vec3<double> d = direction_of(q);
		const double stretch = std::sqrt(dot(d, d));
		sum += std::log(stretch);
		v = (1.0 / stretch) * d;
		p = point_of(q);
	}
	return sum / (steps * dt);
}
