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
double flow_exponent(System system, double dt, int steps, int transient = 10000) {
	Vec3<double> p{ 1.0, 1.0, 1.0 };
	for (int n = 0; n < transient; ++n) {
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

// Largest first
struct Spectrum3 {
	double first;
	double second;
	double third;
};

template <typename System>
Spectrum3 flow_spectrum(System system, double dt, int steps) {
	Vec3<double> p{ 1.0, 1.0, 1.0 };
	for (int n = 0; n < 10000; ++n) {
		p = rk4_step(system, p, dt);
	}

	Vec3<double> u{ 1.0, 0.0, 0.0 };
	Vec3<double> v{ 0.0, 1.0, 0.0 };
	Vec3<double> w{ 0.0, 0.0, 1.0 };
	double sum1 = 0.0;
	double sum2 = 0.0;
	double sum3 = 0.0;
	for (int n = 0; n < steps; ++n) {
		u = direction_of(rk4_step(system, with_direction(p, u), dt));
		v = direction_of(rk4_step(system, with_direction(p, v), dt));
		w = direction_of(rk4_step(system, with_direction(p, w), dt));
		p = rk4_step(system, p, dt);

		const double stretch1 = std::sqrt(dot(u, u));
		u = (1.0 / stretch1) * u;

		v = v + (-dot(v, u)) * u;
		const double stretch2 = std::sqrt(dot(v, v));
		v = (1.0 / stretch2) * v;

		w = w + (-dot(w, u)) * u + (-dot(w, v)) * v;
		const double stretch3 = std::sqrt(dot(w, w));
		w = (1.0 / stretch3) * w;

		sum1 += std::log(stretch1);
		sum2 += std::log(stretch2);
		sum3 += std::log(stretch3);
	}
	const double time = steps * dt;
	return Spectrum3{ sum1 / time, sum2 / time, sum3 / time };
}
