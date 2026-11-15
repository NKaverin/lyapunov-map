#pragma once

#include <cmath>

#include "dual.hpp"
#include "systems.hpp"

template <typename System>
double map_exponent(System system, int steps) {
	double x = 0.3;

	for (int n = 0; n < 1000; ++n) {
		x = system.step(x, n);
	}

	double sum = 0.0;
	for (int n = 0; n < steps; ++n) {
		const Dual y = system.step(Dual(x, 1.0), n);
		sum += std::log(std::abs(y.derivative));
		x = y.value;
	}
	return sum / steps;
}

// A direction (vx, vy) along with the poing, the log is recorded and the direction is scaled to 1
template <typename System>
double plane_exponent(System system, int steps) {
	Vec2<double> p{ 0.1, 0.0 };

	for (int n = 0; n < 1000; ++n) {
		p = system.step(p, n);
	}

	double vx = 1.0;
	double vy = 0.0;
	double sum = 0.0;
	for (int n = 0; n < steps; ++n) {
		const Vec2<Dual> q = system.step(Vec2<Dual>{ Dual(p.x, vx), Dual(p.y, vy) }, n);
		const double dx = q.x.derivative;
		const double dy = q.y.derivative;
		const double stretch = std::sqrt(dx * dx + dy * dy);
		sum += std::log(stretch);
		vx = dx / stretch;
		vy = dy / stretch;
		p = Vec2<double>{ q.x.value, q.y.value };
	}
	return sum / steps;
}


