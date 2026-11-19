#pragma once

#include <cmath>
#include <limits>

#include "device.hpp"
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

// An orbit whose x leaves this range is treated as gone to infinity
constexpr double escape_radius = 1.0e6;
constexpr double escaped = std::numeric_limits<double>::quiet_NaN();

// A direction (vx, vy) along with the poing, the log is recorded and the direction is scaled to 1
template <typename System>
HOST_DEVICE double plane_exponent(System system, int steps) {
	Vec2<double> p{ 0.1, 0.0 };

	for (int n = 0; n < 1000; ++n) {
		p = system.step(p, n);
		if (!(std::abs(p.x) < escape_radius)) {
			return escaped;
		}
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
		if (!(std::abs(p.x) < escape_radius)) {
			return escaped;
		}
	}
	return sum / steps;
}

// Largest first
struct Spectrum {
	double first;
	double second;
};

template <typename System>
Spectrum plane_spectrum(System system, int steps) {
	const double nan = std::numeric_limits<double>::quiet_NaN();
	Vec2<double> p{ 0.1, 0.0 };

	for (int n = 0; n < 1000; ++n) {
		p = system.step(p, n);
		if (!(std::abs(p.x) < escape_radius)) {
			return Spectrum{ nan, nan };
		}
	}
	double ux = 1.0;
	double uy = 0.0;
	double vx = 0.0;
	double vy = 1.0;
	double sum1 = 0.0;
	double sum2 = 0.0;
	for (int n = 0; n < steps; ++n) {
		const Vec2<Dual> qu = system.step(Vec2<Dual>{ Dual(p.x, ux), Dual(p.y, uy) }, n);
		const Vec2<Dual> qv = system.step(Vec2<Dual>{ Dual(p.x, vx), Dual(p.y, vy) }, n);

		ux = qu.x.derivative;
		uy = qu.y.derivative;
		const double stretch1 = std::sqrt(ux * ux + uy * uy);
		ux /= stretch1;
		uy /= stretch1;

		// Remove from v its part along u
		vx = qv.x.derivative;
		vy = qv.y.derivative;
		const double along = vx * ux + vy * uy;
		vx -= along * ux;
		vy -= along * uy;
		const double stretch2 = std::sqrt(vx * vx + vy * vy);
		vx /= stretch2;
		vy /= stretch2;

		sum1 += std::log(stretch1);
		sum2 += std::log(stretch2);
		p = Vec2<double>{ qu.x.value, qu.y.value };
		if (!(std::abs(p.x) < escape_radius)) {
			return Spectrum{ nan, nan };
		}
	}
	return Spectrum{ sum1 / steps, sum2 / steps };
}

