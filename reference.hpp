#pragma once

#include <cmath>

// Hand derived exponents, just a check

// Largest Lyapunov exponent of the map x -> r x (1 - x):
// the average of log|f'(x)| along one orbit, with f'(x) = r ( 1 - 2x)
inline double logistic_exponent(double r, int steps) {
	double x = 0.3;

	// Let the orbit settle on the attractor, before measuring
	for (int n = 0; n < 1000; ++n) {
		x = r * x * (1.0 - x);
	}

	double sum = 0.0;
	for (int n = 0; n < steps; ++n) {
		sum += std::log(std::abs(r * (1.0 - 2.0 * x)));
		x = r * x * (1.0 - x);
	}
	return sum / steps;
}

inline double alternating_exponent(double a, double b, int steps) {
	double x = 0.3;

	for (int n = 0; n < 200; ++n) {
		const double r = (n % 2 == 0) ? a : b;
		x = r * x * (1.0 - x);
	}
	
	double sum = 0.0;
	for (int n = 0; n < steps; ++n) {
		const double r = (n % 2 == 0) ? a : b;
		sum += std::log(std::abs(r * (1.0 - 2.0 * x)));
		x = r * x * (1.0 - x);
	}
	return sum / steps;
}

