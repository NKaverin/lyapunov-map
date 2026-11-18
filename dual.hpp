#pragma once

#include <cmath>

struct Dual {
	double value;
	double derivative;

	// Derivative of a constant 0
	Dual(double v, double d = 0.0) : value(v), derivative(d) {}
};

inline Dual operator+(Dual a, Dual b) {
	return Dual(a.value + b.value, a.derivative + b.derivative);
}

inline Dual operator-(Dual a, Dual b) {
	return Dual(a.value - b.value, a.derivative - b.derivative);
}

// (ab)' = a'b + ab'
inline Dual operator*(Dual a, Dual b) {
	return Dual(a.value * b.value, a.derivative * b.value + a.value * b.derivative);
}

// (cos a)' = - sin (a) a'
inline Dual cos(Dual a) {
	return Dual(std::cos(a.value), -std::sin(a.value) * a.derivative);
}

inline Dual sin(Dual a) {
	return Dual(std::sin(a.value), std::cos(a.value) * a.derivative);
}
