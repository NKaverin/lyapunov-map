#pragma once

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

