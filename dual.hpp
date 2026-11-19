#pragma once

#include <cmath>
#include "device.hpp"

struct Dual {
	double value;
	double derivative;

	// Derivative of a constant 0
	HOST_DEVICE Dual(double v, double d = 0.0) : value(v), derivative(d) {}
};

HOST_DEVICE inline Dual operator+(Dual a, Dual b) {
	return Dual(a.value + b.value, a.derivative + b.derivative);
}

HOST_DEVICE inline Dual operator-(Dual a, Dual b) {
	return Dual(a.value - b.value, a.derivative - b.derivative);
}

// (ab)' = a'b + ab'
HOST_DEVICE inline Dual operator*(Dual a, Dual b) {
	return Dual(a.value * b.value, a.derivative * b.value + a.value * b.derivative);
}

// (cos a)' = - sin (a) a'
HOST_DEVICE inline Dual cos(Dual a) {
	return Dual(std::cos(a.value), -std::sin(a.value) * a.derivative);
}

HOST_DEVICE inline Dual sin(Dual a) {
	return Dual(std::sin(a.value), std::cos(a.value) * a.derivative);
}
