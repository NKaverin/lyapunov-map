#pragma once

#include <cmath>

// Each system has step(state, n), T is double or Dual
 
struct Logistic {
	double r;

	template <typename T>
	T step(T x, int) const {
		return r * x * (1.0 - x);
	}
};

struct Alternating {
	double a;
	double b;

	template <typename T>
	T step(T x, int n) const {
		const double r = (n % 2 == 0) ? a : b;
		return r * x * (1.0 - x);
	}
};

struct Quadratic {
	double a;

	template <typename T>
	T step(T x, int) const {
		return 1.0 - a * x * x;
	}
};

// A double or dual
template <typename T> 
struct Vec2 {
	T x;
	T y;
};

// The Henon map (x, y) -> (1 - a x^2 + y, b x)
struct Henon {
	double a;
	double b;

	template <typename T>
	Vec2<T> step(Vec2<T> p, int) const {
		return Vec2<T>{ 1.0 - a * p.x * p.x + p.y, b* p.x };
	}
};

template <typename T>
struct Vec3 {
	T x;
	T y;
	T z;
};

template <typename T>
Vec3<T> operator+(Vec3<T> a, Vec3<T> b) {
	return Vec3<T>{ a.x + b.x, a.y + b.y, a.z + b.z };
}

template <typename T>
Vec3<T> operator*(double c, Vec3<T> v) {
	return Vec3<T>{ c * v.x, c * v.y, c * v.z };
}

// The Lorenz equations
struct Lorenz {
	double sigma;
	double rho;
	double beta;

	template <typename T>
	Vec3<T> rate(Vec3<T> s) const {
		return Vec3<T>{
			sigma* (s.y - s.x),
			s.x * (rho - s.z) - s.y,
			s.x * s.y - beta * s.z
		};
	}
};

// Forced Duffing oscillator x'' _ delta x' + alpha x + beta x^3 = gamma cos(omega t)
struct Duffing {
	double delta;
	double alpha;
	double beta;
	double gamma;
	double omega;

	template <typename T>
	Vec3<T> rate(Vec3<T> s) const {
		using std::cos;
		return Vec3<T>{
			s.y,
			gamma* cos(s.z) - delta * s.y - alpha * s.x - beta * s.x * s.x * s.x,
			T(omega)
		};
	}
};
