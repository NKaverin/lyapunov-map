#pragma once

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


