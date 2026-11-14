#include <cmath>
#include <iostream>
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// Largest Lyapunov exponent of the map x -> r x (1 - x):
// the average of log|f'(x)| along one orbit, with f'(x) = r ( 1 - 2x)

double logistic_exponent(double r, int steps) {
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

double alternating_exponent(double a, double b, int steps) {
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

struct Rgb {
	unsigned char r;
	unsigned char g;
	unsigned char b;
};

// A rectangle
struct Grid {
	double a_min;
	double a_max;
	double b_min;
	double b_max;
	int width;
	int height;

	double a_at(int px) const {
		return a_min + (a_max - a_min) * px / (width - 1);
	}

	// Row 0 is the top of the image
	double b_at(int py) const {
		return b_max - (b_max - b_min) * py / (height - 1);
	}
};

// Dark red = chaos, dark blue = strongly stable
Rgb exponent_colour(double exponent) {
	if (exponent > 0.0) {
		// ln2	
		double t = exponent / 0.7;
		if (t > 1.0) {
			t = 1.0;
		}
		return Rgb{
			static_cast<unsigned char>(255 - 115 * t),
			static_cast<unsigned char>(255 - 225 * t),
			static_cast<unsigned char>(255 - 225 * t)
		};
	}
	double t = -exponent / 2.0;
	if (t > 1.0) {
		t = 1.0;
	}
	return Rgb{
		static_cast<unsigned char>(255 - 235 * t),
		static_cast<unsigned char>(255 - 200 * t),
		static_cast<unsigned char>(255 - 115 * t)
	};
}

// Colours pixels of the grid, and writes a PNG, returns false otherwise
bool write_map(Grid grid, const char* filename, int steps) {
	std::vector<unsigned char> image(grid.width * grid.height * 3);

	for (int py = 0; py < grid.height; ++py) {
		const double b = grid.b_at(py);
		for (int px = 0; px < grid.width; ++px) {
			const double a = grid.a_at(px);
			const Rgb colour = exponent_colour(alternating_exponent(a, b, steps));
			const int i = (py * grid.width + px) * 3;
			image[i] = colour.r;
			image[i + 1] = colour.g;
			image[i + 2] = colour.b;
		}
	}

	return stbi_write_png(filename, grid.width, grid.height, 3, image.data(), grid.width * 3) != 0;
}

struct Dual {
	double value;
	double derivative;

	// Derivative of a constant 0
	Dual(double v, double d = 0.0) : value(v), derivative(d) {}
};

Dual operator+(Dual a, Dual b) {
	return Dual(a.value + b.value, a.derivative + b.derivative);
}

Dual operator-(Dual a, Dual b) {
	return Dual(a.value - b.value, a.derivative - b.derivative);
}

// (ab)' = a'b + ab'
Dual operator*(Dual a, Dual b) {
	return Dual(a.value * b.value, a.derivative * b.value + a.value * b.derivative);
}

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

// A direction (vx, vy) along with the poing, the log is recorded and the direction is scaled to 1
template <typename System>
double plane_exponent(System system, int steps) {
	Vec2<double> p{ 0.0, 0.0 };

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

int main() {
	const int width = 512;
	const int height = 512;

	std::vector<unsigned char> pixels(width * height * 3);

	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			const int i = (y * width + x) * 3;
			pixels[i] = static_cast<unsigned char>(x * 255 / (width - 1));
			pixels[i + 1] = static_cast<unsigned char>(y * 255 / (height - 1));
			pixels[i + 2] = 128;
		}
	}

	// The last argument is the stride
	if (stbi_write_png("gradient.png", width, height, 3, pixels.data(), width * 3) == 0) {
		std::cerr << "coud not write gradient.png\n";
		return 1;
	}
	std::cout << "wrote gradient.png\n";
	std::cout << "r = 4   " << logistic_exponent(4.0, 1000000) << "\n";
	std::cout << "ln 2    " << std::log(2.0) << "\n";
	// ln(0.5)
	std::cout << "r = 2.5 " << logistic_exponent(2.5, 1000000) << "\n";
	// Alternating
	std::cout << "a = b = 4   " << alternating_exponent(4.0, 4.0, 1000000) << "\n";
	std::cout << "a = b = 2.5 " << alternating_exponent(2.5, 2.5, 1000000) << "\n";

	// Dual
	const Dual x(0.3, 1.0);
	const Dual y = 4.0 * x * (1.0 - x);
	std::cout << "dual       " << y.value << " " << y.derivative << "\n";
	std::cout << "by hand    " << 4.0 * 0.3 * (1.0 - 0.3) << " " << 4.0 * (1.0 - 2.0 * 0.3) << "\n";
	std::cout << "Logistic 4   " << map_exponent(Logistic{ 4.0 }, 1000000) << "\n";
	std::cout << "Alternating  " << map_exponent(Alternating{ 2.5, 3.5 }, 1000000) << "\n";
	std::cout << "by hand      " << alternating_exponent(2.5, 3.5, 1000000) << "\n";
	std::cout << "Quadratic 2  " << map_exponent(Quadratic{ 2.0 }, 1000000) << "\n";

	// Henon
	std::cout << "Henon 1.4 0.3 " << plane_exponent(Henon{ 1.4, 0.3 }, 1000000) << "\n";
	std::cout << "Henon 2 0     " << plane_exponent(Henon{ 2.0, 0.0 }, 1000000) << "\n";

	// Bifurcation diagram
	const int plot_width = 1024;
	const int plot_height = 512;
	const double r_min = 2.5;
	const double r_max = 4.0;

	// 255 = white
	std::vector<unsigned char> plot(plot_width * plot_height * 3, 255);

	for (int px = 0; px < plot_width; ++px) {
		const double r = r_min + (r_max - r_min) * px / (plot_width - 1);
		double x = 0.3;

		for (int n = 0; n < 500; ++n) {
			x = r * x * (1.0 - x);
		}

		for (int n = 0; n < 300; ++n) {
			x = r * x * (1.0 - x);
			const int py = static_cast<int>((1.0 - x) * (plot_height - 1));
			const int i = (py * plot_width + px) * 3;
			plot[i] = 0;
			plot[i + 1] = 0;
			plot[i + 2] = 0;
		}
		// Band alng the bottom, red = positive (chaos), blue = negative
		const double exponent = logistic_exponent(r, 2000);
		for (int py = plot_height - 16; py < plot_height; ++py) {
			const int i = (py * plot_width + px) * 3;
			if (exponent > 0.0) {
				plot[i] = 200;
				plot[i + 1] = 30;
				plot[i + 2] = 30;
			}
			else {
				plot[i] = 30;
				plot[i + 1] = 60;
				plot[i + 2] = 200;
			}
		}
	}

	if (stbi_write_png("bifurcation.png", plot_width, plot_height, 3, plot.data(), plot_width * 3) == 0) {
		std::cerr << "could not write bifurcation.png\n";
		return 1;
	}

	std::cout << "wrote bifurcation.png\n";
	
	// Lyapunov
	const Grid full{ 2.0, 4.0, 2.0, 4.0, 512, 512 };
	if (!write_map(full, "lyapunov.png", 1000)) {
		std::cerr << "could not write lyapunov.png\n";
		return 1;
	}
	std::cout << "wrote lyapunov.png\n";
	// Zooming 
	const Grid zoom{ 3.81, 3.87, 3.81, 3.87, 512, 512 };
	if (!write_map(zoom, "zoom.png", 2000)) {
		std::cerr << "could not write zoom.png\n";
		return 1;
	}
	std::cout << "wrote zoom.png\n";

	return 0;
}