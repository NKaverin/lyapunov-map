#include <chrono>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>

#include "dual.hpp"
#include "flow.hpp"
#include "exponent.hpp"
#include "reference.hpp"
#include "image.hpp"
#include "systems.hpp"

auto seconds_start() {
	const auto start = std::chrono::steady_clock::now();
	return start;
}

double seconds_since(std::chrono::steady_clock::time_point start) {
	const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
	return elapsed.count();
}

// Test image
bool write_gradient() {
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
		return false;
	}

	std::cout << "wrote gradient.png\n";
	return true;
}

// Exponent with known values
void run_checks() {
	// Logistic
	std::cout << "r = 4   " << logistic_exponent(4.0, 1000000) << "\n";
	std::cout << "ln 2    " << std::log(2.0) << "\n";
	// ln(0.5)
	std::cout << "r = 2.5 " << logistic_exponent(2.5, 1000000) << "\n";
	
	// Alternating
	std::cout << "a = b = 4   " << alternating_exponent(4.0, 4.0, 1000000) << "\n";
	std::cout << "a = b = 2.5 " << alternating_exponent(2.5, 2.5, 1000000) << "\n";

	// Dual numbers
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
	
	// The output should be NaN
	std::cout << "Henon 2 0.3   " << plane_exponent(Henon{ 2.0, 0.3 }, 1000000) << "\n";
	
	// Exponents
	const Spectrum both = plane_spectrum(Henon{ 1.4, 0.3 }, 1000000);
	std::cout << "spectrum  " << both.first << " " << both.second << "\n";
	std::cout << "sum       " << both.first + both.second << "\n";
	std::cout << "ln 0.3    " << std::log(0.3) << "\n";
	std::cout << "dimension " << 1.0 - both.first / both.second << "\n";
	
	// Lorenz
	const Lorenz lorenz{ 10.0, 28.0, 8.0 / 3.0 };
	const Vec3<double> rate = lorenz.rate(Vec3<double>{ 1.0, 1.0, 1.0 });
	std::cout << "lorenz rate " << rate.x << " " << rate.y << " " << rate.z << "\n";
	const Vec3<double> start{ 1.0, 1.0, 1.0 };
	const double coarse = flow(lorenz, start, 0.5, 0.02).x;
	const double medium = flow(lorenz, start, 0.5, 0.01).x;
	const double fine = flow(lorenz, start, 0.5, 0.005).x;
	std::cout << "x(0.5)      " << fine << "\n";
	std::cout << "error ratio " << (coarse - medium) / (medium - fine) << "\n";
}

// Bifurcation diagram
bool write_bifurcation() {
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
		return false;
	}

	std::cout << "wrote bifurcation.png\n";
	return true;
}

// Exponent of the alternating logistic map
bool write_maps() {
	const auto full_start = seconds_start();
	const Grid full{ 2.0, 4.0, 2.0, 4.0, 512, 512 };
	if (!write_map(full, "lyapunov.png", 1000, alternating_exponent)) {
		std::cerr << "could not write lyapunov.png\n";
		return false;
	}
	
	std::cout << "wrote lyapunov.png " << seconds_since(full_start) << "s\n";
	
	// Zooming 
	const auto zoom_start = seconds_start();
	const Grid zoom{ 3.81, 3.87, 3.81, 3.87, 512, 512 };
	if (!write_map(zoom, "zoom.png", 2000, alternating_exponent)) {
		std::cerr << "could not write zoom.png\n";
		return false;
	}

	std::cout << "wrote zoom.png " << seconds_since(zoom_start) << "s\n";
	return true;
}

// Henon attractor, one orbit
bool write_henon() {
	const int attractor_width = 768;
	const int attractor_height = 512;
	std::vector<unsigned char> attractor(attractor_width * attractor_height * 3, 255);
	const Henon henon{ 1.4, 0.3 };
	Vec2<double> point{ 0.1, 0.0 };

	for (int n = 0; n < 2000000; ++n) {
		point = henon.step(point, n);
		const int px = static_cast<int>((point.x + 1.5) / 3.0 * (attractor_width - 1));
		const int py = static_cast<int>((0.45 - point.y) / 0.9 * (attractor_height - 1));
		const int i = (py * attractor_width + px) * 3;
		attractor[i] = 0;
		attractor[i + 1] = 0;
		attractor[i + 2] = 0;
	}

	if (stbi_write_png("henon.png", attractor_width, attractor_height, 3, attractor.data(), attractor_width * 3) == 0) {
		std::cerr << "could not write henon.png\n";
		return false;
	}

	std::cout << "wrote henon.png\n";

	return true;
}

// Largest exponent of the Henon map at one pixel
double henon_exponent(double a, double b, int steps) {
	return plane_exponent(Henon{ a, b }, steps);
}

// Exponent of the Henon map. Grey ones = escaped
bool write_henon_map() {
	const auto plane_start = seconds_start();
	const Grid plane{ 0.0, 2.0, -0.5, 0.5, 768, 384 };
	if (!write_map(plane, "henon_map.png", 1000, henon_exponent)) {
		std::cerr << "could not write henon_map.png\n";
		return false;
	}

	std::cout << "wrote_henon_map.png " << seconds_since(plane_start) << "s\n";
	
	// One shrimp, enlarged
	const auto shrimp_start = seconds_start();
	const Grid shrimp{ 1.4, 1.6, 0.1, 0.3, 384, 384 };
	if (!write_map(shrimp, "henon_zoom.png", 1000, henon_exponent)) {
		std::cerr << "could not write henon_zoom.png\n";
		return false;
	}
	std::cout << "wrote henon_zoom.png " << seconds_since(shrimp_start) << "s\n";
	return true;
}

bool write_lorenz() {
	const int width = 768;
	const int height = 512;
	std::vector<unsigned char> image(width * height * 3, 255);
	const Lorenz lorenz{ 10.0, 28.0, 8.0 / 3.0 };
	Vec3<double> s{ 1.0, 1.0, 1.0 };

	for (int n = 0; n < 200000; ++n) {
		s = rk4_step(lorenz, s, 0.005);
		const int px = static_cast<int>((s.x + 25.0) / 50.0 * (width - 1));
		const int py = static_cast<int>((55.0 - s.z) / 55.0 * (height - 1));
		const int i = (py * width + px) * 3;
		image[i] = 0;
		image[i + 1] = 0;
		image[i + 2] = 0;
	}

	if (stbi_write_png("lorenz.png", width, height, 3, image.data(), width * 3) == 0) {
		std::cerr << "could not write lorenz.png\n";
		return false;
	}

	std::cout << "wrote lorenz.png\n";
	return true;
}

int main() {
	std::cout << "threads " << omp_get_max_threads() << "\n";
	if (!write_gradient()) {
		return 1;
	}
	run_checks();
	if (!write_bifurcation() || 
		!write_maps() || 
		!write_henon() || 
		!write_henon_map() ||
		!write_lorenz()) {
		return 1;
	}
	return 0;
}