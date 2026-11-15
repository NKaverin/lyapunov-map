#include <cmath>
#include <iostream>
#include <vector>

#include "dual.hpp"
#include "exponent.hpp"
#include "reference.hpp"
#include "image.hpp"
#include "systems.hpp"

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
	const Grid full{ 2.0, 4.0, 2.0, 4.0, 512, 512 };
	if (!write_map(full, "lyapunov.png", 1000)) {
		std::cerr << "could not write lyapunov.png\n";
		return false;
	}
	
	std::cout << "wrote lyapunov.png\n";
	
	// Zooming 
	const Grid zoom{ 3.81, 3.87, 3.81, 3.87, 512, 512 };
	if (!write_map(zoom, "zoom.png", 2000)) {
		std::cerr << "could not write zoom.png\n";
		return false;
	}

	std::cout << "wrote zoom.png\n";
	return true;
}

// Henon attractor, one orbit
bool write_henon() {
	const int attractor_width = 768;
	const int attractor_height = 512;
	std::vector<unsigned char> attractor(attractor_width* attractor_height * 3, 255);
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
		return 1;
	}

	std::cout << "wrote henon.png\n";

}

int main() {
	if (!write_gradient()) {
		return 1;
	}
	run_checks();
	if (!write_bifurcation() || !write_maps() || !write_henon()) {
		return 1;
	}
	return 0;
}