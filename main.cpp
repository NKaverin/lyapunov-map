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
	std::cout << "r = 2.5 " << logistic_exponent(4.0, 1000000) << "\n";

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
	}

	if (stbi_write_png("bifurcation.png", plot_width, plot_height, 3, plot.data(), plot_width * 3) == 0) {
		std::cerr << "could not write bifurcation.png\n";
		return 1;
	}

	std::cout << "wrote bifurcation.png\n";

	return 0;
}