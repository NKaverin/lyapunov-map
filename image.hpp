#pragma once

#include <cmath>
#include <vector>

#include "device.hpp"
#include "stb_image_write.h"


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

	HOST_DEVICE double a_at(int px) const {
		return a_min + (a_max - a_min) * px / (width - 1);
	}

	// Row 0 is the top of the image
	HOST_DEVICE double b_at(int py) const {
		return b_max - (b_max - b_min) * py / (height - 1);
	}
};

// Dark red = chaos, dark blue = strongly stable
inline Rgb exponent_colour(double exponent) {
	if (std::isnan(exponent)) {
		return Rgb{ 128, 128, 128 };
	}
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
template <typename Exponent>
inline bool write_map(Grid grid, const char* filename, int steps, Exponent exponent) {
	std::vector<unsigned char> image(grid.width * grid.height * 3);

	// Rows are independent
	#pragma omp parallel for schedule(dynamic)
	for (int py = 0; py < grid.height; ++py) {
		const double b = grid.b_at(py);
		for (int px = 0; px < grid.width; ++px) {
			const double a = grid.a_at(px);
			const Rgb colour = exponent_colour(exponent(a, b, steps));
			const int i = (py * grid.width + px) * 3;
			image[i] = colour.r;
			image[i + 1] = colour.g;
			image[i + 2] = colour.b;
		}
	}

	return stbi_write_png(filename, grid.width, grid.height, 3, image.data(), grid.width * 3) != 0;
}

