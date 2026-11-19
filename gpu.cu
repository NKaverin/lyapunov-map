#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "exponent.hpp"
#include "image.hpp"
#include "reference.hpp"
#include "systems.hpp"

void check(cudaError_t result, const char* what) {
	if (result != cudaSuccess) {
		std::printf("%s failed: %s\n", what, cudaGetErrorString(result));
		std::exit(1);
	}
}

__global__ void logistic_kernetl(const double* r_values, double* exponents, int count, int steps) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= count) {
		return;
	}

	const double r = r_values[i];
	double x = 0.3;
	for (int n = 0; n < 1000; ++n) {
		x = r * x * (1.0 - x);
	}

	double sum = 0.0;
	for (int n = 0; n < steps; ++n) {
		sum += log(fabs(r * (1.0 - 2.0 * x)));
		x = r * x * (1.0 - x);
	}
	exponents[i] = sum / steps;
}

__global__ void henon_point_kernel(double a, double b, int steps, double* exponent) {
	*exponent = plane_exponent(Henon{ a, b }, steps);
}

__global__ void henon_map_kernel(Grid grid, int steps, double* exponents) {
	const int px = blockIdx.x * blockDim.x + threadIdx.x;
	const int py = blockIdx.y * blockDim.y + threadIdx.y;
	if (px >= grid.width || py >= grid.height) {
		return;
	}
	exponents[py * grid.width + px] = plane_exponent(Henon{ grid.a_at(px), grid.b_at(py) }, steps);
}

int main() {
	int count = 0;
	check(cudaGetDeviceCount(&count), "cudaGetDeviceCount");
	std::printf("devices %d\n", count);

	cudaDeviceProp props;
	check(cudaGetDeviceProperties(&props, 0), "cudaGetDeviceProperties");
	std::printf("name %s\n", props.name);
	std::printf("compute capability %d.%d\n", props.major, props.minor);
	std::printf("memory %zu MB\n", props.totalGlobalMem / (1024 * 1024));
	std::printf("warp size %d\n", props.warpSize);

	const int samples = 1024;
	const int steps = 100000;
	std::vector<double> r_values(samples);
	for (int i = 0; i < samples; ++i) {
		r_values[i] = 2.5 + 1.5 * i / (samples - 1);
	}

	double* device_r = nullptr;
	double* device_exponents = nullptr;
	check(cudaMalloc(&device_r, samples * sizeof(double)), "cudaMalloc r");
	check(cudaMalloc(&device_exponents, samples * sizeof(double)), "cudaMalloc exponents");
	check(cudaMemcpy(device_r, r_values.data(), samples * sizeof(double), cudaMemcpyHostToDevice), "copy r");

	const int threads = 256;
	const int blocks = (samples + threads - 1) / threads;
	logistic_kernetl <<<blocks, threads >>> (device_r, device_exponents, samples, steps);
	check(cudaGetLastError(), "kernel launch");
	check(cudaDeviceSynchronize(), "kernel run");

	std::vector<double> exponents(samples);
	check(cudaMemcpy(exponents.data(), device_exponents, samples * sizeof(double), cudaMemcpyDeviceToHost), "copy exponents");
	check(cudaFree(device_r), "cudaFree r");
	check(cudaFree(device_exponents), "cudaFree exponents");

	// compare with CPU
	double worst = 0.0;
	double worst_r = 0.0;
	int sign_differs = 0;
	for (int i = 0; i < samples; ++i) {
		const double cpu = logistic_exponent(r_values[i], steps);
		const double difference = std::abs(exponents[i] - cpu);
		if (difference > worst) {
			worst = difference;
			worst_r = r_values[i];
		}
		if ((exponents[i] > 0.0) != (cpu > 0.0)) {
			++sign_differs;
		}
	}
	std::printf("gpu r==4 %f\n", exponents[samples - 1]);
	std::printf("largest difference %g at r=%f\n", worst, worst_r);
	std::printf("sign differs at %d of %d\n", sign_differs, samples);

	// Henon
	double* device_point = nullptr;
	check(cudaMalloc(&device_point, sizeof(double)), "cudaMalloc point");
	henon_point_kernel<<<1, 1 >>> (1.3, 0.3, 100000, device_point);
	check(cudaGetLastError(), "point launch");
	double point = 0.0;
	check(cudaMemcpy(&point, device_point, sizeof(double), cudaMemcpyDeviceToHost), "copy pont");
	check(cudaFree(device_point), "cudaFree point");
	std::printf("henon 1.4 0.3 gpu %f cpu %f\n", point, plane_exponent(Henon{ 1.4, 0.3 }, 100000));
	// The same plane as henon_map
	const Grid plane{ 0.0, 2.0, -0.5, 0.5, 768, 384 };
	const int pixels = plane.width * plane.height;
	double* device_map = nullptr;
	check(cudaMalloc(&device_map, pixels * sizeof(double)), "cudaMallc map");

	const dim3 tile(16, 16);
	const dim3 tiles((plane.width + tile.x - 1) / tile.x, (plane.height + tile.y - 1) / tile.y);
	const auto gpu_start = std::chrono::steady_clock::now();
	henon_map_kernel <<<tiles, tile >>> (plane, 1000, device_map);
	check(cudaGetLastError(), "map launch");
	check(cudaDeviceSynchronize(), "map run");
	const std::chrono::duration<double> gpu_time = std::chrono::steady_clock::now() - gpu_start;
	
	std::vector<double> gpu_map(pixels);
	check(cudaMemcpy(gpu_map.data(), device_map, pixels * sizeof(double), cudaMemcpyDeviceToHost), "copy map");
	check(cudaFree(device_map), "cudaFree map");

	const auto cpu_start = std::chrono::steady_clock::now();
	std::vector<double> cpu_map(pixels);
	for (int py = 0; py < plane.height; ++py) {
		for (int px = 0; px < plane.width; ++px) {
			cpu_map[py * plane.width + px] = plane_exponent(Henon{ plane.a_at(px), plane.b_at(py) }, 1000);
		}
	}
	const std::chrono::duration<double> cpu_time = std::chrono::steady_clock::now() - cpu_start;
	std::printf("map gpu %f s cpu %f s\n", gpu_time.count(), cpu_time.count());

	// Px by px
	int identical = 0;
	int escape_differs = 0;
	int map_sign_differs = 0;
	double map_worst = 0.0;
	for (int i = 0; i < pixels; ++i) {
		const double gpu = gpu_map[i];
		const double cpu = cpu_map[i];
		if (std::isnan(gpu) != std::isnan(cpu)) {
			++escape_differs;
			continue;
		}
		if (std::isnan(gpu)) {
			++identical;
			continue;
		}
		if (gpu == cpu) {
			++identical;
		}
		if ((gpu > 0.0) != (cpu > 0.0)) {
			++map_sign_differs;
		}
		if (std::abs(gpu - cpu) > map_worst) {
			map_worst = std::abs(gpu - cpu);
		}
	}
	std::printf("pixels %d identical %d\n", pixels, identical);
	std::printf("escape differs %d sign differs %d\n", escape_differs, map_sign_differs);
	std::printf("largest difference %g\n", map_worst);
	return 0;
}