#include <cstdio>
#include <cstdlib>

#include <cuda_runtime.h>

void check(cudaError_t result, const char* what) {
	if (result != cudaSuccess) {
		std::printf("%s failed: %s\n", what, cudaGetErrorString(result));
		std::exit(1);
	}
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
	return 0;

}