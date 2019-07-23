#include <stdio.h>

__global__ void cuda_hello(){
    printf("Hello World from GPU!\n");
}

void showCUDAStats( void )
{
	// code from https://devblogs.nvidia.com/how-query-device-properties-and-handle-errors-cuda-cc/
	int nDevices;

	cudaGetDeviceCount(&nDevices);
	for (int i = 0; i < nDevices; i++) {
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		printf("Device Number: %d\n", i);
		printf("  Device name: %s\n", prop.name);
		printf("  Memory Clock Rate (KHz): %d\n",
		       prop.memoryClockRate);
		printf("  Memory Bus Width (bits): %d\n",
		       prop.memoryBusWidth);
		printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
		       2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
	}

	cuda_hello<<<1,1>>>();
	fflush(stdout);
}
