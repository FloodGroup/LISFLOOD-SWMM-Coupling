#include "cuda_util.cuh"
#include <helper_cuda.h>
#include <cstdlib>
#include <cstdio>
#include <mutex>
#include <unordered_set>

namespace
{
std::mutex fallback_mutex;
std::unordered_set<void*> pageable_fallbacks;
}

void lis::cuda::sync()
{
	checkCudaErrors(cudaDeviceSynchronize());
}


void lis::cuda::peek()
{
	checkCudaErrors(cudaPeekAtLastError());
}


//cudaError_t lis::cuda::sync()  // for debug
//{
//	return cudaDeviceSynchronize();
//}
//
//cudaError_t lis::cuda::peek()
//{
//	return cudaPeekAtLastError();
//}


void lis::cuda::copy
(
	void* dst,
	void* src,
	size_t count
)
{
	if (count == 0) return;
	checkCudaErrors(cudaMemcpy(dst, src, count, cudaMemcpyDefault));
}

void* lis::cuda::malloc_unified
(
	size_t size
)
{
	if (size == 0) return nullptr;
	void* ptr;
	checkCudaErrors(cudaMallocManaged(&ptr, size));
	return ptr;
}

void* lis::cuda::malloc_pinned
(
	size_t size
)
{
	if (size == 0) return nullptr;
	void* ptr = nullptr;
	cudaError_t status = cudaMallocHost(&ptr, size);
	if (status == cudaSuccess) return ptr;

	// Large domains can exceed the Windows/WDDM pinned-memory allowance even
	// when ordinary committed memory is available.  Pageable host memory is
	// slower for snapshots but is correct for synchronous cudaMemcpy calls.
	cudaGetLastError();
	ptr = std::malloc(size);
	if (ptr == nullptr)
	{
		std::fprintf(stderr,
			"Host allocation failed after cudaMallocHost fallback (%zu bytes).\n",
			size);
		std::exit(EXIT_FAILURE);
	}
	{
		std::lock_guard<std::mutex> lock(fallback_mutex);
		pageable_fallbacks.insert(ptr);
	}
	std::fprintf(stderr,
		"WARNING: cudaMallocHost failed; using pageable host memory (%zu bytes).\n",
		size);
	return ptr;
}

void* lis::cuda::malloc_device
(
	size_t size
)
{
	if (size == 0) return nullptr;
	void* ptr;
	checkCudaErrors(cudaMalloc(&ptr, size));
	return ptr;
}

void lis::cuda::free_unified
(
	void* ptr
)
{
	if (ptr == nullptr) return;
	checkCudaErrors(cudaFree(ptr));
}

void lis::cuda::free_pinned
(
	void* ptr
)
{
	if (ptr == nullptr) return;
	{
		std::lock_guard<std::mutex> lock(fallback_mutex);
		auto found = pageable_fallbacks.find(ptr);
		if (found != pageable_fallbacks.end())
		{
			pageable_fallbacks.erase(found);
			std::free(ptr);
			return;
		}
	}
	checkCudaErrors(cudaFreeHost(ptr));
}

void lis::cuda::free_device
(
	void* ptr
)
{
	if (ptr == nullptr) return;
	checkCudaErrors(cudaFree(ptr));
}

int lis::cuda::get_device()
{
	int device;
	checkCudaErrors(cudaGetDevice(&device));
	return device;
}

void lis::cuda::get_device_properties
(
	cudaDeviceProp& properties,
	int device
)
{
	checkCudaErrors(cudaGetDeviceProperties(&properties, device));
}
