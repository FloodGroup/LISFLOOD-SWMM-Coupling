#include "unifiedallocator.cuh"
#include "cuda_utils.cuh"

template<typename T>
T* UnifiedAllocator<T>::allocate(std::size_t n)
{
	return static_cast<T*>(malloc_unified(sizeof(T) * n));
}

template<typename T>
void UnifiedAllocator<T>::deallocate(T* ptr, std::size_t unused)
{
	free_unified(ptr);
}
