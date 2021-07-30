#ifndef CUDA_MATRIX_OPERATOR_CUH
#define CUDA_MATRIX_OPERATOR_CUH

#include "CudaHeader.cuh"

__device__ __forceinline__
float getVectorNorm(float* v, int dim = 3);

__device__ __forceinline__
float getVectorDot(float* v1, float* v2, int dim = 3);

__device__ __forceinline__
void getVectorNormalize(float* v, int dim = 3);

__device__ __forceinline__
void getVector3Cross(float* v1, float* v2, float* result);

__device__ __forceinline__
void getVectorSub(float* v1, float* v2, float* result, int dim = 3);

__device__ __forceinline__
void getMutilMatrix(float* lmat, float* rmat, float* result, int dim = 3, int rdim = 3, int cdim = 3);

//un-safe & hard code
__device__ __forceinline__
float getMatrix3Determinant(float* mat3);

//un-safe & hard code
__device__ __forceinline__
void getMatrix3Inverse(float* mat3, float* inv);

//un-safe & hard code
__device__ __forceinline__
float getMatrix4Determinant(float* mat4);

//un-safe & hard code
__device__ __forceinline__
void getMatrix4Inverse(float* mat4, float* inv);

__device__ __forceinline__
void getMutilVVT(float* v, float* mat, int dim = 3);

__device__ __forceinline__
void getMutilMV(float*mat, float* v, float* result, int rdim = 3, int cdim = 3);

__device__ __forceinline__
void getMatrix3EigenValue(float* mat, float* result);

#endif