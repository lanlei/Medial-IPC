#ifndef Cuda_SVD_CUH_
#define Cuda_SVD_CUH_

#include "CudaHeader.cuh"

namespace CudaSVD
{
#define _gamma 5.828427124 // FOUR_GAMMA_SQUARED = sqrt(8)+3;
#define _cstar 0.923879532 // cos(pi/8)
#define _sstar 0.3826834323 // sin(p/8)
#define SVD_EPSILON 1e-6

	__device__ __forceinline__
		float accurateSqrt(float x);

	__device__ __forceinline__
		void condSwap(bool c, float *X, float *Y);

	__device__ __forceinline__
		void condNegSwap(bool c, float *X, float *Y);

	__device__ __forceinline__
		void multAB(float* A, float* B, float* C);

	__device__ __forceinline__
		void multAtB(float* A, float* B, float* C);

	__device__ __forceinline__
		void quatToMat(float* mat, const float* qV);

	__device__ __forceinline__
		void approximateGivensQuaternion(float a11, float a12, float a22, float *ch, float *sh);

	__device__ __forceinline__
		void jacobiConjugation(const uint32_t x, const uint32_t y, const uint32_t z,
			float *s11,
			float *s21, float *s22,
			float *s31, float *s32, float *s33,
			float * qV);

	__device__ __forceinline__
		float dist2(float x, float y, float z);

	// finds transformation that diagonalizes a symmetric matrix
	__device__ __forceinline__
		void jacobiEigenanlysis( // symmetric matrix
			float *s11,
			float *s21, float *s22,
			float *s31, float *s32, float *s33,
			// quaternion representation of V
			float * qV);

	__device__ __forceinline__
		void sortSingularValues(// matrix that we want to decompose
			float* A,
			// sort V simultaneously
			float* v);

	__device__ __forceinline__
		void QRGivensQuaternion(float a1, float a2, float *ch, float *sh);

	__device__ __forceinline__
		void QRDecomposition(// matrix that we want to decompose
			float* A, float* Q, float* R);

	__device__ __forceinline__
		void svd(float* A, float* U, float* S, float* V);

	/// polar decomposition can be reconstructed trivially from SVD result
	/// A = UP
	__device__ __forceinline__
		void pd(float* A,
			// output U
			float* U,
			// output P
			float* P);

	__device__ __forceinline__
		float dotVV(float* v1, float* v2);


	__device__ __forceinline__
		float getMatrixDeterminant(float* mat);

};

#endif // !CUDA_SVD_H
