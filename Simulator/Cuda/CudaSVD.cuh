#ifndef CudaSVD_H__
#define CudaSVD_H__

#include "CudaHeader.cuh"

namespace CudaSVD
{
#define _gamma 5.828427124 // FOUR_GAMMA_SQUARED = sqrt(8)+3;
#define _cstar 0.923879532 // cos(pi/8)
#define _sstar 0.3826834323 // sin(p/8)
#define SVD_EPSILON 1e-12

	__device__ __forceinline__
		qeal accurateSqrt(qeal x);

	__device__ __forceinline__
		void condSwap(bool c, qeal *X, qeal *Y);

	__device__ __forceinline__
		void condNegSwap(bool c, qeal *X, qeal *Y);

	__device__ __forceinline__
		void multAB(qeal* A, qeal* B, qeal* C);

	__device__ __forceinline__
		void multAtB(qeal* A, qeal* B, qeal* C);

	__device__ __forceinline__
		void quatToMat(qeal* mat, const qeal* qV);

	__device__ __forceinline__
		void approximateGivensQuaternion(qeal a11, qeal a12, qeal a22, qeal *ch, qeal *sh);

	__device__ __forceinline__
		void jacobiConjugation(const uint32_t x, const uint32_t y, const uint32_t z,
			qeal *s11,
			qeal *s21, qeal *s22,
			qeal *s31, qeal *s32, qeal *s33,
			qeal * qV);

	__device__ __forceinline__
		qeal dist2(qeal x, qeal y, qeal z);

	// finds transformation that diagonalizes a symmetric matrix
	__device__ __forceinline__
		void jacobiEigenanlysis( // symmetric matrix
			qeal *s11,
			qeal *s21, qeal *s22,
			qeal *s31, qeal *s32, qeal *s33,
			// quaternion representation of V
			qeal * qV);

	__device__ __forceinline__
		void sortSingularValues(// matrix that we want to decompose
			qeal* A,
			// sort V simultaneously
			qeal* v);

	__device__ __forceinline__
		void QRGivensQuaternion(qeal a1, qeal a2, qeal *ch, qeal *sh);

	__device__ __forceinline__
		void QRDecomposition(// matrix that we want to decompose
			qeal* A, qeal* Q, qeal* R);

	__device__ __forceinline__
		void svd(qeal* A, qeal* U, qeal* S, qeal* V);

	/// polar decomposition can be reconstructed trivially from SVD result
	/// A = UP
	__device__ __forceinline__
		void pd(qeal* A,
			// output U
			qeal* U,
			// output P
			qeal* P);

	__device__ __forceinline__
		qeal dotVV(qeal* v1, qeal* v2);


	__device__ __forceinline__
		qeal getMatrixDeterminant(qeal* mat);


	//
	///
	////
	///
	//

		__device__ __forceinline__
		void fastSvd3x3(
			qeal* A,			// input A     
			qeal* U,	// output U      
			qeal* S, 
			qeal* V	// output V
		);
	

#define gone					1065353216
#define gsine_pi_over_eight		1053028117
#define gcosine_pi_over_eight   1064076127
#define gone_half				0.5f
#define gsmall_number			1.e-12f
#define gtiny_number			1.e-20f
#define gfour_gamma_squared		5.8284273147583007813f

	union un { float f; unsigned int ui; };
	__device__ __forceinline__
		void svd(
			float a11, float a12, float a13, float a21, float a22, float a23, float a31, float a32, float a33,			// input A     
			float &u11, float &u12, float &u13, float &u21, float &u22, float &u23, float &u31, float &u32, float &u33,	// output U      
			float &s11,
			//float &s12, float &s13, float &s21, 
			float &s22,
			//float &s23, float &s31, float &s32, 
			float &s33,	// output S
			float &v11, float &v12, float &v13, float &v21, float &v22, float &v23, float &v31, float &v32, float &v33	// output V
		);




};

#endif // !CUDA_SVD_H
