#include "CudaMatrixOperator.cuh"

__device__ __forceinline__
float getVectorNorm(float* v, int dim)
{
	return normf(dim, v);
}

__device__ __forceinline__
float getVectorDot(float* v1, float* v2, int dim)
{
	float sum = 0;
	for (uint32_t i = 0; i < dim; i++)
		sum += v1[i] * v2[i];
	return sum;
}

__device__ __forceinline__
void getVectorNormalize(float* v, int dim)
{
	float len = getVectorNorm(v, dim);
	for (uint32_t i = 0; i < dim; i++)
		v[i] /= len;
}

__device__ __forceinline__
void getVector3Cross(float* v1, float* v2, float* result)
{
	result[0] = v1[1] * v2[2] - v2[1] * v1[2];
	result[1] = v1[2] * v2[0] - v2[2] * v1[0];
	result[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

__device__ __forceinline__
void getVectorSub(float* v1, float* v2, float* result, int dim)
{
	for (uint32_t i = 0; i < dim; i++)
		result[i] = v1[i] - v2[i];
}

__device__ __forceinline__
void getMutilMatrix(float* lmat, float* rmat, float* result, int dim, int rdim, int cdim)
{
	for (uint32_t c = 0; c < cdim; c++)
	{
		for (uint32_t r = 0; r < rdim; r++)
		{
			int idx = c * rdim + r;
			result[idx] = 0;
			for (uint32_t b = 0; b < dim; b++)
					result[idx] += lmat[b * rdim + r] * rmat[c * dim + b];
		}
	}
}

//un-safe & hard code
__device__ __forceinline__
float getMatrix3Determinant(float* mat)
{
	float determinant = mat[0] * mat[4] * mat[8] - mat[2] * mat[4] * mat[6] + mat[1] * mat[5] * mat[6] + mat[2] * mat[3] * mat[7] - mat[0] * mat[5] * mat[7] - mat[1] * mat[3] * mat[8];
	return determinant;
}

//un-safe & hard code
__device__ __forceinline__
void getMatrix3Inverse(float* mat, float* inv)
{
	float len = getMatrix3Determinant(mat);

	inv[0] = (mat[4] * mat[8] - mat[5] * mat[7]) / len; inv[3] = (mat[5] * mat[6] - mat[3] * mat[8]) / len; inv[6] = (mat[3] * mat[7] - mat[4] * mat[6]) / len;
	inv[1] = (mat[2] * mat[7] - mat[1] * mat[8]) / len; inv[4] = (mat[0] * mat[8] - mat[2] * mat[6]) / len; inv[7] = (mat[1] * mat[6] - mat[0] * mat[7]) / len;
	inv[2] = (mat[1] * mat[5] - mat[2] * mat[4]) / len; inv[5] = (mat[2] * mat[3] - mat[0] * mat[5]) / len; inv[8] = (mat[0] * mat[4] - mat[1] * mat[3]) / len;
}

//un-safe & hard code
__device__ __forceinline__
float getMatrix4Determinant(float* mat)
{
	float len = mat[1] * mat[11] * mat[14] * mat[4] - mat[1] * mat[10] * mat[15] * mat[4] -
		mat[11] * mat[13] * mat[2] * mat[4] + mat[10] * mat[13] * mat[3] * mat[4] -
		mat[0] * mat[11] * mat[14] * mat[5] + mat[0] * mat[10] * mat[15] * mat[5] +
		mat[11] * mat[12] * mat[2] * mat[5] - mat[10] * mat[12] * mat[3] * mat[5] -
		mat[1] * mat[11] * mat[12] * mat[6] + mat[0] * mat[11] * mat[13] * mat[6] +
		mat[1] * mat[10] * mat[12] * mat[7] - mat[0] * mat[10] * mat[13] * mat[7] -
		mat[15] * mat[2] * mat[5] * mat[8] + mat[14] * mat[3] * mat[5] * mat[8] + mat[1] * mat[15] * mat[6] * mat[8] -
		mat[13] * mat[3] * mat[6] * mat[8] - mat[1] * mat[14] * mat[7] * mat[8] + mat[13] * mat[2] * mat[7] * mat[8] +
		mat[15] * mat[2] * mat[4] * mat[9] - mat[14] * mat[3] * mat[4] * mat[9] - mat[0] * mat[15] * mat[6] * mat[9] +
		mat[12] * mat[3] * mat[6] * mat[9] + mat[0] * mat[14] * mat[7] * mat[9] - mat[12] * mat[2] * mat[7] * mat[9];
	return len;
}

//un-safe & hard code
__device__ __forceinline__
void getMatrix4Inverse(float* mat, float* inv)
{
	float len = getMatrix4Determinant(mat);
	inv[0] = (-mat[11] * mat[14] * mat[5] + mat[10] * mat[15] * mat[5] + mat[11] * mat[13] * mat[6] - mat[10] * mat[13] * mat[7] - mat[15] * mat[6] * mat[9] + mat[14] * mat[7] * mat[9]) / len;
	inv[1] = (mat[1] * mat[11] * mat[14] - mat[1] * mat[10] * mat[15] - mat[11] * mat[13] * mat[2] + mat[10] * mat[13] * mat[3] + mat[15] * mat[2] * mat[9] - mat[14] * mat[3] * mat[9]) / len;
	inv[2] = (-mat[15] * mat[2] * mat[5] + mat[14] * mat[3] * mat[5] + mat[1] * mat[15] * mat[6] - mat[13] * mat[3] * mat[6] - mat[1] * mat[14] * mat[7] + mat[13] * mat[2] * mat[7]) / len;
	inv[3] = (mat[11] * mat[2] * mat[5] - mat[10] * mat[3] * mat[5] - mat[1] * mat[11] * mat[6] + mat[1] * mat[10] * mat[7] + mat[3] * mat[6] * mat[9] - mat[2] * mat[7] * mat[9]) / len;
	inv[4] = (mat[11] * mat[14] * mat[4] - mat[10] * mat[15] * mat[4] - mat[11] * mat[12] * mat[6] + mat[10] * mat[12] * mat[7] + mat[15] * mat[6] * mat[8] - mat[14] * mat[7] * mat[8]) / len;
	inv[5] = (-mat[0] * mat[11] * mat[14] + mat[0] * mat[10] * mat[15] + mat[11] * mat[12] * mat[2] - mat[10] * mat[12] * mat[3] - mat[15] * mat[2] * mat[8] + mat[14] * mat[3] * mat[8]) / len;
	inv[6] = (mat[15] * mat[2] * mat[4] - mat[14] * mat[3] * mat[4] - mat[0] * mat[15] * mat[6] + mat[12] * mat[3] * mat[6] + mat[0] * mat[14] * mat[7] - mat[12] * mat[2] * mat[7]) / len;
	inv[7] = (-mat[11] * mat[2] * mat[4] + mat[10] * mat[3] * mat[4] + mat[0] * mat[11] * mat[6] - mat[0] * mat[10] * mat[7] - mat[3] * mat[6] * mat[8] + mat[2] * mat[7] * mat[8]) / len;
	inv[8] = (-mat[11] * mat[13] * mat[4] + mat[11] * mat[12] * mat[5] - mat[15] * mat[5] * mat[8] + mat[13] * mat[7] * mat[8] + mat[15] * mat[4] * mat[9] - mat[12] * mat[7] * mat[9]) / len;
	inv[9] = (-mat[1] * mat[11] * mat[12] + mat[0] * mat[11] * mat[13] + mat[1] * mat[15] * mat[8] - mat[13] * mat[3] * mat[8] - mat[0] * mat[15] * mat[9] + mat[12] * mat[3] * mat[9]) / len;
	inv[10] = (-mat[1] * mat[15] * mat[4] + mat[13] * mat[3] * mat[4] + mat[0] * mat[15] * mat[5] - mat[12] * mat[3] * mat[5] + mat[1] * mat[12] * mat[7] - mat[0] * mat[13] * mat[7]) / len;
	inv[11] = (mat[1] * mat[11] * mat[4] - mat[0] * mat[11] * mat[5] + mat[3] * mat[5] * mat[8] - mat[1] * mat[7] * mat[8] - mat[3] * mat[4] * mat[9] + mat[0] * mat[7] * mat[9]) / len;
	inv[12] = (mat[10] * mat[13] * mat[4] - mat[10] * mat[12] * mat[5] + mat[14] * mat[5] * mat[8] - mat[13] * mat[6] * mat[8] - mat[14] * mat[4] * mat[9] + mat[12] * mat[6] * mat[9]) / len;
	inv[13] = (mat[1] * mat[10] * mat[12] - mat[0] * mat[10] * mat[13] - mat[1] * mat[14] * mat[8] + mat[13] * mat[2] * mat[8] + mat[0] * mat[14] * mat[9] - mat[12] * mat[2] * mat[9]) / len;
	inv[14] = (mat[1] * mat[14] * mat[4] - mat[13] * mat[2] * mat[4] - mat[0] * mat[14] * mat[5] + mat[12] * mat[2] * mat[5] - mat[1] * mat[12] * mat[6] + mat[0] * mat[13] * mat[6]) / len;
	inv[15] = (-mat[1] * mat[10] * mat[4] + mat[0] * mat[10] * mat[5] - mat[2] * mat[5] * mat[8] + mat[1] * mat[6] * mat[8] + mat[2] * mat[4] * mat[9] - mat[0] * mat[6] * mat[9]) / len;
}

__device__ __forceinline__
void getMutilVVT(float* v, float* mat, int dim)
{
	for(uint32_t i = 0; i < 3; i++)
		for (uint32_t j = 0; j < 3; j++)
			mat[i * 3 + j] = v[i] * v[j];
}

__device__ __forceinline__
void getMutilMV(float*mat, float* v, float* result, int rdim, int cdim)
{
	for (uint32_t r = 0; r < rdim; r++)
	{
		result[r] = 0;
		for (uint32_t c = 0; c < cdim; c++)
			result[r] += mat[c * rdim + r] * v[c];
	}
}

__device__ __forceinline__
void getMatrix3EigenValue(float* mat, float* result)
{
	float u[3];
	u[0] = 1.0;
	u[1] = 1.0;
	u[2] = 1.0;

	float v[3];
	getMutilMV(mat, u, v);
	float mk = fmaxf(fmaxf(CUDA_ABS(v[0]), CUDA_ABS(v[1])), CUDA_ABS(v[2]));
	float mk_ = 1.0;
	int k = 0;
	while (CUDA_ABS(mk - mk_) > FLT_MIN && k < 20)
	{
		u[0] = v[0] / mk;
		u[1] = v[1] / mk;
		u[2] = v[2] / mk;
		getMutilMV(mat, u, v);
		mk_ = mk;
		mk = fmaxf(fmaxf(CUDA_ABS(v[0]), CUDA_ABS(v[1])), CUDA_ABS(v[2]));
		u[0] = v[0] / mk;
		u[1] = v[1] / mk;
		u[2] = v[2] / mk;
		k++;
	}
	*result = mk;
}
