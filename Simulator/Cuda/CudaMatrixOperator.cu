#include "CudaMatrixOperator.cuh"

__device__ __forceinline__
qeal getVectorNorm(qeal* v, int dim)
{
	qeal x = 0.0;
	for (uint32_t i = 0; i < dim; i++)
		x += v[i] * v[i];
	return CUDA_SQRT(x);
}

__device__ __forceinline__
qeal getVectorDot(qeal* v1, qeal* v2, int dim)
{
	qeal sum = 0;
	for (uint32_t i = 0; i < dim; i++)
		sum += v1[i] * v2[i];
	return sum;
}

__device__ __forceinline__
void getVectorNormalize(qeal* v, int dim)
{
	qeal len = getVectorNorm(v, dim);
	for (uint32_t i = 0; i < dim; i++)
		v[i] /= len;
}

__device__ __forceinline__
void getVector3Cross(qeal* v1, qeal* v2, qeal* result)
{
	result[0] = v1[1] * v2[2] - v2[1] * v1[2];
	result[1] = v1[2] * v2[0] - v2[2] * v1[0];
	result[2] = v1[0] * v2[1] - v2[0] * v1[1];
}

__device__ __forceinline__
void getVectorSub(qeal* v1, qeal* v2, qeal* result, int dim)
{
	for (uint32_t i = 0; i < dim; i++)
		result[i] = v1[i] - v2[i];
}

__device__ __forceinline__
void getMutilMatrix(qeal* lmat, qeal* rmat, qeal* result, int dim, int rdim, int cdim)
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
qeal getMatrix3Determinant(qeal* mat)
{
	qeal determinant = mat[0] * mat[4] * mat[8] - mat[2] * mat[4] * mat[6] + mat[1] * mat[5] * mat[6] + mat[2] * mat[3] * mat[7] - mat[0] * mat[5] * mat[7] - mat[1] * mat[3] * mat[8];
	return determinant;
}

//un-safe & hard code
__device__ __forceinline__
void getMatrix3Inverse(qeal* mat, qeal* inv)
{
	qeal len = getMatrix3Determinant(mat);

	inv[0] = (mat[4] * mat[8] - mat[5] * mat[7]) / len; inv[3] = (mat[5] * mat[6] - mat[3] * mat[8]) / len; inv[6] = (mat[3] * mat[7] - mat[4] * mat[6]) / len;
	inv[1] = (mat[2] * mat[7] - mat[1] * mat[8]) / len; inv[4] = (mat[0] * mat[8] - mat[2] * mat[6]) / len; inv[7] = (mat[1] * mat[6] - mat[0] * mat[7]) / len;
	inv[2] = (mat[1] * mat[5] - mat[2] * mat[4]) / len; inv[5] = (mat[2] * mat[3] - mat[0] * mat[5]) / len; inv[8] = (mat[0] * mat[4] - mat[1] * mat[3]) / len;
}

//un-safe & hard code
__device__ __forceinline__
qeal getMatrix4Determinant(qeal* mat)
{
	qeal len = mat[1] * mat[11] * mat[14] * mat[4] - mat[1] * mat[10] * mat[15] * mat[4] -
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
void getMatrix4Inverse(qeal* mat, qeal* inv)
{
	qeal len = getMatrix4Determinant(mat);
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
void getMutilVVT(qeal* v, qeal* mat, int dim)
{
	for(uint32_t i = 0; i < 3; i++)
		for (uint32_t j = 0; j < 3; j++)
			mat[i * 3 + j] = v[i] * v[j];
}

__device__ __forceinline__
void getMutilMV(qeal*mat, qeal* v, qeal* result, int rdim, int cdim)
{
	for (uint32_t r = 0; r < rdim; r++)
	{
		result[r] = 0;
		for (uint32_t c = 0; c < cdim; c++)
			result[r] += mat[c * rdim + r] * v[c];
	}
}

__device__ __forceinline__
void getMatrix3EigenValue(qeal* mat, qeal* result)
{
	qeal u[3];
	u[0] = 1.0;
	u[1] = 1.0;
	u[2] = 1.0;

	qeal v[3];
	getMutilMV(mat, u, v);
	qeal mk = fmaxf(fmaxf(CUDA_ABS(v[0]), CUDA_ABS(v[1])), CUDA_ABS(v[2]));
	qeal mk_ = 1.0;
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

__device__ __forceinline__
void addMatrix3(qeal* mat0, qeal* mat1, qeal* result, qeal S0, qeal S1)
{
	for (int i = 0; i < 9; i++)
		result[i] = S0 * mat0[i] + S1 * mat1[i];
}

__device__ __forceinline__
bool getMatrix3EigenvalueAndEigenvector(const int dim, qeal* matrix, qeal* eigenvalues, qeal* eigenvectors, qeal tol, int maxIter)
{
	for (int i = 0; i < dim; i++) {
		eigenvectors[i*dim + i] = 1.0f;
		for (int j = 0; j < dim; j++) {
			if (i != j)
				eigenvectors[i*dim + j] = 0.0f;
		}
	}

	int nCount = 0;		//current iteration
	while (1) {
		//find the largest element on the off-diagonal line of the matrix
		double dbMax = matrix[1];
		int nRow = 0;
		int nCol = 1;
		for (int i = 0; i < dim; i++) {			//row
			for (int j = 0; j < dim; j++) {		//column
				double d = fabs(matrix[i*dim + j]);
				if ((i != j) && (d > dbMax)) {
					dbMax = d;
					nRow = i;
					nCol = j;
				}
			}
		}

		if (dbMax < tol)     //precision check 
			break;
		if (nCount > maxIter)       //iterations check
			break;
		nCount++;

		double dbApp = matrix[nRow*dim + nRow];
		double dbApq = matrix[nRow*dim + nCol];
		double dbAqq = matrix[nCol*dim + nCol];
		//compute rotate angle
		double dbAngle = 0.5*atan2(-2 * dbApq, dbAqq - dbApp);
		double dbSinTheta = sin(dbAngle);
		double dbCosTheta = cos(dbAngle);
		double dbSin2Theta = sin(2 * dbAngle);
		double dbCos2Theta = cos(2 * dbAngle);
		matrix[nRow*dim + nRow] = dbApp * dbCosTheta*dbCosTheta +
			dbAqq * dbSinTheta*dbSinTheta + 2 * dbApq*dbCosTheta*dbSinTheta;
		matrix[nCol*dim + nCol] = dbApp * dbSinTheta*dbSinTheta +
			dbAqq * dbCosTheta*dbCosTheta - 2 * dbApq*dbCosTheta*dbSinTheta;
		matrix[nRow*dim + nCol] = 0.5*(dbAqq - dbApp)*dbSin2Theta + dbApq * dbCos2Theta;
		matrix[nCol*dim + nRow] = matrix[nRow*dim + nCol];

		for (int i = 0; i < dim; i++) {
			if ((i != nCol) && (i != nRow)) {
				int u = i * dim + nRow;	//p  
				int w = i * dim + nCol;	//q
				dbMax = matrix[u];
				matrix[u] = matrix[w] * dbSinTheta + dbMax * dbCosTheta;
				matrix[w] = matrix[w] * dbCosTheta - dbMax * dbSinTheta;
			}
		}

		for (int j = 0; j < dim; j++) {
			if ((j != nCol) && (j != nRow)) {
				int u = nRow * dim + j;	//p
				int w = nCol * dim + j;	//q
				dbMax = matrix[u];
				matrix[u] = matrix[w] * dbSinTheta + dbMax * dbCosTheta;
				matrix[w] = matrix[w] * dbCosTheta - dbMax * dbSinTheta;
			}
		}

		//compute eigenvector
		for (int i = 0; i < dim; i++) {
			int u = i * dim + nRow;		//p   
			int w = i * dim + nCol;		//q
			dbMax = eigenvectors[u];
			eigenvectors[u] = eigenvectors[w] * dbSinTheta + dbMax * dbCosTheta;
			eigenvectors[w] = eigenvectors[w] * dbCosTheta - dbMax * dbSinTheta;
		}
	}
	//qeal temp_values[3];
	eigenvalues[0] = matrix[0];
	eigenvalues[1] = matrix[4];
	eigenvalues[2] = matrix[8];
	qeal pdbTmpVec[9];
	for (int i = 0; i < 9; i++)
		pdbTmpVec[i] = eigenvectors[i];
	eigenvectors[0] = pdbTmpVec[6];
	eigenvectors[1] = pdbTmpVec[3];
	eigenvectors[2] = pdbTmpVec[0];

	eigenvectors[3] = pdbTmpVec[7];
	eigenvectors[4] = pdbTmpVec[4];
	eigenvectors[5] = pdbTmpVec[1];

	eigenvectors[6] = pdbTmpVec[8];
	eigenvectors[7] = pdbTmpVec[5];
	eigenvectors[8] = pdbTmpVec[2];
	return true;

}

