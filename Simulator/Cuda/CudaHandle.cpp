#include "CudaHandle.h"

cusolverDnHandle_t dnHandle;
cublasFillMode_t dnUplo;
cusolverSpHandle_t cusolverSpH;
cusparseHandle_t cusparseH;
cudaStream_t stream;
cusparseMatDescr_t spdescrA;
csrcholInfo_t sp_chol_info;
cublasHandle_t blasHandle;
qeal cublas_pos_one = 1.0;
qeal cublas_neg_one = -1.0;
qeal cublas_zero = 0.0;

void setCublasAndCuSparse()
{
	cusolverDnCreate(&dnHandle);
	dnUplo = CUBLAS_FILL_MODE_LOWER;
	cublasCreate(&blasHandle);
	cusolverSpCreate(&cusolverSpH);
	cudaStreamCreate(&stream);
	cusolverSpSetStream(cusolverSpH, stream);
	cusparseCreate(&cusparseH);
	cusparseSetStream(cusparseH, stream);
	cusparseCreateMatDescr(&spdescrA);
	cusolverSpCreateCsrcholInfo(&sp_chol_info);
	cusparseSetMatType(spdescrA, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(spdescrA, CUSPARSE_INDEX_BASE_ZERO);
}

void freeCublasAndCusparse()
{
	cusolverDnDestroy(dnHandle);
	cublasDestroy(blasHandle);
	cusolverSpDestroy(cusolverSpH);
	cudaStreamDestroy(stream);
	cusparseDestroy(cusparseH);
	cusolverSpDestroyCsrcholInfo(sp_chol_info);
}