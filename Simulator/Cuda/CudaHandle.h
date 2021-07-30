#pragma once
#ifndef CudaHanlde_H__
#define CudaHanlde_H__
#include "DataCore.h"
#include <cuda_runtime.h>
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>
#include "cublas.h"
#include "cublas_v2.h"
#include "cusolverDn.h"
#include "cusparse.h"
#include "cusolverSp.h"
#include "cusolverSp_LOWLEVEL_PREVIEW.h"

extern cusolverDnHandle_t dnHandle;
extern cublasFillMode_t dnUplo;
extern cusolverSpHandle_t cusolverSpH;
extern cusparseHandle_t cusparseH;
extern cudaStream_t stream;
extern cusparseMatDescr_t spdescrA;

extern csrcholInfo_t sp_chol_info;

extern cublasHandle_t blasHandle;
extern qeal cublas_pos_one;// done
extern qeal cublas_neg_one;// done
extern qeal cublas_zero;// done
void setCublasAndCuSparse();
void freeCublasAndCusparse();

#endif