#pragma once
#ifndef CUDA_HEADER_H
#define CUDA_HEADER_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda_runtime_api.h>
#include "math.h" // CUDA math library

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/count.h>

#include "cublas.h"
#include "cublas_v2.h"
#include "cusolverDn.h"
#include "cusparse.h"
#include "cusolverSp.h"
#include "cusolverSp_LOWLEVEL_PREVIEW.h"

#define CUDA_CALL(x) { const cudaError_t a = (x); if(a != cudaSuccess){ printf("\nCUDA Error: %s (err_num = %d) \n", cudaGetErrorString(a), a); cudaDeviceReset(); assert(0);}}

#define MAX_THREADS_NUM 1024
#define MAX_BLOCKS_NUM 1024
#define THREADS_NUM 512

#define MAX_CELLS_NUM 10000
#define MAX_CELLS_CONTAuint32_tS_TRIANGLE_NUM 500

#define CUDA_SQRT(d) sqrtf(d)
#define CUDA_ABS(d) fabsf(d)
#define IS_CUDA_FLOAT_ZERO(d) (fabsf(d) < FLT_MIN)



#endif