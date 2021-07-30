#include "GpuFunc.cuh"
#include "Simulator\Cuda\CudaMatrixOperator.cu"
#include "Simulator\Cuda\CudaSVD.cu"

namespace MIPC
{	
	__host__ void updatePredictivePos
	(
		int hostDim,
		int* devDim,
		qeal* devXn,
		qeal* devV,
		qeal* devTimeStep,
		qeal* devXtilde
	)
	{
		dim3 blockSize(THREADS_NUM);
		uint32_t size = hostDim;
		uint32_t num_block = (size + (THREADS_NUM - 1)) / THREADS_NUM;
		dim3 gridSize(num_block);

		updatePredictivePos << <gridSize, blockSize >> >
			(
				devDim,
				devXn,
				devV,
				devTimeStep,
				devXtilde
				);
		cudaDeviceSynchronize();
	}

	__global__ void updatePredictivePos
	(
		int* devDim,
		qeal* devXn,
		qeal* devV,
		qeal* devTimeStep,
		qeal* devXtilde
	)
	{
		__shared__ int dim;
		__shared__ qeal timeStep;
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			dim = *devDim;
			timeStep = *devTimeStep;
		}
		__syncthreads();

		for (; tid < dim; tid += length)
		{
			devXtilde[tid] = devXn[tid] + timeStep * devV[tid] /*+ timeStep * timeStep * (devForce[tid] / devTetPointsMass[tid])*/;
		}
	}

	__host__ void updateLineSearchX
	(
		int hostDim,
		int* devDim,
		qeal* devXn,
		qeal* devDir,
		qeal alpha,
		qeal* devX
	)
	{
		qeal* devAlpha;
		cudaMalloc((void**)&devAlpha, sizeof(qeal));
		cudaMemcpy(devAlpha, &alpha, sizeof(qeal), cudaMemcpyHostToDevice);

		dim3 blockSize(THREADS_NUM);
		uint32_t size = hostDim;
		uint32_t num_block = (size + (THREADS_NUM - 1)) / THREADS_NUM;
		dim3 gridSize(num_block);

		updateLineSearchX << <gridSize, blockSize >> >
			(
				devDim,
				devXn,
				devDir,
				devAlpha,
				devX
				);
		cudaDeviceSynchronize();
		cudaFree(devAlpha);
	}

	__global__ void updateLineSearchX
	(
		int* devDim,
		qeal* devXn,
		qeal* devDir,
		qeal* devAlpha,
		qeal* devX
	)
	{
		__shared__ int dim;
		__shared__ qeal alpha;
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			dim = *devDim;
			alpha = *devAlpha;
		}
		__syncthreads();

		for (; tid < dim; tid += length)
		{
			devX[tid] = devXn[tid] + alpha * devDir[tid];
		}
	}

	__host__ void computeInertialDiffX
	(
		int hostDim,
		int* devDim,
		qeal* devX,
		qeal* devXtilde,
		qeal* devTetPointsMass,
		qeal* devDiffX,
		qeal* devMassMulDiffX
	)
	{
		dim3 blockSize(THREADS_NUM);
		uint32_t size = hostDim;
		uint32_t num_block = (size + (THREADS_NUM - 1)) / THREADS_NUM;
		dim3 gridSize(num_block);

		computeInertialDiffX << <gridSize, blockSize >> >
			(
				devDim,
				devX,
				devXtilde,
				devTetPointsMass,
				devDiffX,
				devMassMulDiffX
				);
		cudaDeviceSynchronize();
	}

	__global__ void computeInertialDiffX
	(
		int* devDim,
		qeal* devX,
		qeal* devXtilde,
		qeal* devTetPointsMass,
		qeal* devDiffX,
		qeal* devMassMulDiffX
	)
	{
		__shared__ int dim;
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			dim = *devDim;
		}
		__syncthreads();

		for (; tid < dim; tid += length)
		{
			qeal v = devX[tid] - devXtilde[tid];
			devDiffX[tid] = v;
			devMassMulDiffX[tid] = devTetPointsMass[tid] * v;
		}
	}

	__host__ void computeElementsEnergy
	(
		int hostTetElementNum,
		int* devTetElementNum,
		qeal* devTetElementDisplacement,
		qeal* devTetElementDm,
		qeal* devTetElementInvDm,
		qeal* devTetElementAttri,
		qeal* devTetElementVol,
		qeal* devTimeStep,
		qeal* devTetElementPotentialEnergy
	)
	{
		dim3 blockSize(THREADS_NUM);
		uint32_t size = hostTetElementNum;
		uint32_t num_block = (size + (THREADS_NUM - 1)) / THREADS_NUM;
		dim3 gridSize(num_block);

		computeElementsEnergy << <gridSize, blockSize >> >
			(
				devTetElementNum,
				devTetElementDisplacement,
				devTetElementDm,
				devTetElementInvDm,
				devTetElementAttri,
				devTetElementVol,
				devTimeStep,
				devTetElementPotentialEnergy
				);
		cudaDeviceSynchronize();
	}

	__global__ void computeElementsEnergy
	(
		int* devTetElementNum,
		qeal* devTetElementDisplacement,
		qeal* devTetElementDm,
		qeal* devTetElementInvDm,
		qeal* devTetElementAttri,
		qeal* devTetElementVol,
		qeal* devTimeStep,
		qeal* devTetElementPotentialEnergy
	)
	{
		__shared__ int elementNum;
		__shared__ qeal timeStep;
		const int length = gridDim.x *  blockDim.x;
		int eleId = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			elementNum = *devTetElementNum;
			timeStep = *devTimeStep;
		}
		__syncthreads();

		for (; eleId < elementNum; eleId += length)
		{
			qeal* eleDisplacement = devTetElementDisplacement + 12 * eleId;
			qeal* Dm = devTetElementDm + 9 * eleId;
			qeal* invDm = devTetElementInvDm + 9 * eleId;
			qeal v = devTetElementVol[eleId];
			qeal Ds[9];
			computeElementDs(eleDisplacement, Dm, Ds);
			qeal F[9];
			getMutilMatrix(Ds, invDm, F);

			qeal Ic = F[0] * F[0] + F[3] * F[3] + F[6] * F[6] +
				F[1] * F[1] + F[4] * F[4] + F[7] * F[7] +
				F[2] * F[2] + F[5] * F[5] + F[8] * F[8];

			qeal J = getMatrix3Determinant(F);

			qeal lameMu = devTetElementAttri[3 * eleId + 1], lameLamda = devTetElementAttri[3 * eleId + 2];

			//qeal miu = (4.0 * lameMu) / 3.0;
			//qeal lamda = lameLamda + (5.0 * lameMu) / 6.0;
			qeal mu = lameMu;
			qeal lamda = lameLamda;

			qeal alpha = 1.0 + (3.0 * mu) / (4.0 * lamda);
			qeal e = 0.5 * mu * (Ic - 3.0) + 0.5 * lamda * (J - alpha) *  (J - alpha) - 0.5 * mu * log(Ic + 1);

			devTetElementPotentialEnergy[eleId] = timeStep * timeStep * v * e;
		}
	}

	__host__ void computeReducedInertia
	(
		int hostDim,
		int* devDim,
		qeal* devReducedX,
		qeal* devReducedXn,
		qeal* devReducedVn,
		qeal* devTimeStep,
		qeal* devInertia
	)
	{
		dim3 blockSize(THREADS_NUM);
		uint32_t size = hostDim;
		uint32_t num_block = (size + (THREADS_NUM - 1)) / THREADS_NUM;
		dim3 gridSize(num_block);

		computeReducedInertia << <gridSize, blockSize >> >
			(
				devDim,
				devReducedX,
				devReducedXn,
				devReducedVn,
				devTimeStep,
				devInertia
				);
		cudaDeviceSynchronize();
	}


	__global__ void computeReducedInertia
	(
		int* devDim,
		qeal* devReducedX,
		qeal* devReducedXn,
		qeal* devReducedVn,
		qeal* devTimeStep,
		qeal* devInertia
	)
	{
		__shared__ int dim;
		__shared__ qeal timeStep;
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			dim = *devDim;
			timeStep = *devTimeStep;
		}
		__syncthreads();

		for (; tid < dim; tid += length)
		{
			devInertia[tid] = devReducedX[tid] - devReducedXn[tid] - timeStep * devReducedVn[tid];
		}
	}

	__host__ void updatedVelocity
	(
		int hostDim,
		int* devDim,
		qeal* devX,
		qeal* devXtilde,
		qeal* devTimeStep,
		qeal* devVelocity
	)
	{
		dim3 blockSize(THREADS_NUM);
		uint32_t size = hostDim;
		uint32_t num_block = (size + (THREADS_NUM - 1)) / THREADS_NUM;
		dim3 gridSize(num_block);

		updatedVelocity << <gridSize, blockSize >> >
			(
				devDim,
				devX,
				devXtilde,
				devTimeStep,
				devVelocity
				);
		cudaDeviceSynchronize();
	}

	__global__ void updatedVelocity
	(
		int* devDim,
		qeal* devX,
		qeal* devXtilde,
		qeal* devTimeStep,
		qeal* devVelocity
	)
	{
		__shared__ int dim;
		__shared__ qeal timeStep;
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			dim = *devDim;
			timeStep = 1.0 / (*devTimeStep);
		}
		__syncthreads();

		for (; tid < dim; tid += length)
		{
			devVelocity[tid] += 1.0 * timeStep * (devX[tid] - devXtilde[tid]);
		}
	}


	__host__ void assembleTetELementX
	(
		int hostTetElementNum,
		int* devTetElementNum,
		int*devTetElementIndices,
		qeal* devTetPointsX,
		qeal* devTetElementX
	)
	{
		const int threadsNum = 12 * 64;
		dim3 blockSize(threadsNum); //
		uint32_t size = hostTetElementNum * 12;
		uint32_t num_block = (size + (threadsNum - 1)) / threadsNum;
		dim3 gridSize(num_block);

		assembleTetELementX << <gridSize, blockSize >> >
			(
				devTetElementNum,
				devTetElementIndices,
				devTetPointsX,
				devTetElementX
				);
		cudaDeviceSynchronize();

	}

	__global__ void assembleTetELementX
	(
		int* devTetElementNum,
		int*devTetElementIndices,
		qeal* devTetPointsX,
		qeal* devTetElementX
	)
	{
		__shared__ int dim;
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			dim = 12 * *devTetElementNum;
		}
		__syncthreads();

		for (; tid < dim; tid += length)
		{
			int offset = tid % 12;
			int eleId = (tid - offset) / 12;
			int vIndicesOffset = offset % 3;
			int vIndices = (offset - vIndicesOffset) / 3;

			int vid = devTetElementIndices[4 * eleId + vIndices];
			devTetElementX[tid] = devTetPointsX[3 * vid + vIndicesOffset];
		}

	}
	
	__host__ void assembleTetPointsForceFromElementForce
	(
		int hostTetPointsNum,
		int* devTetPointsNum,
		int* devTetPointsSharedElementNum,
		int* devTetPointsSharedElementOffset,
		int* devTetPointsSharedElementList,
		qeal* devTetElementForce,
		qeal* devTetPointInternalForce
	)
	{
		const int threadsNum = 3 * 256;
		dim3 blockSize(threadsNum); //
		uint32_t size = hostTetPointsNum * 3;
		uint32_t num_block = (size + (threadsNum - 1)) / threadsNum;
		dim3 gridSize(num_block);

		assembleTetPointsForceFromElementForce << <gridSize, blockSize >> >
			(
				devTetPointsNum,
				devTetPointsSharedElementNum,
				devTetPointsSharedElementOffset,
				devTetPointsSharedElementList,
				devTetElementForce,
				devTetPointInternalForce
				);
		CUDA_CALL(cudaDeviceSynchronize());

	}


	__global__ void assembleTetPointsForceFromElementForce
	(
		int* devTetPointsNum,
		int* devTetPointsSharedElementNum,
		int* devTetPointsSharedElementOffset,
		int* devTetPointsSharedElementList,
		qeal* devTetElementForce,
		qeal* devTetPointInternalForce
	)
	{
		__shared__ int dim;
		__shared__ qeal pointForce[768];
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		int vOffset = tid % 3;
		int vid = (tid - vOffset) / 3;

		if (threadIdx.x == 0)
		{
			dim = 3 * *devTetPointsNum;
		}

		__syncthreads();

		for (; tid < dim; tid += length)
		{
			int vOffset = tid % 3;
			int vid = (tid - vOffset) / 3;
			int num = devTetPointsSharedElementNum[vid];
			int offset = devTetPointsSharedElementOffset[vid];

			pointForce[threadIdx.x] = 0;
			for (int i = 0; i < num; i++)
			{
				int eleId = devTetPointsSharedElementList[offset + 2 * i];
				int id = devTetPointsSharedElementList[offset + 2 * i + 1];

				pointForce[threadIdx.x] += devTetElementForce[12 * eleId + 3 * id + vOffset];
			}
			devTetPointInternalForce[tid] = pointForce[threadIdx.x];
		}
	}

	__host__ void computeTetElementInternalForce
	(
		int hostTetElementNum,
		int* devTetElementNum,
		qeal* devTetElementDisplacement,
		qeal* devTetElementDm,
		qeal* devTetElementInvDm,
		qeal* devTetElementdFPK,
		qeal* devTetElementAttri,
		qeal* devTetElementInternalForce
	)
	{
		dim3 blockSize(THREADS_NUM);
		uint32_t size = hostTetElementNum;
		uint32_t num_block = (size + (THREADS_NUM - 1)) / THREADS_NUM;
		dim3 gridSize(num_block);

		computeTetElementDeformationGradient << <gridSize, blockSize >> >
			(
				devTetElementNum,
				devTetElementDisplacement,
				devTetElementDm,
				devTetElementInvDm,
				devTetElementdFPK
				);
		cudaDeviceSynchronize();

		computeTetElementFPK << <gridSize, blockSize >> >
			(
				devTetElementNum,
				devTetElementDisplacement,
				devTetElementDm,
				devTetElementInvDm,
				devTetElementAttri,
				devTetElementdFPK
				);
		cudaDeviceSynchronize();

		blockSize = dim3(576); // 9 * 64
		size = 9 * hostTetElementNum;
		num_block = (size + (576 - 1)) / 576;
		gridSize = dim3(num_block);

		computeTetElementInternalForce << <gridSize, blockSize >> >
			(
				devTetElementNum,
				devTetElementInvDm,
				devTetElementdFPK,
				devTetElementAttri,
				devTetElementInternalForce
				);
		cudaDeviceSynchronize();
	}

	__global__ void computeTetElementFPK
	(
		int* devTetElementNum,
		qeal* devTetElementDisplacement,
		qeal* devTetElementDm,
		qeal* devTetElementInvDm,
		qeal* devTetElementAttri,
		qeal* devTetElementdFPK
	)
	{
		__shared__ int elementNum;
		const int length = gridDim.x *  blockDim.x;
		int eleId = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			elementNum = *devTetElementNum;
		}
		__syncthreads();

		for (; eleId < elementNum; eleId += length)
		{
			/*qeal* eleDisplacement = devTetElementDisplacement + 12 * eleId;
			qeal* Dm = devTetElementDm + 9 * eleId;
			qeal* invDm = devTetElementInvDm + 9 * eleId;*/
			qeal* F = devTetElementdFPK + 9 * eleId;

			/*qeal Ds[9];
			computeElementDs(eleDisplacement, Dm, Ds);
			qeal F[9];
			getMutilMatrix(Ds, invDm, F);*/

			qeal volume = devTetElementAttri[3 * eleId];
			qeal lameMu = devTetElementAttri[3 * eleId + 1], lameLamda = devTetElementAttri[3 * eleId + 2];

			qeal Ic = F[0] * F[0] + F[3] * F[3] + F[6] * F[6] +
				F[1] * F[1] + F[4] * F[4] + F[7] * F[7] +
				F[2] * F[2] + F[5] * F[5] + F[8] * F[8];

			qeal J = getMatrix3Determinant(F);

			qeal dJdF[9];
			dJdF[0] = F[4] * F[8] - F[7] * F[5];   dJdF[3] = F[7] * F[2] - F[1] * F[8];	  dJdF[6] = F[1] * F[5] - F[4] * F[2];
			dJdF[1] = F[6] * F[5] - F[3] * F[8];   dJdF[4] = F[0] * F[8] - F[6] * F[2];   dJdF[7] = F[3] * F[2] - F[0] * F[5];
			dJdF[2] = F[3] * F[7] - F[6] * F[4];   dJdF[5] = F[6] * F[1] - F[0] * F[7];   dJdF[8] = F[0] * F[4] - F[3] * F[1];

			qeal alpha = 1.0 + (3.0 * lameMu) / (4.0 * lameLamda);
			qeal a0 = lameMu * (1.0 - 1.0 / (Ic + 1.0));
			qeal a3 = lameLamda * (J - alpha);

			for (int i = 0; i < 9; i++)
				F[i] = a0 * F[i] + a3 * dJdF[i];
		}
	}

	__global__ void computeTetElementInternalForce
	(
		int* devTetElementNum,
		qeal* devTetElementInvDm,
		qeal* devTetElementdFPK,
		qeal* devTetElementAttri,
		qeal* devTetElementInternalForce
	)
	{
		__shared__ int elementNum;
		__shared__ qeal FPK[576]; // 9 * 64
		__shared__ qeal invDm[576]; // 9 * 64
		__shared__ qeal volume[64]; // 64
		const int length = gridDim.x *  blockDim.x;
		int idx = threadIdx.x % 9;
		int eid = (threadIdx.x - idx) / 9;
		int r = idx % 3;
		int c = (idx - r) / 3;
		int eleId = (blockIdx.x  * blockDim.x) / 9 + eid;

		if (threadIdx.x == 0)
		{
			elementNum = *devTetElementNum;
		}
		FPK[9 * eid + idx] = devTetElementdFPK[9 * eleId + idx];
		invDm[9 * eid + idx] = devTetElementInvDm[9 * eleId + idx];
		if (idx == 0)
			volume[eid] = devTetElementAttri[3 * eleId];
		__syncthreads();

		for (; eleId < elementNum; eleId += length)
		{
			devTetElementInternalForce[12 * eleId + idx] = volume[eid] * (FPK[9 * eid + r] * invDm[9 * eid + c] + FPK[9 * eid + r + 3] * invDm[9 * eid + c + 3] + FPK[9 * eid + 6 + r] * invDm[9 * eid + 6 + c]);
			__syncthreads();
			if (idx < 3)
			{
				devTetElementInternalForce[12 * eleId + idx + 9] = -1.0 * (devTetElementInternalForce[12 * eleId + idx] + devTetElementInternalForce[12 * eleId + idx + 3] + devTetElementInternalForce[12 * eleId + idx + 6]);
			}
		}
	}

	__global__ void computeTetElementStiffnessMatrix
	(
		int* devTetElementNum,
		qeal* devTetElementdFdu,
		qeal* devTetElementdPdF,
		qeal* devTetElementAttri,
		qeal* devTetElementStiffness
	)
	{
		__shared__ int elementNum;
		__shared__ qeal dPdu[108];
		__shared__ qeal dFdu[108];
		__shared__ qeal dPdF[81];
		__shared__ qeal volume; // 128

		int ridx = threadIdx.x % 3;
		int cidx = threadIdx.y % 3;
		int r = (threadIdx.x - ridx) / 3;
		int c = (threadIdx.y - cidx) / 3;

		if (threadIdx.x == 0)
		{
			elementNum = *devTetElementNum;
			volume = devTetElementAttri[3 * blockIdx.x];
		}

		if (threadIdx.x < 9)
		{
			dFdu[9 * threadIdx.y + threadIdx.x] = devTetElementdFdu[108 * blockIdx.x + 9 * threadIdx.y + threadIdx.x];
			if (threadIdx.y < 9)
				dPdF[9 * threadIdx.y + threadIdx.x] = devTetElementdPdF[81 * blockIdx.x + 9 * threadIdx.y + threadIdx.x];
		}

		__syncthreads();


		if (threadIdx.x < 9)
		{
			dPdu[9 * threadIdx.y + threadIdx.x] = dPdF[threadIdx.x] * dFdu[9 * threadIdx.y] +
				dPdF[9 + threadIdx.x] * dFdu[9 * threadIdx.y + 1] +
				dPdF[18 + threadIdx.x] * dFdu[9 * threadIdx.y + 2] +
				dPdF[27 + threadIdx.x] * dFdu[9 * threadIdx.y + 3] +
				dPdF[36 + threadIdx.x] * dFdu[9 * threadIdx.y + 4] +
				dPdF[45 + threadIdx.x] * dFdu[9 * threadIdx.y + 5] +
				dPdF[54 + threadIdx.x] * dFdu[9 * threadIdx.y + 6] +
				dPdF[63 + threadIdx.x] * dFdu[9 * threadIdx.y + 7] +
				dPdF[72 + threadIdx.x] * dFdu[9 * threadIdx.y + 8];
		}
		__syncthreads();

		qeal k = volume *
			(dPdu[9 * threadIdx.x] * dFdu[9 * threadIdx.y] +
				dPdu[9 * threadIdx.x + 1] * dFdu[9 * threadIdx.y + 1] +
				dPdu[9 * threadIdx.x + 2] * dFdu[9 * threadIdx.y + 2] +
				dPdu[9 * threadIdx.x + 3] * dFdu[9 * threadIdx.y + 3] +
				dPdu[9 * threadIdx.x + 4] * dFdu[9 * threadIdx.y + 4] +
				dPdu[9 * threadIdx.x + 5] * dFdu[9 * threadIdx.y + 5] +
				dPdu[9 * threadIdx.x + 6] * dFdu[9 * threadIdx.y + 6] +
				dPdu[9 * threadIdx.x + 7] * dFdu[9 * threadIdx.y + 7] +
				dPdu[9 * threadIdx.x + 8] * dFdu[9 * threadIdx.y + 8]);
		devTetElementStiffness[144 * blockIdx.x + 12 * threadIdx.y + threadIdx.x] = k;
	}

	__global__ void computeTetElementDeformationGradient
	(
		int* devTetElementNum,
		qeal* devTetElementDisplacement,
		qeal* devTetElementDm,
		qeal* devTetElementInvDm,
		qeal* devTetElementdF
	)
	{
		__shared__ int elementNum;
		const int length = gridDim.x *  blockDim.x;
		int eleId = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			elementNum = *devTetElementNum;
		}
		__syncthreads();

		for (; eleId < elementNum; eleId += length)
		{
			qeal* eleDisplacement = devTetElementDisplacement + 12 * eleId;
			qeal* Dm = devTetElementDm + 9 * eleId;
			qeal* invDm = devTetElementInvDm + 9 * eleId;
			qeal* F = devTetElementdF + 9 * eleId;

			qeal Ds[9];
			computeElementDs(eleDisplacement, Dm, Ds);
			getMutilMatrix(Ds, invDm, F);
		}
	}



	__host__ void computeTetElementStiffness
	(
		int hostTetElementNum,
		int* devTetElementNum,
		qeal* devTetElementDisplacement,
		qeal* devTetElementDm,
		qeal* devTetElementInvDm,
		qeal* devTetElementdFdu,
		qeal* devTetElementFPK, // as eigen value of dPdF
		qeal* devTetElementdPdF, // as eigen vector of dPdF
		qeal* devTetElementAttri,
		qeal* devTetElementStiffness
	)
	{
		dim3 blockSize(THREADS_NUM_32);
		uint32_t size = hostTetElementNum;
		uint32_t num_block = (size + (THREADS_NUM_32 - 1)) / THREADS_NUM_32;
		dim3 gridSize(num_block);

		computeTetElementdPdFEigenSystem << <gridSize, blockSize >> >
			(
				devTetElementNum,
				devTetElementDisplacement,
				devTetElementDm,
				devTetElementInvDm,
				devTetElementAttri,
				devTetElementFPK,
				devTetElementdPdF
				);
		cudaDeviceSynchronize();

		blockSize = dim3(9, 9);
		gridSize = dim3(hostTetElementNum);

		SPDdPdF << <gridSize, blockSize >> >
			(
				devTetElementNum,
				devTetElementFPK,
				devTetElementdPdF,
				devTetElementdPdF
				);
		cudaDeviceSynchronize();

		blockSize = dim3(12, 12);
		gridSize = dim3(hostTetElementNum);

		computeTetElementStiffnessMatrix << <gridSize, blockSize >> >
			(
				devTetElementNum,
				devTetElementdFdu,
				devTetElementdPdF,
				devTetElementAttri,
				devTetElementStiffness
				);
		cudaDeviceSynchronize();
	}

	__global__ void computeTetElementdPdFEigenSystem
	(
		int* devTetElementNum,
		qeal* devTetElementDisplacement,
		qeal* devTetElementDm,
		qeal* devTetElementInvDm,
		qeal* devTetElementAttri,
		qeal* devTetElementEigenValue,
		qeal* devTetElementEigenVector
	)
	{
		__shared__ int elementNum;
		__shared__ qeal eigenValue[288]; // 32 * 9
		__shared__ qeal eigenVector[2592]; // 32 * 81
		__shared__ qeal mF[288];
		__shared__ qeal mU[288];
		__shared__ qeal mV[288];
		__shared__ qeal Singular[288];

		const int length = gridDim.x *  blockDim.x;
		int eleId = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			elementNum = *devTetElementNum;
		}
		__syncthreads();

		for (; eleId < elementNum; eleId += length)
		{
			qeal* eleDisplacement = devTetElementDisplacement + 12 * eleId;
			qeal* Dm = devTetElementDm + 9 * eleId;
			qeal* invDm = devTetElementInvDm + 9 * eleId;
			qeal* vector = devTetElementEigenVector + 81 * eleId;
			qeal* value = devTetElementEigenValue + 9 * eleId;

			qeal* F = mF + 9 * threadIdx.x;
			qeal* U = mU + 9 * threadIdx.x;
			qeal* V = mV + 9 * threadIdx.x;
			qeal* S = Singular + 9 * threadIdx.x;

			qeal Ds[9];
			computeElementDs(eleDisplacement, Dm, Ds);

			getMutilMatrix(Ds, invDm, F);
			//CudaSVD::svd(F, U, S, V);
			CudaSVD::fastSvd3x3(F, U, S, V);
			qeal volume = devTetElementAttri[3 * eleId];
			qeal lameMu = devTetElementAttri[3 * eleId + 1], lameLamda = devTetElementAttri[3 * eleId + 2];
			qeal Ic = F[0] * F[0] + F[3] * F[3] + F[6] * F[6] +
				F[1] * F[1] + F[4] * F[4] + F[7] * F[7] +
				F[2] * F[2] + F[5] * F[5] + F[8] * F[8];

			qeal J = S[0] * S[4] * S[8];

			qeal dJdF[9];
			dJdF[0] = F[4] * F[8] - F[7] * F[5];   dJdF[3] = F[7] * F[2] - F[1] * F[8];	  dJdF[6] = F[1] * F[5] - F[4] * F[2];
			dJdF[1] = F[6] * F[5] - F[3] * F[8];   dJdF[4] = F[0] * F[8] - F[6] * F[2];   dJdF[7] = F[3] * F[2] - F[0] * F[5];
			dJdF[2] = F[3] * F[7] - F[6] * F[4];   dJdF[5] = F[6] * F[1] - F[0] * F[7];   dJdF[8] = F[0] * F[4] - F[3] * F[1];

			qeal alpha = 1.0 + (3.0 * lameMu) / (4.0 * lameLamda);
			qeal a0 = lameMu * (1.0 - 1.0 / (Ic + 1.0));
			qeal a1 = 2.0 * lameMu / ((Ic + 1.0) * (Ic + 1.0));
			qeal a2 = lameLamda;
			qeal a3 = lameLamda * (J - alpha);

			qeal s0 = S[0], s1 = S[4], s2 = S[8];
			qeal s0s0 = s0 * s0, s1s1 = s1 * s1, s2s2 = s2 * s2;
			qeal s0s1 = s0 * s1, s0s2 = s0 * s2, s1s2 = s1 * s2;

			qeal Ax[9];
			Ax[0] = a0 + a1 * s0s0 + a2 * s1s1 * s2s2;
			Ax[4] = a0 + a1 * s1s1 + a2 * s0s0 * s2s2;
			Ax[8] = a0 + a1 * s2s2 + a2 * s0s0 * s1s1;

			Ax[1] = Ax[3] = a1 * s0s1 + a2 * s0s1 * s2s2 + a3 * s2;
			Ax[2] = Ax[6] = a1 * s0s2 + a2 * s0s2 * s1s1 + a3 * s1;
			Ax[5] = Ax[7] = a1 * s1s2 + a2 * s1s2 * s0s0 + a3 * s0;

			qeal AxEvalues[3];
			qeal AxEvectors[9];
			getMatrix3EigenvalueAndEigenvector(3, Ax, AxEvalues, AxEvectors);

			qeal* FEvalues = eigenValue + 9 * threadIdx.x;

			FEvalues[0] = (AxEvalues[0] > 1e-13) ? AxEvalues[0] : 0.0;
			FEvalues[1] = (AxEvalues[1] > 1e-13) ? AxEvalues[1] : 0.0;
			FEvalues[2] = (AxEvalues[2] > 1e-13) ? AxEvalues[2] : 0.0;
			FEvalues[3] = lameLamda * (J - alpha) * s0 + a0;
			FEvalues[3] = (FEvalues[3] > 1e-13) ? FEvalues[3] : 0.0;
			FEvalues[4] = lameLamda * (J - alpha) * s1 + a0;
			FEvalues[4] = (FEvalues[4] > 1e-13) ? FEvalues[4] : 0.0;
			FEvalues[5] = lameLamda * (J - alpha) * s2 + a0;
			FEvalues[5] = (FEvalues[5] > 1e-13) ? FEvalues[5] : 0.0;
			FEvalues[6] = -lameLamda * (J - alpha) * s0 + a0;
			FEvalues[6] = (FEvalues[6] > 1e-13) ? FEvalues[6] : 0.0;
			FEvalues[7] = -lameLamda * (J - alpha) * s1 + a0;
			FEvalues[7] = (FEvalues[7] > 1e-13) ? FEvalues[7] : 0.0;
			FEvalues[8] = -lameLamda * (J - alpha) * s2 + a0;
			FEvalues[8] = (FEvalues[8] > 1e-13) ? FEvalues[8] : 0.0;

			qeal* FEvectors = eigenVector + 81 * threadIdx.x;

			FEvectors[0] = U[0] * AxEvectors[0] * V[0] + U[3] * AxEvectors[1] * V[3] + U[6] * AxEvectors[2] * V[6];
			FEvectors[1] = U[1] * AxEvectors[0] * V[0] + U[4] * AxEvectors[1] * V[3] + U[7] * AxEvectors[2] * V[6];
			FEvectors[2] = U[2] * AxEvectors[0] * V[0] + U[5] * AxEvectors[1] * V[3] + U[8] * AxEvectors[2] * V[6];
			FEvectors[3] = U[0] * AxEvectors[0] * V[1] + U[3] * AxEvectors[1] * V[4] + U[6] * AxEvectors[2] * V[7];
			FEvectors[4] = U[1] * AxEvectors[0] * V[1] + U[4] * AxEvectors[1] * V[4] + U[7] * AxEvectors[2] * V[7];
			FEvectors[5] = U[2] * AxEvectors[0] * V[1] + U[5] * AxEvectors[1] * V[4] + U[8] * AxEvectors[2] * V[7];
			FEvectors[6] = U[0] * AxEvectors[0] * V[2] + U[3] * AxEvectors[1] * V[5] + U[6] * AxEvectors[2] * V[8];
			FEvectors[7] = U[1] * AxEvectors[0] * V[2] + U[4] * AxEvectors[1] * V[5] + U[7] * AxEvectors[2] * V[8];
			FEvectors[8] = U[2] * AxEvectors[0] * V[2] + U[5] * AxEvectors[1] * V[5] + U[8] * AxEvectors[2] * V[8];

			FEvectors[9] = U[0] * AxEvectors[3] * V[0] + U[3] * AxEvectors[4] * V[3] + U[6] * AxEvectors[5] * V[6];
			FEvectors[10] = U[1] * AxEvectors[3] * V[0] + U[4] * AxEvectors[4] * V[3] + U[7] * AxEvectors[5] * V[6];
			FEvectors[11] = U[2] * AxEvectors[3] * V[0] + U[5] * AxEvectors[4] * V[3] + U[8] * AxEvectors[5] * V[6];
			FEvectors[12] = U[0] * AxEvectors[3] * V[1] + U[3] * AxEvectors[4] * V[4] + U[6] * AxEvectors[5] * V[7];
			FEvectors[13] = U[1] * AxEvectors[3] * V[1] + U[4] * AxEvectors[4] * V[4] + U[7] * AxEvectors[5] * V[7];
			FEvectors[14] = U[2] * AxEvectors[3] * V[1] + U[5] * AxEvectors[4] * V[4] + U[8] * AxEvectors[5] * V[7];
			FEvectors[15] = U[0] * AxEvectors[3] * V[2] + U[3] * AxEvectors[4] * V[5] + U[6] * AxEvectors[5] * V[8];
			FEvectors[16] = U[1] * AxEvectors[3] * V[2] + U[4] * AxEvectors[4] * V[5] + U[7] * AxEvectors[5] * V[8];
			FEvectors[17] = U[2] * AxEvectors[3] * V[2] + U[5] * AxEvectors[4] * V[5] + U[8] * AxEvectors[5] * V[8];

			FEvectors[18] = U[0] * AxEvectors[6] * V[0] + U[3] * AxEvectors[7] * V[3] + U[6] * AxEvectors[8] * V[6];
			FEvectors[19] = U[1] * AxEvectors[6] * V[0] + U[4] * AxEvectors[7] * V[3] + U[7] * AxEvectors[8] * V[6];
			FEvectors[20] = U[2] * AxEvectors[6] * V[0] + U[5] * AxEvectors[7] * V[3] + U[8] * AxEvectors[8] * V[6];
			FEvectors[21] = U[0] * AxEvectors[6] * V[1] + U[3] * AxEvectors[7] * V[4] + U[6] * AxEvectors[8] * V[7];
			FEvectors[22] = U[1] * AxEvectors[6] * V[1] + U[4] * AxEvectors[7] * V[4] + U[7] * AxEvectors[8] * V[7];
			FEvectors[23] = U[2] * AxEvectors[6] * V[1] + U[5] * AxEvectors[7] * V[4] + U[8] * AxEvectors[8] * V[7];
			FEvectors[24] = U[0] * AxEvectors[6] * V[2] + U[3] * AxEvectors[7] * V[5] + U[6] * AxEvectors[8] * V[8];
			FEvectors[25] = U[1] * AxEvectors[6] * V[2] + U[4] * AxEvectors[7] * V[5] + U[7] * AxEvectors[8] * V[8];
			FEvectors[26] = U[2] * AxEvectors[6] * V[2] + U[5] * AxEvectors[7] * V[5] + U[8] * AxEvectors[8] * V[8];

			qeal eScalar = 1.0 / CUDA_SQRT(2.0);
#pragma unroll
			for (int i = 0; i < 9; i++)
			{
				FEvectors[3 * 9 + i] = eScalar * (-U[6 + i % 3] * V[3 + (i - (i % 3)) / 3] + U[3 + i % 3] * V[6 + (i - (i % 3)) / 3]);
			}

#pragma unroll
			for (int i = 0; i < 9; i++)
			{
				FEvectors[4 * 9 + i] = eScalar * (-U[6 + i % 3] * V[0 + (i - (i % 3)) / 3] + U[0 + i % 3] * V[6 + (i - (i % 3)) / 3]);
			}

#pragma unroll
			for (int i = 0; i < 9; i++)
			{
				FEvectors[5 * 9 + i] = eScalar * (-U[3 + i % 3] * V[0 + (i - (i % 3)) / 3] + U[0 + i % 3] * V[3 + (i - (i % 3)) / 3]);
			}
#pragma unroll
			for (int i = 0; i < 9; i++)
			{
				FEvectors[6 * 9 + i] = eScalar * (U[6 + i % 3] * V[3 + (i - (i % 3)) / 3] + U[3 + i % 3] * V[6 + (i - (i % 3)) / 3]);
			}
#pragma unroll
			for (int i = 0; i < 9; i++)
			{
				FEvectors[7 * 9 + i] = eScalar * (U[6 + i % 3] * V[0 + (i - (i % 3)) / 3] + U[0 + i % 3] * V[6 + (i - (i % 3)) / 3]);
			}
#pragma unroll
			for (int i = 0; i < 9; i++)
			{
				FEvectors[8 * 9 + i] = eScalar * (U[3 + i % 3] * V[0 + (i - (i % 3)) / 3] + U[0 + i % 3] * V[3 + (i - (i % 3)) / 3]);
			}

			for (int i = 0; i < 81; i++)
				vector[i] = FEvectors[i];

			for (int i = 0; i < 9; i++)
				value[i] = FEvalues[i];
		}
	}

	__global__ void SPDdPdF
	(
		int* devTetElementNum,
		qeal* devTetElementEigenValue,
		qeal* devTetElementEigenVector,
		qeal* devTetElementdPdF
	)
	{
		__shared__ int elementNum;
		__shared__ qeal eigenVector[81];
		__shared__ qeal eigenValue[9];
		__shared__ qeal dPdF[81];

		int ridx = threadIdx.x % 3;
		int cidx = threadIdx.y % 3;
		int r = (threadIdx.x - ridx) / 3;
		int c = (threadIdx.y - cidx) / 3;

		if (threadIdx.x == 0)
		{
			elementNum = *devTetElementNum;
		}

		if (threadIdx.x == threadIdx.y)
			eigenValue[threadIdx.x] = devTetElementEigenValue[9 * blockIdx.x + threadIdx.x];
		eigenVector[9 * threadIdx.y + threadIdx.x] = devTetElementEigenVector[81 * blockIdx.x + 9 * threadIdx.y + threadIdx.x];
		__syncthreads();

		dPdF[9 * threadIdx.y + threadIdx.x] = 0;
		for (int i = 0; i < 9; i++)
			dPdF[9 * threadIdx.y + threadIdx.x] += eigenVector[9 * i + threadIdx.x] * eigenValue[i] * eigenVector[9 * i + threadIdx.y];

		devTetElementdPdF[81 * blockIdx.x + 9 * threadIdx.y + threadIdx.x] = dPdF[9 * threadIdx.y + threadIdx.x];
	}



	__device__ __forceinline__
		void computeU12x12q(qeal * Mat12x12, qeal * q, qeal* result)
	{
		for (int i = 0; i < 12; i++)
		{
			result[i] += Mat12x12[0 + i] * q[0] + Mat12x12[12 + i] * q[1] + Mat12x12[24 + i] * q[2] + Mat12x12[36 + i] * q[3] + Mat12x12[48 + i] * q[4] + Mat12x12[60 + i] * q[5] + Mat12x12[72 + i] * q[6] + Mat12x12[84 + i] * q[7] + Mat12x12[96 + i] * q[8] + Mat12x12[108 + i] * q[9] + Mat12x12[120 + i] * q[10] + Mat12x12[132 + i] * q[11];
		}
	}

	__device__ __forceinline__
		void computeUT12x12q(qeal * Mat12x12, qeal * q, qeal* result)
	{
		for (int i = 0; i < 12; i++)
		{
			result[i] = Mat12x12[12 * i] * q[0] + Mat12x12[12 * i + 1] * q[1] + Mat12x12[12 * i + 2] * q[2] + Mat12x12[12 * i + 3] * q[3] + Mat12x12[12 * i + 4] * q[4] + Mat12x12[12 * i + 5] * q[5] + Mat12x12[12 * i + 6] * q[6] + Mat12x12[12 * i + 7] * q[7] + Mat12x12[12 * i + 8] * q[8] + Mat12x12[12 * i + 9] * q[9] + Mat12x12[12 * i + 10] * q[10] + Mat12x12[12 * i + 11] * q[11];
		}
	}

	__device__ __forceinline__
		void computeElementDs(qeal* displacement, qeal* Dm, qeal* result)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				result[j * 3 + i] = displacement[3 * j + i] - displacement[3 * 3 + i];
			}
		}

		for (int i = 0; i < 9; i++)
			result[i] += Dm[i];
	}

	__device__ __forceinline__
		void computeStableNeoHookean(int eleId, qeal * F, qeal * FPK, qeal * PF, qeal & lameMu, qeal & lameLamda)
	{
		qeal Ic = F[0] * F[0] + F[3] * F[3] + F[6] * F[6] +
			F[1] * F[1] + F[4] * F[4] + F[7] * F[7] +
			F[2] * F[2] + F[5] * F[5] + F[8] * F[8];

		qeal J = getMatrix3Determinant(F);

		qeal dJdF[9];
		dJdF[0] = F[4] * F[8] - F[7] * F[5];   dJdF[3] = F[7] * F[2] - F[1] * F[8];	  dJdF[6] = F[1] * F[5] - F[4] * F[2];
		dJdF[1] = F[6] * F[5] - F[3] * F[8];   dJdF[4] = F[0] * F[8] - F[6] * F[2];   dJdF[7] = F[3] * F[2] - F[0] * F[5];
		dJdF[2] = F[3] * F[7] - F[6] * F[4];   dJdF[5] = F[6] * F[1] - F[0] * F[7];   dJdF[8] = F[0] * F[4] - F[3] * F[1];

		qeal alpha = 1.0 + (3.0 * lameMu) / (4.0 * lameLamda);
		qeal a0 = lameMu * (1.0 - 1.0 / (Ic + 1.0));
		qeal a1 = 2.0 * lameMu / ((Ic + 1.0) * (Ic + 1.0));
		qeal a2 = lameLamda;
		qeal a3 = lameLamda * (J - alpha);

		addMatrix3(F, dJdF, FPK, a0, a3);

	}


	__host__ void assembleReducedStiffness
	(
		int hostAssembleBlockIndexNum,
		int* devAssembleBlockIndexNum,
		int* devAssembleTask,
		int* devStiffnessBlockSharedTetElementList,
		int* devAssembleTaskSharedTetElementNum,
		int* devAssembleTaskSharedTetElementoffset,
		int* devPojectionStiffnessList,
		int* devTetElementSharedFrameList,
		int* devTetElementSharedFrameOffset,
		qeal* devTetElementFrameProjectionBuffer,
		int* devTetElementFrameProjectionNum,
		int* devTetElementFrameProjectionOffset,
		qeal* devTetElementStiffness,
		int* devReducedDim,
		qeal* devReducedStiffness
	)
	{
		dim3 blockSize(12, 12);
		dim3 gridSize(hostAssembleBlockIndexNum);

		assembleReducedStiffness << <gridSize, blockSize >> >
			(
				devAssembleBlockIndexNum,
				devAssembleTask,
				devStiffnessBlockSharedTetElementList,
				devAssembleTaskSharedTetElementNum,
				devAssembleTaskSharedTetElementoffset,
				devPojectionStiffnessList,
				devTetElementSharedFrameList,
				devTetElementSharedFrameOffset,
				devTetElementFrameProjectionBuffer,
				devTetElementFrameProjectionNum,
				devTetElementFrameProjectionOffset,
				devTetElementStiffness,
				devReducedDim,
				devReducedStiffness
				);
		cudaDeviceSynchronize();
	}

	__global__ void assembleReducedStiffness
	(
		int* devAssembleBlockIndexNum,
		int* devAssembleTask,
		int* devStiffnessBlockSharedTetElementList,
		int* devAssembleTaskSharedTetElementNum,
		int* devAssembleTaskSharedTetElementoffset,
		int* devPojectionStiffnessList,
		int* devTetElementSharedFrameList,
		int* devTetElementSharedFrameOffset,
		qeal* devTetElementFrameProjectionBuffer,
		int* devTetElementFrameProjectionNum,
		int* devTetElementFrameProjectionOffset,
		qeal* devTetElementStiffness,
		int* devReducedDim,
		qeal* devReducedStiffness
	)
	{
		__shared__ int sharedInteger[6];
		//
		__shared__ qeal reducedStiffness[144];
		__shared__ qeal SP[144];
		__shared__ qeal S[144];
		__shared__ qeal Uit[16];
		__shared__ qeal Uj[16];

		if (threadIdx.x == 0 && threadIdx.y == 0)
		{
			sharedInteger[0] = *devAssembleBlockIndexNum; // blockNum
			sharedInteger[1] = *devReducedDim; // reduced dim
			sharedInteger[2] = 12 * devAssembleTask[2 * blockIdx.x]; // iFrameId * 12
			sharedInteger[3] = 12 * devAssembleTask[2 * blockIdx.x + 1]; // jFrameId * 12
			sharedInteger[4] = devAssembleTaskSharedTetElementoffset[blockIdx.x]; // offset
			sharedInteger[5] = devAssembleTaskSharedTetElementNum[blockIdx.x]; // num;
		}
		reducedStiffness[12 * threadIdx.y + threadIdx.x] = 0;
		__syncthreads();

#pragma unroll
		for (int i = 0; i < sharedInteger[5]; i++)
		{
			int index = devStiffnessBlockSharedTetElementList[sharedInteger[4] + i];

			int eleId = devPojectionStiffnessList[3 * index];
			int UitIndex = devPojectionStiffnessList[3 * index + 1];
			int UjIndex = devPojectionStiffnessList[3 * index + 2];
			int offset = devTetElementSharedFrameOffset[eleId];
			qeal* localUi = devTetElementFrameProjectionBuffer + 16 * (offset + UitIndex);
			qeal* localUj = devTetElementFrameProjectionBuffer + 16 * (offset + UjIndex);
			qeal* stiffness = devTetElementStiffness + 144 * eleId;

			int ridx = threadIdx.x % 3;
			int cidx = threadIdx.y % 3;
			int r = (threadIdx.x - ridx) / 3;
			int c = (threadIdx.y - cidx) / 3;
			__syncthreads();

			S[12 * threadIdx.y + threadIdx.x] = stiffness[12 * threadIdx.y + threadIdx.x];
			__syncthreads();

			if (ridx == 0 && cidx == 0)
			{
				Uj[4 * c + r] = localUj[4 * c + r];
				Uit[4 * r + c] = localUi[4 * c + r];
			}
			__syncthreads();

			SP[12 * threadIdx.y + threadIdx.x] = S[12 * (3 * 0 + cidx) + threadIdx.x] * Uj[4 * c + 0] + S[12 * (3 * 1 + cidx) + threadIdx.x] * Uj[4 * c + 1] + S[12 * (3 * 2 + cidx) + threadIdx.x] * Uj[4 * c + 2] + S[12 * (3 * 3 + cidx) + threadIdx.x] * Uj[4 * c + 3];
			__syncthreads();

			reducedStiffness[12 * threadIdx.y + threadIdx.x] += Uit[4 * 0 + r] * SP[12 * (3 * c + cidx) + 3 * 0 + ridx] + Uit[4 * 1 + r] * SP[12 * (3 * c + cidx) + 3 * 1 + ridx] + Uit[4 * 2 + r] * SP[12 * (3 * c + cidx) + 3 * 2 + ridx] + Uit[4 * 3 + r] * SP[12 * (3 * c + cidx) + 3 * 3 + ridx];
		}
		__syncthreads();

		qeal* stiffness = devReducedStiffness + (sharedInteger[3] + threadIdx.y) * sharedInteger[1] + sharedInteger[2] + threadIdx.x;
		stiffness[0] = reducedStiffness[12 * threadIdx.y + threadIdx.x];
		if (sharedInteger[2] != sharedInteger[3])
		{
			stiffness = devReducedStiffness + (sharedInteger[2] + threadIdx.x) * sharedInteger[1] + sharedInteger[3] + threadIdx.y;
			stiffness[0] = reducedStiffness[12 * threadIdx.y + threadIdx.x];
		}
	}


	__host__ void computeMedialPointsMovingDir
	(
		int hostMedialPointsNum,
		int* devMedialPointsNum,
		qeal* devMedialOriPointPosition,
		qeal* devReducedXtilde,
		qeal* devMedialPointMovingDir
	)
	{
		dim3 blockSize(THREADS_NUM);
		uint32_t num_block = (hostMedialPointsNum + (THREADS_NUM - 1)) / THREADS_NUM;
		dim3 gridSize(num_block);
		computeMedialPointsMovingDir << <gridSize, blockSize >> >
			(
				devMedialPointsNum,
				devMedialOriPointPosition,
				devReducedXtilde,
				devMedialPointMovingDir
				);
		cudaDeviceSynchronize();
	}

	__global__ void computeMedialPointsMovingDir
	(
		int* devMedialPointsNum,
		qeal* devMedialOriPointPosition,
		qeal* devReducedXtilde,
		qeal* devMedialPointMovingDir
	)
	{
		__shared__ qeal oriXYZ[1536]; // 3 * THREADS_NUM
		__shared__ int num;
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			num = *devMedialPointsNum;
		}
		__syncthreads();
		if (tid < num)
		{
			oriXYZ[3 * threadIdx.x] = devMedialOriPointPosition[3 * tid]; // x
			oriXYZ[3 * threadIdx.x + 1] = devMedialOriPointPosition[3 * tid + 1]; //y
			oriXYZ[3 * threadIdx.x + 2] = devMedialOriPointPosition[3 * tid + 2]; // z
		}
		__syncthreads();
		for (; tid < num; tid += length)
		{
			devMedialPointMovingDir[3 * tid] = oriXYZ[3 * threadIdx.x] * devReducedXtilde[12 * tid] + oriXYZ[3 * threadIdx.x + 1] * devReducedXtilde[12 * tid + 3] + oriXYZ[3 * threadIdx.x + 2] * devReducedXtilde[12 * tid + 6] + devReducedXtilde[12 * tid + 9];

			devMedialPointMovingDir[3 * tid + 1] = oriXYZ[3 * threadIdx.x] * devReducedXtilde[12 * tid + 1] + oriXYZ[3 * threadIdx.x + 1] * devReducedXtilde[12 * tid + 4] + oriXYZ[3 * threadIdx.x + 2] * devReducedXtilde[12 * tid + 7] + devReducedXtilde[12 * tid + 10];

			devMedialPointMovingDir[3 * tid + 2] = oriXYZ[3 * threadIdx.x] * devReducedXtilde[12 * tid + 2] + oriXYZ[3 * threadIdx.x + 1] * devReducedXtilde[12 * tid + 5] + oriXYZ[3 * threadIdx.x + 2] * devReducedXtilde[12 * tid + 8] + devReducedXtilde[12 * tid + 11];
		}

	}

	__host__ void fillTetElementFrameProjectionBuffer
	(
		int hostTetElementNum,
		int* devTetElementNum,
		int* devTetElementIndices,
		qeal* devTetPoints,
		qeal* devTetElementFrameWeightList,
		int* devTetElementFrameProjectionNum,
		int* devTetElementFrameProjectionOffset,
		qeal* devTetElementFrameProjectionBuffer
	)
	{
		dim3 blockSize(THREADS_NUM);
		uint32_t num_block = (hostTetElementNum + (THREADS_NUM - 1)) / THREADS_NUM;
		dim3 gridSize(num_block);
		fillTetElementFrameProjectionBuffer << <gridSize, blockSize >> >
			(
				devTetElementNum,
				devTetElementIndices,
				devTetPoints,
				devTetElementFrameWeightList,
				devTetElementFrameProjectionNum,
				devTetElementFrameProjectionOffset,
				devTetElementFrameProjectionBuffer
				);
		CUDA_CALL(cudaDeviceSynchronize());
	}

	__global__ void fillTetElementFrameProjectionBuffer
	(
		int* devTetElementNum,
		int* devTetElementIndices,
		qeal* devTetPoints,
		qeal* devTetElementFrameWeightList,
		int* devTetElementFrameProjectionNum,
		int* devTetElementFrameProjectionOffset,
		qeal* devTetElementFrameProjectionBuffer
	)
	{
		__shared__ int num;
		const int length = gridDim.x *  blockDim.x;
		int tid = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
		{
			num = *devTetElementNum;
		}
		__syncthreads();

		for (; tid < num; tid += length)
		{
			int offset = devTetElementFrameProjectionOffset[tid];
			int frameNum = devTetElementFrameProjectionNum[tid];
			qeal* weight = devTetElementFrameWeightList + 4 * offset;
			qeal* buffer = devTetElementFrameProjectionBuffer + 16 * offset;

			int v0 = devTetElementIndices[4 * tid];
			int v1 = devTetElementIndices[4 * tid + 1];
			int v2 = devTetElementIndices[4 * tid + 2];
			int v3 = devTetElementIndices[4 * tid + 3];

			qeal v0x = devTetPoints[3 * v0];
			qeal v0y = devTetPoints[3 * v0 + 1];
			qeal v0z = devTetPoints[3 * v0 + 2];

			qeal v1x = devTetPoints[3 * v1];
			qeal v1y = devTetPoints[3 * v1 + 1];
			qeal v1z = devTetPoints[3 * v1 + 2];

			qeal v2x = devTetPoints[3 * v2];
			qeal v2y = devTetPoints[3 * v2 + 1];
			qeal v2z = devTetPoints[3 * v2 + 2];

			qeal v3x = devTetPoints[3 * v3];
			qeal v3y = devTetPoints[3 * v3 + 1];
			qeal v3z = devTetPoints[3 * v3 + 2];



			for (int i = 0; i < frameNum; i++)
			{
				qeal w0 = weight[4 * i + 0];
				qeal w1 = weight[4 * i + 1];
				qeal w2 = weight[4 * i + 2];
				qeal w3 = weight[4 * i + 3];

				buffer[16 * i + 0] = w0 * v0x;
				buffer[16 * i + 1] = w1 * v1x;
				buffer[16 * i + 2] = w2 * v2x;
				buffer[16 * i + 3] = w3 * v3x;

				buffer[16 * i + 4] = w0 * v0y;
				buffer[16 * i + 5] = w1 * v1y;
				buffer[16 * i + 6] = w2 * v2y;
				buffer[16 * i + 7] = w3 * v3y;

				buffer[16 * i + 8] = w0 * v0z;
				buffer[16 * i + 9] = w1 * v1z;
				buffer[16 * i + 10] = w2 * v2z;
				buffer[16 * i + 11] = w3 * v3z;

				buffer[16 * i + 12] = w0;
				buffer[16 * i + 13] = w1;
				buffer[16 * i + 14] = w2;
				buffer[16 * i + 15] = w3;
			}
		}
	}

}

