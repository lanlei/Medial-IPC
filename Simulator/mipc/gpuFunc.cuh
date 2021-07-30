#ifndef MIPC_CUDA_FUNC_CUH
#define MIPC_CUDA_FUNC_CUH
#include "Simulator\Cuda\CudaHeader.cuh"

namespace MIPC
{	
	__device__ __forceinline__
		void computeElementDs(qeal* displacement, qeal* Dm, qeal* result);
	
	__host__ void updatePredictivePos
	(
		int hostDim,
		int* devDim,
		qeal* devXn,
		qeal* devV,
		qeal* devTimeStep,
		qeal* devXtilde
	);

	__global__ void updatePredictivePos
	(
		int* devDim,
		qeal* devXn,
		qeal* devV,
		qeal* devTimeStep,
		qeal* devXtilde
	);

	__host__ void updateLineSearchX
	(
		int hostDim,
		int* devDim,
		qeal* devXn,
		qeal* devDir,
		qeal alpha,
		qeal* devX
	);

	__global__ void updateLineSearchX
	(
		int* devDim,
		qeal* devXn,
		qeal* devDir,
		qeal* devAlpha,
		qeal* devX
	);

	__host__ void computeInertialDiffX
	(
		int hostDim,
		int* devDim,
		qeal* devX,
		qeal* devXtilde,
		qeal* devTetPointsMass,
		qeal* devDiffX,
		qeal* devMassMulDiffX
	);

	__global__ void computeInertialDiffX
	(
		int* devDim,
		qeal* devX,
		qeal* devXtilde,
		qeal* devTetPointsMass,
		qeal* devDiffX,
		qeal* devMassMulDiffX
	);

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
	);

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
	);

	__host__ void computeReducedInertia
	(
		int hostDim,
		int* devDim,
		qeal* devReducedX,
		qeal* devReducedXn,
		qeal* devReducedVn,
		qeal* devTimeStep,
		qeal* devInertia
	);


	__global__ void computeReducedInertia
	(
		int* devDim,
		qeal* devReducedX,
		qeal* devReducedXn,
		qeal* devReducedVn,
		qeal* devTimeStep,
		qeal* devInertia
	);

	__host__ void updatedVelocity
	(
		int hostDim,
		int* devDim,
		qeal* devX,
		qeal* devXtilde,
		qeal* devTimeStep,
		qeal* devVelocity
	);

	__global__ void updatedVelocity
	(
		int* devDim,
		qeal* devX,
		qeal* devXtilde,
		qeal* devTimeStep,
		qeal* devVelocity
	);

	__host__ void assembleTetELementX
	(
		int hostTetElementNum,
		int* devTetElementNum,
		int*devTetElementIndices,
		qeal* devTetPointsX,
		qeal* devTetElementX
	);

	__global__ void assembleTetELementX
	(
		int* devTetElementNum,
		int*devTetElementIndices,
		qeal* devTetPointsX,
		qeal* devTetElementX
	);

	__host__ void assembleTetPointsForceFromElementForce
	(
		int hostTetPointsNum,
		int* devTetPointsNum,
		int* devTetPointsSharedElementNum,
		int* devTetPointsSharedElementOffset,
		int* devTetPointsSharedElementList,
		qeal* devTetElementForce,
		qeal* devTetPointInternalForce
	);

	__global__ void assembleTetPointsForceFromElementForce
	(
		int* devTetPointsNum,
		int* devTetPointsSharedElementNum,
		int* devTetPointsSharedElementOffset,
		int* devTetPointsSharedElementList,
		qeal* devTetElementForce,
		qeal* devTetPointInternalForce
	);

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
	);

	__global__ void computeTetElementFPK
	(
		int* devTetElementNum,
		qeal* devTetElementDisplacement,
		qeal* devTetElementDm,
		qeal* devTetElementInvDm,
		qeal* devTetElementAttri,
		qeal* devTetElementdFPK
	);

	__global__ void computeTetElementInternalForce
	(
		int* devTetElementNum,
		qeal* devTetElementInvDm,
		qeal* devTetElementdFPK,
		qeal* devTetElementAttri,
		qeal* devTetElementInternalForce
	);

	__global__ void computeTetElementStiffnessMatrix
	(
		int* devTetElementNum,
		qeal* devTetElementdFdu,
		qeal* devTetElementdPdF,
		qeal* devTetElementAttri,
		qeal* devTetElementStiffness
	);

	__global__ void computeTetElementDeformationGradient
	(
		int* devTetElementNum,
		qeal* devTetElementDisplacement,
		qeal* devTetElementDm,
		qeal* devTetElementInvDm,
		qeal* devTetElementdFPK
	);

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
	);


	//
	__global__ void computeTetElementdPdFEigenSystem
	(
		int* devTetElementNum,
		qeal* devTetElementDisplacement,
		qeal* devTetElementDm,
		qeal* devTetElementInvDm,
		qeal* devTetElementAttri,
		qeal* devTetElementEigenValue,
		qeal* devTetElementEigenVector
	);

	__global__ void SPDdPdF
	(
		int* devTetElementNum,
		qeal* devTetElementEigenValue,
		qeal* devTetElementEigenVector,
		qeal* devTetElementdPdF
	);


	__device__ __forceinline__
		void computeU12x12q(qeal* Mat12x12, qeal* q, qeal* result);

	__device__ __forceinline__
		void computeUT12x12q(qeal* Mat12x12, qeal* q, qeal* result);

	__device__ __forceinline__
		void computeElementDs(qeal* displacement, qeal* Dm, qeal* result);

	__device__ __forceinline__
		void computeStableNeoHookean(int eleid, qeal* F, qeal* FPK, qeal* PF, qeal& lameMiu, qeal& lameLamda);

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
	);

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
	);

	__host__ void computeMedialPointsMovingDir
	(
		int hostMedialPointsNum,
		int* devMedialPointsNum,
		qeal* devMedialOriPointPosition,
		qeal* devReducedXtilde,
		qeal* devMedialPointMovingDir
	);

	__global__ void computeMedialPointsMovingDir
	(
		int* devMedialPointsNum,
		qeal* devMedialOriPointPosition,
		qeal* devReducedXtilde,
		qeal* devMedialPointMovingDir
	);

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
	);

	__global__ void fillTetElementFrameProjectionBuffer
	(
		int* devTetElementNum,
		int* devTetElementIndices,
		qeal* devTetPoints,
		qeal* devTetElementFrameWeightList,
		int* devTetElementFrameProjectionNum,
		int* devTetElementFrameProjectionOffset,
		qeal* devTetElementFrameProjectionBuffer
	);
}



#endif