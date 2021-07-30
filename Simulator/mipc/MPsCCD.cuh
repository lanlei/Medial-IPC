#ifndef Medial_Primitives_CCD_CUH
#define Medial_Primitives_CCD_CUH
#include "Simulator\Cuda\CudaHeader.cuh"
#include "Simulator\Cuda\CudaMatrixOperator.cuh"

namespace MIPC
{
#define COLLISION_CC 0
#define COLLISION_SS 1
#define COLLISION_DEFORMABLE_WITH_STATIC_CC 2
#define COLLISION_DEFORMABLE_WITH_STATIC_SS 3
#define COLLISION_STATIC_WITH_DEFORMABLE_CC 4
#define COLLISION_STATIC_WITH_DEFORMABLE_SS 5

#define RANGE_FLAG 0
#define RANGE_MIN_INDEX 1
#define RANGE_MAX_INDEX 2

	__host__ void MPsCCD
	(
		int hostCollisionEventNum,
		int* devCollisionEventNum,
		qeal* devMedialPointPosition,
		qeal* devMedialPointRadius,
		qeal* devStaticMedialPointPosition,
		qeal* devStaticMedialPointRadius,
		qeal* devMedialPointMovingDir,
		int* devCollisionEventList,
		qeal* devCCD
	);

	__global__ void MPsCCD
	(
		int* devCollisionEventNum,
		qeal* devMedialPointPosition,
		qeal* devMedialPointRadius,
		qeal* devStaticMedialPointPosition,
		qeal* devStaticMedialPointRadius,
		qeal* devMedialPointMovingDir,
		int* devCollisionEventList,
		qeal* devCCD
	);

	__device__ __forceinline__
		bool dcd(qeal* C1, qeal* V1, qeal* C2, qeal* V2, qeal* C3, qeal* V3, qeal& R1, qeal& R2, qeal& R3, bool is_ss, qeal norrow);

	__device__ __forceinline__
		bool dcd(qeal* C1, qeal* C2, qeal* C3, qeal& R1, qeal& R2, qeal& R3, bool is_ss, qeal norrow);

	__device__ __forceinline__
		bool dcd(qeal A, qeal B, qeal C, qeal D, qeal E, qeal F, bool is_ss, qeal norrow);

	__device__ __forceinline__
		bool ccd(qeal * C1, qeal * V1, qeal * C2, qeal * V2, qeal * C3, qeal * V3, qeal& R1, qeal& R2, qeal& R3, bool ss, qeal& ftc);

	__device__ __forceinline__
		void computeCollisionDistance(qeal * C1, qeal * C2, qeal * C3, qeal& R1, qeal& R2, qeal& R3, bool ss, qeal& dist, qeal& alpha, qeal& beta);

	__device__ __forceinline__
		bool checkEndpointAlphaBetaCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow);

	__device__ __forceinline__
		bool checkAlphaIsZeroCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow);

	__device__ __forceinline__
		bool checkBetaIsZeroCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow);

	__device__ __forceinline__
		bool checkAlphaIsOneCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow);

	__device__ __forceinline__
		bool checkBetaIsOneCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow);

	__device__ __forceinline__
		bool checkAlphaPlusBetaIsOneCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow);

	__device__ __forceinline__
		bool checkAlphaBetaCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, int collideType, bool is_ss, qeal norrow);

	__device__ __forceinline__
		qeal valueOfQuadircSurface2D(qeal x, qeal y, qeal& A, qeal& B, qeal& C, qeal& D, qeal& E, qeal& F);

	__device__ __forceinline__
		qeal computeQuadricEquation(qeal x, qeal& a, qeal& b, qeal& c);

	__device__ __forceinline__
		qeal computeQuarticEquation(qeal x, qeal& a, qeal& b, qeal& c, qeal& d, qeal& e);

	__device__ __forceinline__
		qeal computeSexticEquation(qeal x, qeal& a, qeal& b, qeal& c, qeal& d, qeal& e, qeal& f, qeal& g);

	__device__ __forceinline__
		void solveQuadricNEq(qeal& a, qeal& b, qeal& c, qeal * x, int& countRoot, qeal* posRange, qeal* negRange, qeal mini = 0.0, qeal maxi = 1.0);

	__device__ __forceinline__
		bool solveQuadricEquation(qeal& a, qeal& b, qeal& c, qeal* x, int& countRoot, bool banSame, qeal mini, qeal maxi);

	__device__ __forceinline__
		void solveQuarticNEq(qeal&  a, qeal& b, qeal& c, qeal& d, qeal& e, qeal * x, int& countRoot, qeal* posRange, qeal* negRange, qeal mini = 0.0, qeal maxi = 1.0);

	__device__ __forceinline__
		bool solveCubicEquation(qeal& a, qeal& b, qeal& c, qeal& d, qeal* x, int& countRoot, bool banSame, qeal mini, qeal maxi);

	__device__ __forceinline__
		bool solveQuarticEquation(qeal& a, qeal& b, qeal& c, qeal& d, qeal& e, qeal* x, int& countRoot, bool banSame, qeal mini, qeal maxi);

	__device__ __forceinline__
		bool overlapRange(qeal* a, qeal mini, qeal maxi);

	__device__ __forceinline__
		void genQuarticCoeffs(qeal& P_4, qeal& P_3, qeal& P_2, qeal& P_1, qeal& P_0, qeal& A_2, qeal& A_1, qeal& A_0, qeal& B_2, qeal& B_1, qeal& B_0, qeal& C_2, qeal& C_1, qeal& C_0, qeal& D_2, qeal& D_1, qeal& D_0, qeal S = 1.0);

	__device__ __forceinline__
		void genSexticCoeffs(qeal& P_6, qeal& P_5, qeal& P_4, qeal& P_3, qeal& P_2, qeal& P_1, qeal& P_0, qeal& A_2, qeal& A_1, qeal& A_0, qeal& B_2, qeal& B_1, qeal& B_0, qeal& C_2, qeal& C_1, qeal& C_0, qeal S = 1.0);

	__device__ __forceinline__
		qeal SolveSexticEqForFTC(qeal& W_6, qeal& W_5, qeal& W_4, qeal& W_3, qeal& W_2, qeal& W_1, qeal& W_0, int& banRootsNum, qeal* banRoots, qeal* range, qeal norrow);

	__device__ __forceinline__
		qeal NewtonSolverSexticEqForFTC(qeal W_6, qeal W_5, qeal W_4, qeal W_3, qeal W_2, qeal W_1, qeal W_0, qeal* r, qeal norrow);

	__device__ __forceinline__
		qeal NewtonSolverForIntersectionToi(qeal& W_6, qeal& W_5, qeal& W_4, qeal& W_3, qeal& W_2, qeal& W_1, qeal& W_0, qeal mini, qeal maxi);

	__device__ __forceinline__
		qeal VectorDot(qeal* v1, qeal* v2);
}



#endif