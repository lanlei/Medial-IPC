#pragma once
#ifndef COLLISION_DETECTION_ON_MEDIAL_MESH_H
#define COLLISION_DETECTION_ON_MEDIAL_MESH_H
#include "MatrixCore.h"
#include "Commom\PolynomialSolver.h"
#include <queue>

namespace CDMM
{
#define CHECK_MEDIAL_PRIMITIVE_VALID
#define CDMM_NORROW 1e-8
#define CDMM_NewtonSolverPrecision 1e-15
	
	using namespace PolynoimalSolver;
	typedef ValidRange RootRange;

	qeal valueOfQuadircSurface2D(const qeal x, const qeal y, const qeal A, const qeal B, const qeal C, const qeal D, const qeal E, const qeal F);

	void genQuarticCoeffs(qeal & P_4, qeal & P_3, qeal & P_2, qeal & P_1, qeal & P_0, qeal A_2, qeal A_1, qeal A_0, qeal B_2, qeal B_1, qeal B_0, qeal C_2, qeal C_1, qeal C_0, qeal D_2, qeal D_1, qeal D_0, qeal S = 1.0);

	void genSexticCoeffs(qeal & P_6, qeal & P_5, qeal & P_4, qeal & P_3, qeal & P_2, qeal & P_1, qeal & P_0, qeal A_2, qeal A_1, qeal A_0, qeal B_2, qeal B_1, qeal B_0, qeal C_2, qeal C_1, qeal C_0, qeal S = 1.0);

	qeal SolveSexticEqForFTC(qeal & W_6, qeal & W_5, qeal & W_4, qeal & W_3, qeal & W_2, qeal & W_1, qeal & W_0, int banRootsNum, qeal* banRoots, RootRange & range, qeal norrow = CDMM_NORROW);

	qeal NewtonSolverSexticEqForFTC(qeal & W_6, qeal & W_5, qeal & W_4, qeal & W_3, qeal & W_2, qeal & W_1, qeal & W_0, Vector2 & range, qeal norrow = CDMM_NORROW);

	bool dcd(Vector3 c11, qeal r11, Vector3 c12, qeal r12, Vector3 c21, qeal r21, Vector3 c22, qeal r22, bool space_triangle = false, qeal norrow = CDMM_NORROW);
	bool dcd(Vector3 C1, Vector3  C2, Vector3 C3, qeal&  R1, qeal& R2, qeal& R3, bool space_triangle = false, qeal norrow = CDMM_NORROW);
	bool dcd(qeal A, qeal B, qeal C, qeal D, qeal E, qeal F, bool space_triangle = false, qeal norrow = CDMM_NORROW);

	bool ccd(Vector3& C1, Vector3& V1, Vector3& C2, Vector3& V2, Vector3& C3, Vector3& V3, qeal& R1, qeal& R2, qeal& R3, qeal& ftc, bool space_triangle, qeal norrow = CDMM_NORROW);

	bool ccd2(Vector3& C1, Vector3& V1, Vector3& C2, Vector3& V2, Vector3& C3, Vector3& V3, qeal& R1, qeal& R2, qeal& R3, qeal& ftc, bool space_triangle, qeal norrow = CDMM_NORROW);

	bool checkEndpointAlphaBetaCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, bool space_triangle = false, qeal norrow = CDMM_NORROW);

	bool checkAlphaIsZeroCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, qeal norrow = CDMM_NORROW);

	bool checkAlphaIsOneCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, qeal norrow = CDMM_NORROW);

	bool checkBetaIsZeroCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, qeal norrow = CDMM_NORROW);

	bool checkBetaIsOneCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, qeal norrow = CDMM_NORROW);

	bool checkAlphaPlusBetaIsOneCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, qeal norrow = CDMM_NORROW);

	bool checkAlphaBetaCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, bool space_triangle = false, qeal norrow = CDMM_NORROW);

	bool checkAlphaBetaCCD2(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, int collideType,  bool space_triangle = false, qeal norrow = CDMM_NORROW);

	bool checkAlphaBetaCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, int banRootsNum, qeal* banRoots, RootRange & pos, RootRange & neg, bool space_triangle = false, qeal norrow = CDMM_NORROW);

	bool checkAlphaBetaCCD2(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, int banRootsNum, qeal* banRoots, RootRange & pos, RootRange & neg, int collideType, bool space_triangle = false, qeal norrow = CDMM_NORROW);
}


#endif