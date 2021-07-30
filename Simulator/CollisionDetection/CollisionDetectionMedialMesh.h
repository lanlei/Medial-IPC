#pragma once
#ifndef COLLISION_DETECTION_ON_MEDIAL_MESH_H
#define COLLISION_DETECTION_ON_MEDIAL_MESH_H
#include "DataCore.h"
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
	bool dcd(Vector3 c11, qeal r11, Vector3 c12, qeal r12, Vector3 c21, qeal r21, Vector3 c22, qeal r22, bool space_triangle = false, qeal norrow = CDMM_NORROW);
	bool dcd(Vector3 C1, Vector3  C2, Vector3 C3, qeal&  R1, qeal& R2, qeal& R3, bool space_triangle = false, qeal norrow = CDMM_NORROW);
	bool dcd(qeal A, qeal B, qeal C, qeal D, qeal E, qeal F, bool space_triangle = false, qeal norrow = CDMM_NORROW);

}


#endif