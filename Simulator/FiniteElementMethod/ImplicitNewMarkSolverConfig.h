#pragma once
#ifndef IMPLICIT_NEW_MARK_SOLVER_CONFIG_H
#define IMPLICIT_NEW_MARK_SOLVER_CONFIG_H 
#include "MatrixCore.h"

namespace FiniteElementMethod
{
	struct ImplicitNewMarkSolverConfig
	{
		ImplicitNewMarkSolverConfig(qeal timeStep, qeal massCoef = 0.01, qeal stiffnessCoef = 0.01, int maxIterations = 10, qeal epsilon = 1e-6, qeal newMarkBeta = 0.25, qeal newMarkGamma = 0.5);

		qeal newmarkBeta, newmarkGamma;
		qeal dampingMassCoef, dampingStiffnessCoef;
		qeal alpha[6];
		qeal epsilon;
		int maxIterations;

		void updateAlphas(qeal timeStep);
	};
}



#endif
