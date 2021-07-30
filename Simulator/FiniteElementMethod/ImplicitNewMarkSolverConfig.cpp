#include "ImplicitNewMarkSolverConfig.h"
#include <iomanip>
namespace FiniteElementMethod
{

	ImplicitNewMarkSolverConfig::ImplicitNewMarkSolverConfig(qeal timeStep, qeal massCoef, qeal stiffnessCoef, int maxIterNum, qeal ep, qeal markBeta, qeal markGamma) :
		newmarkBeta(markBeta),
		newmarkGamma(markGamma),
		dampingMassCoef(massCoef),
		dampingStiffnessCoef(stiffnessCoef),
		maxIterations(maxIterNum),
		epsilon(ep)
	{
		updateAlphas(timeStep);
	}

	void ImplicitNewMarkSolverConfig::updateAlphas(qeal timeStep)
	{
		alpha[0] = 1.0 / (newmarkBeta * timeStep * timeStep);
		alpha[1] = 1.0 / (newmarkBeta * timeStep);
		alpha[2] = (1.0 - 2.0 * newmarkBeta) / (2.0 * newmarkBeta);
		alpha[3] = newmarkGamma / (newmarkBeta * timeStep);
		alpha[4] = 1 - newmarkGamma / newmarkBeta;
		alpha[5] = (1.0 - newmarkGamma / (2.0 * newmarkBeta)) * timeStep;
	}

}

