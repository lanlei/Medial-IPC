#pragma once
#ifndef FINITE_ELEMENT_METHOD_MODEL_H
#define FINITE_ELEMENT_METHOD_MODEL_H
#include "MatrixCore.h"
#include "Model\BaseModel.h"
#include "Model/VolumetricMesh/IsotropicMaterial/IsotropicMaterial.h"
#include "Model/VolumetricMesh/IsotropicMaterial/NeoHookeanIsotropicMaterial.h"
#include "Model/VolumetricMesh/IsotropicMaterial/StvkIsotropicMaterial.h"
#include "Eigen/src/Geometry/AlignedBox.h"

namespace FiniteElementMethod
{
	void computeSVD(Matrix3 F, Matrix3& mU, Matrix3& mV, Matrix3& mSingularF, qeal tol, int flags);

	class FemModel : public BaseModel
	{
	public:
		FemModel() {}
	};

}


#endif