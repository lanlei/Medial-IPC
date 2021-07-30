#include "FemModel.h"
//#include "IsotropicMaterial\StvkIsotropicMaterial.h"
//#include "IsotropicMaterial\NeoHookeanIsotropicMaterial.h"
#include "Commom\SparseMatrixTopology.h"

namespace FiniteElementMethod
{
	void computeSVD(Matrix3 F, Matrix3 & mU, Matrix3 & mV, Matrix3 & mSingularF, qeal tol, int flags)
	{
		Eigen::JacobiSVD<Matrix3> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		mU = svd.matrixU();
		mV = svd.matrixV();
		Vector3 sv = svd.singularValues();
		mSingularF.setZero();
		mSingularF.data()[0] = sv.data()[0];
		mSingularF.data()[4] = sv.data()[1];
		mSingularF.data()[8] = sv.data()[2];
	}
}
