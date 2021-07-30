#pragma once
#ifndef  STVK_ISOTROPIC_MATERIAL_H_
#define STVK_ISOTROPIC_MATERIAL_H_

#include "IsotropicMaterialWithCompressionResistance.h"

class StvkIsotropicMaterial : public IsotropicMaterialWithCompressionResistance
	{
	public:
		StvkIsotropicMaterial(BaseTetMeshHandle* tetHandle, int enableCompressionResistance = 0, qeal compressionResistance = 0.0);
		
		virtual ~StvkIsotropicMaterial();

		virtual void updateMaterial(BaseTetMeshHandle* tetHandle);

		virtual qeal computeEnergy(int eid, qeal * invariants);
		// invariants and gradient are 3-vectors
		virtual void computeEnergyGradient(int eid, qeal * invariants, qeal * gradient); 
		// invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).
		virtual void computeEnergyHessian(int eid, qeal * invariants, qeal * hessian); 

		virtual void getLambdaLame(int eleid, qeal &lambda);
		virtual void getMuLame(int eleid, qeal &mu);

	protected:
		qeal *_lambdaLame;
		qeal *_muLame;

		qeal _compressionResistance;
		qeal* _EdivNuFactor;
		virtual qeal getCompressionResistanceFactor(int eid);
	};

#endif