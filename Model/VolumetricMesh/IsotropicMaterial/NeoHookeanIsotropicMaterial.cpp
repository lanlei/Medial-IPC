#include "NeoHookeanIsotropicMaterial.h"
#include "VolumetricMeshENuMaterial.h"
#include <math.h>

NeoHookeanIsotropicMaterial::NeoHookeanIsotropicMaterial(BaseTetMeshHandle* tetHandle, int enableCompressionResistance, qeal compressionResistance)
	{
		int elementNum = tetHandle->getTetMesh()->tetElementNum;
		_lambdaLame = (qeal*)malloc(sizeof(qeal) * elementNum);
		_muLame = (qeal*)malloc(sizeof(qeal*) * elementNum);

		if (enableCompressionResistance)
			_EdivNuFactor = (qeal*)malloc(sizeof(qeal) * elementNum);
		else
			_EdivNuFactor = NULL;

		for (int eid = 0; eid < elementNum; eid++)
		{
			BaseTetMeshMaterial* material = tetHandle->getElementMaterial(eid);
			ENuMaterial * eNuMaterial = downcastENuMaterial(material);
			if (eNuMaterial == NULL)
			{
				printf("Error: NeoHookeanIsotropicMaterial: mesh does not consist of E, nu materials.\n");
				throw 1;
			}

			_lambdaLame[eid] = eNuMaterial->getLambda();
			_muLame[eid] = eNuMaterial->getMu();

			if (enableCompressionResistance)
				_EdivNuFactor[eid] = compressionResistance * eNuMaterial->getE() / (1.0 - 2.0 * eNuMaterial->getNu());
		}
	}

	NeoHookeanIsotropicMaterial::~NeoHookeanIsotropicMaterial()
	{
		free(_EdivNuFactor);
		free(_muLame);
		free(_lambdaLame);
	}

	void NeoHookeanIsotropicMaterial::updateMaterial(BaseTetMeshHandle* tetHandle)
	{
		int elementNum = tetHandle->getTetMesh()->tetElementNum;
		for (int eid = 0; eid < elementNum; eid++)
		{
			BaseTetMeshMaterial* material = tetHandle->getElementMaterial(eid);
			ENuMaterial * eNuMaterial = downcastENuMaterial(material);
			if (eNuMaterial == NULL)
			{
				printf("Error: NeoHookeanIsotropicMaterial: mesh does not consist of E, nu materials.\n");
				throw 1;
			}

			_lambdaLame[eid] = eNuMaterial->getLambda();
			_muLame[eid] = eNuMaterial->getMu();

			if (_enableCompressionResistance)
			{
				_EdivNuFactor[eid] = _compressionResistance * eNuMaterial->getE() / (1.0 - 2.0 * eNuMaterial->getNu());
			}
		}
	}

	qeal NeoHookeanIsotropicMaterial::computeEnergy(int eid, qeal * invariants)
	{
		qeal IC = invariants[0];
		qeal IIIC = invariants[2];
		qeal J = sqrt(IIIC);
		qeal logJ = log(J);

		// Note: computation of J and logJ will fail for an inverted element.
		// The IsotropicHyperelasticFEM class will prevent inversions (assuming proper
		// threshold was set), so normally this is not an issue.

		qeal energy = 0.5 * _muLame[eid] * (IC - 3.0) - _muLame[eid] * logJ + 0.5 * _lambdaLame[eid] * logJ * logJ;

		addCompressionResistanceEnergy(eid, invariants, &energy);

		return energy;
	}

	void NeoHookeanIsotropicMaterial::computeEnergyGradient(int eid, qeal * invariants, qeal * gradient)
	{
		qeal IIIC = invariants[2];
		gradient[0] = 0.5 * _muLame[eid];
		gradient[1] = 0.0;
		gradient[2] = (-0.5 * _muLame[eid] + 0.25 * _lambdaLame[eid] * log(IIIC)) / IIIC;

		addCompressionResistanceGradient(eid, invariants, gradient);
	}

	void NeoHookeanIsotropicMaterial::computeEnergyHessian(int eid, qeal * invariants, qeal * hessian)
	{
		qeal IIIC = invariants[2];
		// 11
		hessian[0] = 0.0;
		// 12
		hessian[1] = 0.0;
		// 13
		hessian[2] = 0.0;
		// 22
		hessian[3] = 0.0;
		// 23
		hessian[4] = 0.0;
		// 33
		hessian[5] = (0.25 * _lambdaLame[eid] + 0.5 * _muLame[eid] - 0.25 * _lambdaLame[eid] * log(IIIC)) / (IIIC * IIIC);

		addCompressionResistanceHessian(eid, invariants, hessian);
	}

	void NeoHookeanIsotropicMaterial::getLambdaLame(int eleid, qeal & lambda)
	{
		lambda = _lambdaLame[eleid];
	}

	void NeoHookeanIsotropicMaterial::getMuLame(int eleid, qeal & mu)
	{
		mu = _muLame[eleid];
	}
	
	qeal NeoHookeanIsotropicMaterial::getCompressionResistanceFactor(int eid)
	{
		return _EdivNuFactor[eid];
	}
