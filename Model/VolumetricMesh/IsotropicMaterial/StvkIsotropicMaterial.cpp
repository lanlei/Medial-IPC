#include "StvkIsotropicMaterial.h"
#include "VolumetricMeshENuMaterial.h"
#include <math.h>

StvkIsotropicMaterial::StvkIsotropicMaterial(BaseTetMeshHandle* tetHandle, int enableCompressionResistance, qeal compressionResistance) : IsotropicMaterialWithCompressionResistance(enableCompressionResistance), _compressionResistance(compressionResistance)
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
	
	StvkIsotropicMaterial::~StvkIsotropicMaterial()
	{
		free(_EdivNuFactor);
		free(_muLame);
		free(_lambdaLame);
	}

	void StvkIsotropicMaterial::updateMaterial(BaseTetMeshHandle* tetHandle)
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

	qeal StvkIsotropicMaterial::computeEnergy(int eid, qeal * invariants)
	{
		qeal IC = invariants[0];
		qeal IIC = invariants[1];
		//qeal IIIC = invariants[2]; // not needed for StVK
		qeal energy = 0.125 * _lambdaLame[eid] * (IC - 3.0) * (IC - 3.0) + 0.25 * _muLame[eid] * (IIC - 2.0 * IC + 3.0);

		addCompressionResistanceEnergy(eid, invariants, &energy);

		return energy;
	}

	void StvkIsotropicMaterial::computeEnergyGradient(int eid, qeal * invariants, qeal * gradient)
	{
		qeal IC = invariants[0];
		gradient[0] = 0.25 * _lambdaLame[eid] * (IC - 3.0) - 0.5 * _muLame[eid];
		gradient[1] = 0.25 * _muLame[eid];
		gradient[2] = 0.0;

		addCompressionResistanceGradient(eid, invariants, gradient);
	}

	void StvkIsotropicMaterial::computeEnergyHessian(int eid, qeal * invariants, qeal * hessian)
	{
		// 11
		hessian[0] = 0.25 * _lambdaLame[eid];
		// 12
		hessian[1] = 0.0;
		// 13
		hessian[2] = 0.0;
		// 22
		hessian[3] = 0.0;
		// 23
		hessian[4] = 0.0;
		// 33
		hessian[5] = 0.0;

		addCompressionResistanceHessian(eid, invariants, hessian);
	}

	void StvkIsotropicMaterial::getLambdaLame(int eleid, qeal & lambda)
	{
		lambda = _lambdaLame[eleid];
	}

	void StvkIsotropicMaterial::getMuLame(int eleid, qeal & mu)
	{
		mu = _muLame[eleid];
	}

	qeal StvkIsotropicMaterial::getCompressionResistanceFactor(int eid)
	{
		return _EdivNuFactor[eid];
	}