#include "VolumetricMeshENuMaterial.h"

ENuMaterial::ENuMaterial(std::string name, qeal density, qeal E, qeal nu) :BaseTetMeshMaterial(name, density), _E(E), _nu(nu) {}

ENuMaterial::ENuMaterial(const ENuMaterial & eNuMaterial)
	: BaseTetMeshMaterial(eNuMaterial.getName(), eNuMaterial.getDensity()), _E(eNuMaterial.getE()), _nu(eNuMaterial.getNu()) {}

BaseTetMeshMaterial * ENuMaterial::clone() const
{
	return new ENuMaterial(*this);
}

BaseTetMeshMaterial::materialType ENuMaterial::getType()
{
	return BaseTetMeshMaterial::ENU;
}

qeal ENuMaterial::getE() const
{
	return _E;
}

qeal ENuMaterial::getNu() const
{
	return _nu;
}

qeal ENuMaterial::getLambda() const
{
	return (_nu * _E) / ((1 + _nu) * (1 - 2 * _nu));
}

qeal ENuMaterial::getMu() const
{
	return _E / (2 * (1 + _nu));
}

void ENuMaterial::setE(qeal E)
{
	_E = E;
}

void ENuMaterial::setNu(qeal nu)
{
	_nu = nu;
}


ENuMaterial * downcastENuMaterial(BaseTetMeshMaterial * material)
{
	if (material->getType() != BaseTetMeshMaterial::ENU)
		return NULL;

	return (ENuMaterial*)material;
}