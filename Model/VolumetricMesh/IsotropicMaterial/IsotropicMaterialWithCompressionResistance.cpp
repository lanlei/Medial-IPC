#include "IsotropicMaterialWithCompressionResistance.h"
#include <math.h>

IsotropicMaterialWithCompressionResistance::IsotropicMaterialWithCompressionResistance(int enableCompressionResistance) : IsotropicMaterial(), _enableCompressionResistance(enableCompressionResistance)
{
}

IsotropicMaterialWithCompressionResistance::~IsotropicMaterialWithCompressionResistance() {}


qeal IsotropicMaterialWithCompressionResistance::getCompressionResistanceFactor(int eid)
{
	return 1.0; // generic
}

void IsotropicMaterialWithCompressionResistance::addCompressionResistanceEnergy(int eid, qeal* invariants, qeal * energy)
{
	if (_enableCompressionResistance)
	{
		qeal IIIC = invariants[2];
		qeal J = sqrt(IIIC);

		if (J < 1)
		{
			qeal compressionResistanceFactor = getCompressionResistanceFactor(eid);
			*energy += -compressionResistanceFactor * (J - 1.0) * (J - 1.0) * (J - 1.0) / 2592.0;
		}
	}
}

void IsotropicMaterialWithCompressionResistance::addCompressionResistanceGradient(int eid, qeal * invariants, qeal * gradient)
{
	if (_enableCompressionResistance)
	{
		double IIIC = invariants[2];
		double J = sqrt(IIIC);

		if (J < 1)
		{
			double compressionResistanceFactor = getCompressionResistanceFactor(eid);
			gradient[2] += -compressionResistanceFactor * (J - 1.0) * (J - 1.0) / (1728.0 * J);
		}
	}
}

void IsotropicMaterialWithCompressionResistance::addCompressionResistanceHessian(int eid, qeal * invariants, qeal * hessian)
{
	if (_enableCompressionResistance)
	{
		qeal IIIC = invariants[2];
		qeal J = sqrt(IIIC);

		if (J < 1.0)
		{
			qeal compressionResistanceFactor = getCompressionResistanceFactor(eid);
			hessian[5] += compressionResistanceFactor * (1.0 - J) * (1.0 + J) / (3456.0 * J * J * J);
		}
	}
}