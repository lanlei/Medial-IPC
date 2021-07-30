#include "IsotropicMaterial.h"
#include "StvkIsotropicMaterial.h"
#include "NeoHookeanIsotropicMaterial.h"

IsotropicMaterial::IsotropicMaterial() {}

IsotropicMaterial::~IsotropicMaterial() {}

std::shared_ptr<IsotropicMaterial> IsotropicMaterial::isotropicMaterialFactory(const std::string & materialType, BaseTetMeshHandle * tetHandle, int enableCompressionResistance, qeal compressionResistance)
{
	if (materialType.compare("NeoHookeanIsotropicMaterial") == 0)
	{
		return std::make_shared<NeoHookeanIsotropicMaterial>(tetHandle, enableCompressionResistance, compressionResistance);
	}
	else if (materialType.compare("StvkIsotropicMaterial") == 0)
	{
		return std::make_shared<StvkIsotropicMaterial>(tetHandle, enableCompressionResistance, compressionResistance);
	}
	else return std::make_shared<NeoHookeanIsotropicMaterial>(tetHandle, enableCompressionResistance, compressionResistance);
}
