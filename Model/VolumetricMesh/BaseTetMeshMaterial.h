#pragma once
#ifndef BASE_TET_MESH_MATERIAL_H
#define BASE_TET_MESH_MATERIAL_H
#include "MatrixCore.h"

class BaseTetMeshMaterial
{
public:		
	BaseTetMeshMaterial(const std::string name, qeal density) :_name(name), _density(density) {}
	BaseTetMeshMaterial(const BaseTetMeshMaterial& material) :_name(material._name), _density(material._density) {}
	BaseTetMeshMaterial(BaseTetMeshMaterial* material) :_name(material->_name), _density(material->_density) {}
	BaseTetMeshMaterial(): _density(0.0){}

	virtual ~BaseTetMeshMaterial(){}

	virtual std::string getName() const { return _name; }
	virtual qeal getDensity() const { return _density; }
	virtual void setName(const std::string name) { _name = name; }
	virtual void setDensity(qeal density) { _density = density; }

	// ENU = any isotropic material parameterized by E (Young's modulus), nu (Poisson's ratio)
	// ORTHOTROPIC = orthotropic anisotropic material
	// MOONEYRIVLIN = Mooney-Rivlin material
	typedef enum { INVALID, ENU } materialType;
	virtual materialType getType() { return type; }
protected:
	std::string _name;
	qeal _density;
	materialType type;
};


#endif