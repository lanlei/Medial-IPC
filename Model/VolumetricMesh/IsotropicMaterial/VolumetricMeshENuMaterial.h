#pragma once
/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 3.1                               *
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2016 USC *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 *                                                                       *
 * Research: Jernej Barbic, Fun Shing Sin, Daniel Schroeder,             *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC                 *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#ifndef VOLUMETRIC_MESH_ENU_MATERIAL_H_
#define VOLUMETRIC_MESH_ENU_MATERIAL_H_

#include "Model\VolumetricMesh\BaseTetMeshMaterial.h"

 // stores an isotropic material specified by E (Young's modulus), nu (Poisson's ratio), and density
	 // such a material specification is very common: (corotational) linear FEM, StVK, etc.

class ENuMaterial : public BaseTetMeshMaterial
{
public:
	ENuMaterial(std::string name, qeal density, qeal E, qeal nu);
	ENuMaterial(const ENuMaterial & eNuMaterial);
	virtual ~ENuMaterial() {}
	virtual BaseTetMeshMaterial * clone() const;
	virtual BaseTetMeshMaterial::materialType getType();

	qeal getE() const; // Young's modulus
	qeal getNu() const; // Poisson's ratio
	qeal getLambda() const; // Lame's lambda coefficient
	qeal getMu() const; // Lame's mu coefficient
	void setE(qeal E);
	void setNu(qeal nu);
protected:
	qeal _E;
	qeal _nu;
};

// obtain pointer to ENuMaterial (necessary inside classes that assume ENu material)
ENuMaterial * downcastENuMaterial(BaseTetMeshMaterial * material); // performs a check via getType and returns NULL if material is not ENU



#endif