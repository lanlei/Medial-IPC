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

#ifndef  ISOTROPIC_MATERIAL_H_
#define ISOTROPIC_MATERIAL_H_

#include "MatrixCore.h"
#include "Model/VolumetricMesh/BaseTetMeshHandle.h"
#include "Model/VolumetricMesh/BaseTetMeshHandle.h"

class IsotropicMaterial
{
public:
	IsotropicMaterial();
	virtual~IsotropicMaterial();

	static std::shared_ptr<IsotropicMaterial> isotropicMaterialFactory(const std::string &materialType, BaseTetMeshHandle* tetHandle, int enableCompressionResistance = 0, qeal compressionResistance = 0.0);

	virtual void updateMaterial(BaseTetMeshHandle* tetHandle) {};
	virtual qeal computeEnergy(int eid, qeal* invariants) { return 0.0; };
	// invariants and gradient are 3-vectors
	virtual void computeEnergyGradient(int eid, qeal* invariants, qeal* gradinet) {};
	// invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).
	virtual void computeEnergyHessian(int eid, qeal* invariants, qeal* hessian) {};
	virtual void getLambdaLame(int eid, qeal& lambda) {}
	virtual void getMuLame(int eid, qeal& mu) {}
};



#endif 