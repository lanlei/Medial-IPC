#pragma once
#ifndef BASE_MESH_BUFFER_POOL_H
#define BASE_MESH_BUFFER_POOL_H
/*
	Only Triangle Mesh
*/
#include "MatrixCore.h"
#include "Commom\BufferSerialization.h"
#include "BaseRenderMaterial.h"

class BaseSurfaceMeshBufferPool
{
public:
	BaseSurfaceMeshBufferPool() :
		meshNum(0),
		totalPointsNum(0),
		totalFacesNum(0),
		totalRenderGroupNum(0)
	{
	}

public:
	int meshNum;
	int totalPointsNum, totalFacesNum, totalRenderGroupNum;
	OverallVectorBuffer<qeal> pointsBuffer;
	OverallVectorBuffer<qeal> texCoordsBuffer;
	OverallVectorBuffer<qeal> pointsNormalBuffer;
	OverallVectorBuffer<qeal> facesNormalBuffer;
	OverallVectorBuffer<int> faceIndicesBuffer;
	OverallVectorBuffer<int> groupFacesNumBuffer;
	OverallVectorBuffer<int> pointFaceListBuffer;
};


#endif