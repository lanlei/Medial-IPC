#pragma once
#ifndef BASE_MEDIAL_MESH_BUFFER_POOL_H
#define BASE_MEDIAL_MESH_BUFFER_POOL_H
#include "MatrixCore.h"
#include "Commom\BufferSerialization.h"

class BaseMedialMeshBufferPool
{
public:
	BaseMedialMeshBufferPool() :
		totalMedialPoinsNum(0),
		totalMedialPrimitivesNum(0),
		totalMedialConesNum(0),
		totalMedialSlabsNum(0)
	{
	}

public:
	int totalMedialPoinsNum, totalMedialPrimitivesNum, totalMedialConesNum, totalMedialSlabsNum;
	OverallVectorBuffer<qeal> medialPointsBuffer;
	OverallVectorBuffer<qeal> medialRadiusBuffer;
	OverallVectorBuffer<int> medialPrimitiveIndicesBuffer;
	OverallVectorBuffer<int> medialPointsNeighborListBuffer;
	OverallVectorBuffer<int> medialPointsConesBuffer;
	OverallVectorBuffer<int> medialPointsSlabsBuffer;
};


#endif