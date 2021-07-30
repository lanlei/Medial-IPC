#pragma once
#ifndef BASE_TET_MESH_BUFFER_POOL_H
#define BASE_TET_MESH_BUFFER_POOL_H
/*
	Only Tetgen File
*/
#include "MatrixCore.h"
#include "Commom\BufferSerialization.h"

class BaseTetMeshBufferPool
{
public:
	BaseTetMeshBufferPool() :
		totalTetPointsNum(0),
		totalTetElementNum(0)
	{
	}

public:
	int totalTetPointsNum, totalTetElementNum;
	OverallVectorBuffer<qeal> tetPointsBuffer;
	OverallVectorBuffer<int> tetElementIndicesBuffer;
	OverallVectorBuffer<int> tetPointElementListBuffer;
	OverallVectorBuffer<int> tetPointNeighborListBuffer;
};


#endif