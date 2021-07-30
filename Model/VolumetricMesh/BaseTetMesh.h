#pragma once
#ifndef BASE_TET_MESH_H
#define BASE_TET_MESH_H

#include <memory>
#include <vector>
#include "Commom\BufferSerialization.h"
#include "BaseTetMeshBufferPool.h"
#include <QGLViewer/vec.h>
#include <QGLViewer/Quaternion.h>
#include "MatrixCore.h"

class BaseTetMeshBuffer
{
	typedef FragmentVectorBufferPtr <qeal> Points;
	typedef FragmentVectorBufferPtr <int> Indices;
public:
	BaseTetMeshBuffer() :
		tetPointsNum(0),
		tetElementNum(0)
	{}
	virtual void refresh(BaseTetMeshBufferPool* pool)
	{
		tetPoints.buffer = pool->tetPointsBuffer.buffer.data() + tetPoints.offset;
		tetElementIndices.buffer = pool->tetElementIndicesBuffer.buffer.data() + tetElementIndices.offset;

		for (int i = 0; i < tetPointElementList.size(); i++)
			tetPointElementList[i].buffer = pool->tetPointElementListBuffer.buffer.data() + tetPointElementList[i].offset;
		for (int i = 0; i < tetPointNeighborList.size(); i++)
			tetPointNeighborList[i].buffer = pool->tetPointNeighborListBuffer.buffer.data() + tetPointNeighborList[i].offset;

		tetPointsIdOffset = tetPoints.offset / 3;
		tetElementsIdOffset = tetElementIndices.offset / 4;
	}

	int tetPointsNum, tetElementNum;
	Points tetPoints;
	Indices tetElementIndices;
	std::vector <Indices> tetPointElementList;
	std::vector <Indices> tetPointNeighborList;

	void copyTetPointsToBuffer(qeal* buffer, int size = 0);
	void copyBufferToTetPoints(qeal* buffer, int size = 0);

	Vector3 getTetPoint(int nid);
	void setTetPoint(int nid, qeal* p);
	Vector4i getTetElement(int eleid);
	void setTetElement(int eleid, int* indices);

	VectorX getTetElementPoint(int eid);

	Vector4i getTetElementOverall(int eid);
	int getTetElementOverallId(int eid);	
	int getTetPointOverallId(int nid);
protected:
	int tetPointsIdOffset;
	int tetElementsIdOffset;
};

class BaseTetMesh : public BaseTetMeshBuffer
{
public:
	BaseTetMesh() : _hide(true), _valid(false), _hasNodes(false), _hasElements(false), _transparent(false)
	{
	}
	bool readMeshFromTetFormat(const std::string filename, BaseTetMeshBufferPool* pool);
	bool readNodeFromTetFormat(const std::string filename, BaseTetMeshBufferPool* pool);
	bool readElementFromTetFormat(const std::string filename, BaseTetMeshBufferPool* pool);

	bool writeMeshToTetFormat(const std::string filename);
	bool writeNodeToTetFormat(const std::string filename);
	bool writeElementToTetFormat(const std::string filename);

	bool isInsideTetElement(const int eid, const qeal* point, qeal* weight);
	void computeTetElementBarycentricWeights(const int eid, const qeal* point, qeal* weight);
	int searchCloseTetNode(const qeal* point);

	//
	virtual void uniform(const qeal div);
	void translateMesh(const qeal x, const qeal y, const qeal z);
	void scaleMesh(qeal s, const qeal cx, const qeal cy, const qeal cz);
	void rotateMesh(const qglviewer::Quaternion oldR, const qglviewer::Quaternion R, const qeal cx, const qeal cy, const qeal cz);

	void enableHide(bool enable) { _hide = enable; }
	bool isHide() { return _hide; }
	void enableTransparent(bool enable) { _transparent = enable; }
	bool isTransparent() { return _transparent; }
	void setTetMeshValid(bool valid) {_valid = valid;}
	bool isTetMeshValid() { return _valid; }

	friend class BaseTetMeshHandle;
protected:
	bool _hide;
	bool _valid;
	bool _hasNodes;
	bool _hasElements;
	bool _transparent;
};


#endif