#pragma once
#ifndef BASE_MODEL_H
#define BASE_MODEL_H
#include "TriangleMesh/BaseSurfaceHandle.h"
#include "VolumetricMesh/BaseTetMeshHandle.h"
#include "MedialMesh/BaseMedialMedialHandle.h"
#include "Commom\FileIO.h"
#include <unordered_set>

class BaseModel : public BaseSurfaceMesh, public BaseTetMesh, public BaseMedialMesh
{
public:
	BaseModel();
	~BaseModel();

	virtual void refreshBuffer(BaseSurfaceMeshBufferPool* sPool, BaseTetMeshBufferPool* tPool, BaseMedialMeshBufferPool* mPool);

	virtual void init();
	virtual void initMeshesHandel();
	virtual qeal uniform();

	virtual void conbineTetAndMat(BaseTetMeshBufferPool * tetPool, BaseMedialMeshBufferPool* matPool, const std::string newTetNodeFilename, const std::string newTetElementFilename, const std::string newMedialMeshFilename);

	virtual void alignAllMesh(qeal* tetPointsPtr);
	virtual void computeTetMeshInterpolation();
	virtual void updateSurfaceFromTetMesh();
	virtual void updateSurfaceFromTetMesh(qeal* currentTetPoints);
	virtual void updateMedialMeshFromTetMesh(qeal* currentTetPoints);
	//
	virtual void translateModel(const qeal x, const qeal y, const qeal z);
	virtual void scaleModel(qeal s);
	virtual void rotateModel(const qeal rx, const qeal ry, const qeal rz, const qeal sita, bool oldClear = false);

	void getTranslation(qeal& x, qeal& y, qeal& z) { x = _tx; y = _ty; z = _tz; }
	void getRotation(qeal& x, qeal& y, qeal& z, qeal& sita) {
		x = _rx; y = _ry; z = _rz; sita = _rsita;
	}
	void getScale(qeal& scale) { scale =  _scale; }

	void setValid(bool flag) { _valid = flag; }
	bool isValid() { return _valid; }

	void enableGravityForce(bool flag) { enableGravity = flag; }

	BaseSurfaceHandle* getSurfaceHandle() {
		return _smHandle;
	}

	BaseTetMeshHandle* getTetMeshHandle() {
		return _tmHandle;
	}

	BaseMedialMeshHandle* getMedialMeshHandle() {
		return _mmHandle;
	}

	//
	// UI
	virtual void clickModelList(){}
	virtual void clickTetElementSetList(int id) {}
	virtual void clickMedialMeshSetList(int id) {}
	//


	std::vector<int> tetMeshInterpolationIndices;
	std::vector<qeal> tetMeshInterpolationWeight;
	bool enableGravity;
protected:
	BaseSurfaceHandle* _smHandle; // surface mesh
	BaseTetMeshHandle* _tmHandle; // tet mesh
	BaseMedialMeshHandle* _mmHandle; // medial mesh


	qeal _tx, _ty, _tz;
	qeal _scale;
	qeal _rx, _ry, _rz, _rsita;

	bool _valid;

};
#endif