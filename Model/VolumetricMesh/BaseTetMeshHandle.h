#pragma once
#ifndef BASE_TET_MESH_HANDLE_H
#define BASE_TET_MESH_HANDLE_H
#include "BaseTetMesh.h"
#include "IsotropicMaterial/VolumetricMeshENuMaterial.h"
#include "Commom\DataStructure\DataSet.h"
#include "Commom\GeometryComputation.h"
#include <QOpenGLFunctions>

typedef DataSet<int> BaseTetElementSet;
typedef DataSet<int> BaseTetNodeSet;
class BaseTetElementRegion;
class BaseTetElementParam;
class BaseTetNodeParam;

class BaseTetMeshHandle
{
public:
	BaseTetMeshHandle(): _tetMesh(nullptr), renderElementSetId(0){}
	BaseTetMeshHandle(BaseTetMesh* tetMesh);
	BaseTetMesh* getTetMesh() const;

	virtual void init();
	virtual void correctElementNodeOrder();
	virtual void computeVolume();
	virtual void computeElementFaceArea();
	virtual void computeTetMeshNormal();
	void initElementABCD();
	virtual void computeElementShapeFunctionDerivate();
	virtual void markBoundaryElement();
	//
	BaseTetElementParam* getTetElementParam(const int eid);
	BaseTetNodeParam* getTetNodeParam(const int nid);

	BaseTetMeshMaterial* getElementMaterialById(int id) const;
	BaseTetMeshMaterial* getElementMaterial(int eid) const;
	BaseTetMeshMaterial* getNodeMaterial(int nid) const;

	bool readTetMeshSetInfo(std::string& filename);
	void writeTetMeshSetInfo(std::string& filename);

	void writeTetMeshSetInfoDebug(std::string& filename);

	void setTetMeshSet(std::vector<BaseTetElementSet>& eleSet) {
		_elementSet.clear();
		_elementSet = eleSet;
	}

	std::vector<BaseTetElementSet>& getTetMeshElementSet() {return _elementSet;}
	std::vector<BaseTetNodeSet>& getTetMeshNodeSet() { return _nodeSet; }

	virtual void renderTetMesh(QOpenGLFunctions* f);
	void setRenderElementSetId(int id) { renderElementSetId = id; }
//protected:
public:
	BaseTetMesh* _tetMesh;
	std::vector<bool> _boundaryElementFlag;
	std::vector<std::vector<int>> _elementNeighbor;
	std::vector<BaseTetElementSet> _elementSet;
	std::vector<BaseTetNodeSet> _nodeSet;
	std::vector<BaseTetElementRegion> _elementRegion;
	std::vector<BaseTetElementParam*> _elementsParam;
	std::vector<BaseTetNodeParam*> _nodesParam;

	std::vector<BaseTetMeshMaterial*> _materials;
	std::vector<BaseTetElementRegion*> _regions;
	std::vector<int> _elementMaterials;
	std::vector<int> _nodeMaterials;

	//reder 
	int renderElementSetId;

};

class BaseTetElementRegion
{
public:
	BaseTetElementRegion(int materialIndex, int setIndex)
	{
		_materialIndex = materialIndex;
		_setIndex = setIndex;
	}

	inline int getMaterialIndex() const {
		return _materialIndex;
	}
	inline int getSetIndex() const {
		return _setIndex;
	}

	virtual bool readBinary(std::ifstream fin) {}

	virtual bool writeBinary(std::ifstream fin) {}

protected:
	int _setIndex, _materialIndex;
};

class BaseTetElementParam
{
public:
	BaseTetElementParam(int eleid = -1);
	virtual void computeElementDeformationByShapeMatrix(Matrix3 &Ds, qeal* displacement, int* indices);
	virtual void computeShapeFunctionDerivate(qeal* points);
	virtual void computeVolume(qeal* points);
	virtual void computeCenter(qeal* points);
	virtual void computeFaceArea(qeal* points);
	virtual void computeElementNormal(qeal* points);
	virtual void computeABCD();

	virtual bool readBinary(std::ifstream fin) { return true; }

	virtual bool writeBinary(std::ifstream fin) { return true; }

	int eleid;

	Matrix4 shapeFunction;
	Matrix4 invShapeFunction;
	Matrix3 Dm;
	Matrix3 invDm;
	Matrix3 invDmT;
	Vector3 phiGradient[4];
	MatrixX phiDerivate;

	MatrixX dFdu;

	int faceIndices[12];
	qeal area[4];
	Vector3 internalForce[4];
	Vector3 nodesNormal[4];
	Vector3 facesNormal[4];

	Matrix3 A[4][4];
	qeal B[4][4];
	Vector3 C[4][4][4];
	qeal D[4][4][4][4];

	Vector3 center;
	qeal volume;
	int region;
};

class BaseTetNodeParam
{
public:
	BaseTetNodeParam() {};

	virtual bool readBinary(std::ifstream fin) { return true; }

	virtual bool writeBinary(std::ifstream fin) { return true; }

	qeal volume;
	int setIndex;
};

#endif