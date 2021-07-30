#pragma once
#ifndef BASE_MESH_H
#define BASE_MESH_H
/*
	Only Triangle Mesh
*/
#include <memory>
#include <vector>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include "MatrixCore.h"
#include "BaseRenderMaterial.h"
#include "Commom\BaseFrame.h"
#include "Commom\BufferSerialization.h"
#include "BaseSurfaceMeshBufferPool.h"

class BaseSurfaceMeshBuffer
{
	typedef FragmentVectorBufferPtr <qeal> Points;
	typedef FragmentVectorBufferPtr <qeal> Texcoords;
	typedef FragmentVectorBufferPtr <qeal> Normals;
	typedef FragmentVectorBufferPtr <int> Indices;
public:
	BaseSurfaceMeshBuffer() :
		pointsNum(0),
		facesNum(0),
		renderGroupNum(0)
	{}
	virtual void refresh(BaseSurfaceMeshBufferPool* pool)
	{
		points.buffer = pool->pointsBuffer.buffer.data() + points.offset;
		texCoords.buffer = pool->texCoordsBuffer.buffer.data() + texCoords.offset;
		pointNormals.buffer = pool->pointsNormalBuffer.buffer.data() + pointNormals.offset;
		faceNormals.buffer = pool->facesNormalBuffer.buffer.data() + faceNormals.offset;
		faceIndices.buffer = pool->faceIndicesBuffer.buffer.data() + faceIndices.offset;

		for (int i = 0; i < pointFaceIndices.size(); i++)
			pointFaceIndices[i].buffer = pool->pointFaceListBuffer.buffer.data() + pointFaceIndices[i].offset;
		renderGroupFaceNum.buffer = pool->groupFacesNumBuffer.buffer.data() + renderGroupFaceNum.offset;

		pointsOffset = points.offset / 3;
		faceOffset = faceIndices.offset / 3;
	}

	int pointsNum, facesNum, renderGroupNum;
	Points points;
	Texcoords texCoords;
	Normals pointNormals;
	Normals faceNormals;
	Indices faceIndices;
	std::vector <Indices> pointFaceIndices;
	std::vector<BaseRenderMaterial*> renderMaterials;
	Indices renderGroupFaceNum;

	void copySurfacePointsToBuffer(qeal* buffer, int size = 0);
	void copyBufferToSurfacePoints(qeal* buffer, int size = 0);
	
	Vector3 getSurfacePoint(int nid);
	void setSurfacePoint(int nid, qeal* p);

	Vector3 getSurfacePointNormal(int nid);
	void setSurfacePointNormal(int nid, qeal* n);

	Vector3 getSurfaceFaceNormal(int nid);
	void setSurfaceFaceNormal(int nid, qeal* n);

	Vector3i getSurfaceFace(int fid);

	int getSurfacePointOverallId(int nid);
	int getSurfaceFaceOverallId(int fid);
protected:
	int pointsOffset;
	int faceOffset;
};

class BaseSurfaceMesh : public BaseSurfaceMeshBuffer, public BaseFrame
{
	enum MeshType
	{
		DYNAMIC_MESH = 0,
		STATICS_MESH = 1
	};
public:
	BaseSurfaceMesh(MeshType type = DYNAMIC_MESH)
	{
		hasTextureCoords = false;
		hasNormals = false;
		shadow = true;
		_vertexArrayBuf = (std::unique_ptr<QOpenGLBuffer>) new QOpenGLBuffer();
		if (!_vertexArrayBuf->isCreated())
			_vertexArrayBuf->create();
		_hide = false;
		_type = type;
	}

	virtual void setNameAndDir(const std::string filename);

	bool readMeshFromObjFormat(const std::string filename, BaseSurfaceMeshBufferPool* pool);

	bool writeMeshToObjFormat(const std::string filename);

	virtual void render(QOpenGLShaderProgram* program, QOpenGLFunctions* f, bool drawEdge = false);
	virtual void updateVBO(); // only running in render();

	virtual qeal uniform();
	virtual void computeBBox();
	virtual void getCenter(qeal& cx, qeal& cy, qeal& cz);

	//render status
	void enableHide(bool enable) { _hide = enable; }
	bool isHide() { return _hide; }
	MeshType getType() { return _type; }
	void setType(MeshType type) { _type = type;}
	bool isDynamic() { return _type == DYNAMIC_MESH ? true : false; }

	//
	void translateMesh(const qeal x, const qeal y, const qeal z);
	void scaleMesh(qeal s, const qeal cx, const qeal cy, const qeal cz);
	void rotateMesh(const qglviewer::Quaternion oldR, const qglviewer::Quaternion R, const qeal cx, const qeal cy, const qeal cz);
protected:
	virtual void initVBO();
public:
	std::string name, dir, format, nickName;
	std::array<qeal, 6> bbox;

	bool hasNormals;
	bool hasTextureCoords;
	bool shadow;

protected:
	// shader buffer
	std::unique_ptr<QOpenGLBuffer> _vertexArrayBuf;
	std::vector<qeal> _renderVertexBuffer;
	bool _hide;
	MeshType _type;
};

#endif