#pragma once
#ifndef BASE_MEDIAL_MESH_Handle_H
#define BASE_MEDIAL_MESH_Handle_H
#include "BaseMedialMesh.h"
#include "Commom\DataStructure\DataSet.h"
#include <QOpenGLFunctions>

typedef DataSet<int> BaseMedialPointSet;
typedef DataSet<int> BaseMedialPrimitiveSet;

class BaseMedialMeshHandle
{
public:
	BaseMedialMeshHandle():_medialMesh(nullptr), renderMedialMeshSetId(0){}
	BaseMedialMeshHandle(BaseMedialMesh* mesh);
	BaseMedialMesh* getMesh() const;
	virtual void init();

	bool readMedialMeshSetInfo(std::string& filename);
	void writeMedialMeshSetInfo(std::string& filename);
	void setMedialMeshSet(std::vector<BaseMedialPrimitiveSet>& primitiveSet, std::vector<BaseMedialPointSet>& pointSet)
	{
		_primitiveSet = primitiveSet;
		_pointSet = pointSet;

	}

	void getMedialMeshSet(std::vector<BaseMedialPrimitiveSet>& primitiveSet, std::vector<BaseMedialPointSet>& pointSet)
	{
		primitiveSet = _primitiveSet;
		pointSet = _pointSet;
	}

	virtual void renderMedialMesh(QOpenGLFunctions* f);
	void setRenderMedialMeshSetId(int id) { renderMedialMeshSetId = id; }
protected:
	BaseMedialMesh* _medialMesh;
	std::vector<BaseMedialPrimitiveSet> _primitiveSet;
	std::vector<BaseMedialPointSet> _pointSet;

	int renderMedialMeshSetId;
};

#endif