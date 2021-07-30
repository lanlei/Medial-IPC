#pragma once
#ifndef BASE_SURFACE_HANDLE_H
#define BASE_SURFACE_HANDLE_H
#include "BaseSurfaceMesh.h"
#include "Commom\DataStructure\DataSet.h"

typedef DataSet<int> BaseSurfaceFaceSet;
typedef DataSet<int> BaseSurfacePointSet;

class BaseSurfaceHandle
{
public:
	BaseSurfaceHandle(): _surface(nullptr){}
	BaseSurfaceHandle(BaseSurfaceMesh* mesh);
	BaseSurfaceMesh* getMesh() const;
	virtual void init();
	virtual void computeFaceNormal();
	virtual void computePointNormal();
protected:
	BaseSurfaceMesh* _surface;
	std::vector<BaseSurfaceFaceSet> _faceSet;
	std::vector<BaseSurfacePointSet> _pointSet;
};
#endif