#include "BaseSurfaceHandle.h"

BaseSurfaceHandle::BaseSurfaceHandle(BaseSurfaceMesh * mesh): _surface(mesh)
{
}

BaseSurfaceMesh * BaseSurfaceHandle::getMesh() const
{
	return _surface;
}

void BaseSurfaceHandle::init()
{
}

void BaseSurfaceHandle::computeFaceNormal()
{
}

void BaseSurfaceHandle::computePointNormal()
{
}
