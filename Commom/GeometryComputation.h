#pragma once
#ifndef GEOMETRY_COMPUTATION_H
#define GEOMETRY_COMPUTATION_H
#include "MatrixCore.h"

qeal TriangleArea(const Vector3& v0, const Vector3& v1, const Vector3& v2);

void TriangleNormal(const Vector3& v0, const Vector3& v1, const Vector3& v2, Vector3& n);

bool DistanceToLine(const Vector3 & p, const Vector3 & v0, const Vector3 & v1, qeal& dist, qeal & t, Vector3 & fp);

#endif