#include "GeometryComputation.h"

qeal TriangleArea(const Vector3& v0, const Vector3& v1, const Vector3& v2)
{
	return 0.5 * fabs(((v1 - v0).cross(v2 - v0)).norm());
}

void TriangleNormal(const Vector3& v0, const Vector3& v1, const Vector3& v2, Vector3& n)
{
	Vector3 vec1 = v1 - v0;
	Vector3 vec2 = v2 - v0;
	n = vec1.cross(vec2);
	n.normalize();
}

bool DistanceToLine(const Vector3 & p, const Vector3 & v0, const Vector3 & v1, qeal& dist, qeal & t, Vector3 & fp)
{
	Vector3 v0v1 = v1 - v0;
	Vector3 pv0 = v0 - p;
	Vector3 pv1 = v1 - p;

	qeal area = abs((v0v1.cross(pv0)).norm());

	if (!IS_QEAL_ZERO(v0v1.norm()))
	{
		dist = area / v0v1.norm();
		t = (pv0.dot(pv0) - pv0.dot(pv1)) / (pv0.dot(pv0) + pv1.dot(pv1) - 2 * pv0.dot(pv1));
		fp = (1.0 - t) * v0 + t * v1;
		return true;
	}
	else return false;
}
