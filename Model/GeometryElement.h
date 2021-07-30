#pragma once
#ifndef Geometry_Element_H
#define Geometry_Element_H

#include "TriangleMesh/BaseSurfaceMesh.h"
#include "Commom/DataConversion.h"
#include "Commom/BaseFrame.h"
#include "MatrixCore.h"
#include "Commom/GeometryComputation.h"

class Box : public BaseSurfaceMesh
{
public:
	Box(qeal halfX = 0.5, qeal halfY = 0.5, qeal halfZ = 0.5, QVector3D center = QVector3D(0.0, 0.0, 0.0), QColor uniformColor = QColor(75, 75, 217));
	void scale(qeal sx, qeal sy, qeal sz);
	void setFromBounding(const QVector3D min, const QVector3D max);
	void setCenter(QVector3D center);
protected:
	QVector3D _center;
	qeal _halfZ;
	qeal _halfY;
	qeal _halfX;
};

class Floor : public BaseSurfaceMesh
{
public:
	Floor(int xz = 30.0, QColor uniformColor = QColor(217, 217, 217));
	void scaleXZ(int xz);
	void setYPlane(qeal y);

	int getScaleXZ() {return _scaleXZ;}
	qeal getYPlane() { return _y;}
protected:
	int _scaleXZ;
	qeal _y;
};

class SphereElement
{
public:
	SphereElement() {}
	SphereElement(Vector3 cent, qeal r, localFrame* f = nullptr)
	{
		radius = r;
		center = cent;
		if (f != nullptr)
			frame.setReferenceFrame(f);
	}

	void operator=(const SphereElement& b);

	void show();
	qeal radius;
	Vector3 center;
	OrentationFrame frame;
};

typedef SphereElement Sphere;

class SplintElement
{
public:
	SplintElement() {}
	SplintElement(const Vector3  v0, const Vector3  v1, const Vector3  v2, Vector3  n = Vector3(0, 0, 0));

	void updateNormal(bool reverse = false);
	Vector3 vt[3];
	Vector3 nt;
};

class TwoSplintElements
{
public:
	TwoSplintElements(Vector3 c0, qeal r0, Vector3 c1, qeal r1, Vector3 c2, qeal r2)
	{
		Vector3 c0c1 = c1 - c0;
		Vector3 c0c2 = c2 - c0;
		Vector3 c1c2 = c2 - c1;
		qeal c0c1len = c0c1.norm();
		qeal c0c2len = c0c2.norm();
		qeal c1c2len = c1c2.norm();
		qeal dr0r1 = fabs(r0 - r1);
		qeal dr0r2 = fabs(r0 - r2);
		qeal dr1r2 = fabs(r1 - r2);

		// some spheres are concentric and there are no triangles.
		if ((c0c1len < QEAL_MIN) || (c0c2len < QEAL_MIN) || (c1c2len < QEAL_MIN))
			return;
		Vector3 norm;
		norm = c0c1.cross(c0c2);
		norm.normalize();

		// equal-radius spheres
		if ((dr0r1 < QEAL_MIN) && (dr0r2 < QEAL_MIN) && (dr1r2 < QEAL_MIN))
		{
			st1.vt[0] = c0 + norm * r0;
			st1.vt[1] = c1 + norm * r1;
			st1.vt[2] = c2 + norm * r2;
			st1.updateNormal();

			st2.vt[0] = c0 - norm * r0;
			st2.vt[1] = c1 - norm * r1;
			st2.vt[2] = c2 - norm * r2;
			st2.updateNormal();
		}
		else
		{
			// two points on the tangent plane
			Vector3 apex0, apex1;
			// two spheres are equal-radius
			if (dr0r1 < QEAL_MIN)
			{
				apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
				apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
			}
			else if (dr0r2 < QEAL_MIN)
			{
				apex0 = (r1 * c0 - r0 * c1) / (r1 - r0);
				apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
			}
			else if (dr1r2 < QEAL_MIN)
			{
				apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
				apex1 = (r0 * c1 - r1 * c0) / (r0 - r1);
			}
			else
			{
				apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
				apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
			}

			qeal distc0;
			Vector3 fp;
			qeal t;
			DistanceToLine(c0, apex0, apex1, distc0, t, fp);

			qeal sangle = r0 / distc0;
			if (fabs(sangle) > 1.0)
				return;

			qeal cangle = sqrt(1. - r0 * r0 / distc0 / distc0);
			Vector3 norfpc0(c0 - fp);
			norfpc0.normalize();
			Vector3 newnorm[2];
			newnorm[0] = norm * cangle - norfpc0 * sangle;
			newnorm[1] = -norm * cangle - norfpc0 * sangle;

			st1.vt[0] = c0 + r0 * newnorm[0];
			st1.vt[1] = c1 + r1 * newnorm[0];
			st1.vt[2] = c2 + r2 * newnorm[0];
			st1.updateNormal();

			st2.vt[0] = c0 + r0 * newnorm[1];
			st2.vt[1] = c1 + r1 * newnorm[1];
			st2.vt[2] = c2 + r2 * newnorm[1];
			st2.updateNormal(true);

		}
	}
	void show();
	SplintElement st1, st2;
};

class ConeElement
{
public:
	ConeElement(Vector3 sp0, qeal r0, Vector3 sp1, qeal r1):
		_sp0(sp0), _sp1(sp1), _r0(r0), _r1(r1)
	{
		computeCone();
	}
	ConeElement(Vector3 sp0, qeal r0, Vector3 sp1, qeal r1, Vector3 _apex, Vector3 _bpex, Vector3 _axis, qeal _base, qeal _top, qeal _height, qeal _sita, Vector3 _rot_axis, qeal _rot_angle) :
		_sp0(sp0), _sp1(sp1), _r0(r0), _r1(r1),
		apex(_apex),
		bpex(_bpex),
		axis(_axis),
		base(_base),
		top(_top),
		height(_height),
		sita(_sita),
		rot_axis(_rot_axis),
		rot_angle(_rot_angle) {}

	void computeCone();
	void show();

	Vector3 apex;
	Vector3 bpex;
	Vector3 axis;
	qeal base;
	qeal top;
	qeal height;
	qeal sita;
	Vector3 rot_axis;
	qeal rot_angle;

protected:
	Vector3 _sp0, _sp1;
	qeal _r0, _r1;
};

#endif