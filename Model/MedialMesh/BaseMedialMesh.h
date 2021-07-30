#pragma once
#ifndef BASE_MEDIAL_MESH_H
#define BASE_MEDIAL_MESH_H

#include <memory>
#include <vector>
#include "Commom\BufferSerialization.h"
#include "Model\GeometryElement.h"
#include "BaseMedialMeshBufferPool.h"
#include <QGLViewer/vec.h>
#include <QGLViewer/Quaternion.h>
#include "MatrixCore.h"

class BaseMedialMeshBuffer
{
	typedef FragmentVectorBufferPtr <qeal> Points;
	typedef FragmentVectorBufferPtr <int> Indices;

public:
	BaseMedialMeshBuffer() :
		medialPointsNum(0),
		medialConesNum(0),
		medialSlabsNum(0),
		medialPrimitivesNum(0)
	{

	}

	virtual void refresh(BaseMedialMeshBufferPool* pool)
	{
		medialPoints.buffer = pool->medialPointsBuffer.buffer.data() + medialPoints.offset;
		medialRadius.buffer = pool->medialRadiusBuffer.buffer.data() + medialRadius.offset;
		medialPrimitives.buffer = pool->medialPrimitiveIndicesBuffer.buffer.data() + medialPrimitives.offset;

		for (int i = 0; i < medialPointsNeighborList.size(); i++)
			medialPointsNeighborList[i].buffer = pool->medialPointsNeighborListBuffer.buffer.data() + medialPointsNeighborList[i].offset;

		for (int i = 0; i < medialPointsCones.size(); i++)
			medialPointsCones[i].buffer = pool->medialPointsConesBuffer.buffer.data() + medialPointsCones[i].offset;

		for (int i = 0; i < medialPointsSlabs.size(); i++)
			medialPointsSlabs[i].buffer = pool->medialPointsConesBuffer.buffer.data() + medialPointsSlabs[i].offset;

		medialPointIdOffset = medialPoints.offset / 3;
		medialPrimitiveIdOffset = medialPrimitives.offset / 3;
	}

	int medialPointsNum, medialConesNum, medialSlabsNum, medialPrimitivesNum;
	Points medialPoints;
	Points medialRadius;
	Indices medialPrimitives;
	std::vector<Indices> medialPointsNeighborList;
	std::vector<Indices> medialPointsCones;
	std::vector<Indices> medialPointsSlabs;

	void copyMedialPointsToBuffer(qeal* buffer, int size = 0);
	void copyBufferToMedialPoints(qeal* buffer, int size = 0);

	void copyMedialRadiusToBuffer(qeal* buffer, int size = 0);
	void copyBufferToMedialRadius(qeal* buffer, int size = 0);

	Vector3 getMedialPoint(int nid);
	void setMedialPoint(int nid, qeal* p);

	qeal getMedialRadius(int nid);
	void setMedialRadius(int nid, qeal& r);
	
	Vector3i getMedialPrimitive(int pid);
	Vector2i getMedialCone(int cid);
	Vector3i getMedialSlab(int sid);

	int getMedialPointOverallId(int id);
	int getMedialPrimitiveOverallId(int pid);
protected:
	int medialPointIdOffset;
	int medialPrimitiveIdOffset;
};

class BaseMedialMesh : public BaseMedialMeshBuffer
{
public:
	BaseMedialMesh() :
		_hide(true),
		_transparent(false),
		_valid(false),
		hasBindedTetMesh(false)
	{}

	bool readMeshFromMatFormat(const std::string filename, BaseMedialMeshBufferPool* pool);
	bool writeMeshToMatFormat(const std::string filename);

	virtual void uniform(const qeal div);
	void translateMesh(const qeal x, const qeal y, const qeal z);
	void scaleMesh(qeal s, const qeal cx, const qeal cy, const qeal cz);
	void rotateMesh(const qglviewer::Quaternion oldR, const qglviewer::Quaternion R, const qeal cx, const qeal cy, const qeal cz);

	virtual void encloseObject(qeal* objPos, int dim);
	//
	void enableHide(bool enable) { _hide = enable; }
	bool isHide() { return _hide; }
	void enableTransparent(bool enable) { _transparent = enable; }
	bool isTransparent() { return _transparent; }
	void setMedialMeshValid(bool valid) { _valid = valid; }
	bool isMedialMeshValid() { return _valid; }
	void setBindedTetMesh(bool bind) { hasBindedTetMesh = bind; }
	bool isBindedTetMesh() {return hasBindedTetMesh;}

	std::vector<int> bindedTM;
	std::vector<int> bindedInverseTM;

	std::vector<Vector2i> edgeList;
	std::vector<Vector3i> faceList;

	friend class BaseMedialMeshHandle;
protected:
	bool _hide;
	bool _valid;
	bool hasBindedTetMesh;
	bool _transparent;
};

static void getNearestSphereOnMedialCone(Vector3& sc, qeal & sr, Vector3 & c11, qeal & r11, Vector3 & c12, qeal & r12, qeal & t, qeal & dist)
{
	bool inverse = false;
	if (r11 > r12)
	{
		inverse = true;
		Vector3 ctemp = c11;
		qeal rtemp = r11;
		c11 = c12;
		r11 = r12;
		c12 = ctemp;
		r12 = rtemp;
	}

	Vector3 cq = sc;
	qeal rq = sr;

	Vector3 c12c11 = c11 - c12;
	Vector3 cqc12 = c12 - cq;
	qeal R1 = r11 - r12;
	qeal A = c12c11.dot(c12c11);
	qeal D = 2.0 * (c12c11.dot(cqc12));
	qeal F = cqc12.dot(cqc12);

	t = (-1.0*(A*D - R1 * R1*D)) - sqrt((D*D - 4.0*A*F)*(R1*R1 - A)*R1*R1);
	t = t / (2.0*(A*A - A * R1*R1));

	if (t < 0.0) t = 0.0;
	if (t > 1.0) t = 1.0;

	Vector3 ct = t * c11 + (1.0 - t)*c12;
	double rt = t * r11 + (1.0 - t)*r12;

	dist = (ct - cq).norm() - (rq + rt);

	if (inverse)
	{
		t = 1.0 - t;
	}
}

static void getNearestSphereOnMedialSlab(Vector3 & sc, qeal & sr, Vector3 & c11, qeal & r11, Vector3 & c12, qeal & r12, Vector3 & c13, qeal & r13, qeal & t1, qeal & t2, qeal & dist)
{
	Vector3 cq = sc;
	qeal rq = sr;
	Vector3 c13c11 = c11 - c13;
	Vector3 c13c12 = c12 - c13;
	Vector3 cqc13 = c13 - cq;
	qeal R1 = r11 - r13;
	qeal R2 = r12 - r13;
	qeal A = c13c11.dot(c13c11);
	qeal B = c13c11.dot(c13c12);
	qeal C = c13c12.dot(c13c12);
	qeal D = 2.0 * (c13c11.dot(cqc13));
	qeal E = 2.0 * (c13c12.dot(cqc13));
	qeal F = cqc13.dot(cqc13);

	if (R1 == 0 && R2 == 0)
	{
		t1 = (B*E - 2.0*C*D) / (4.0*A*C - B * B);
		t2 = (B*D - 2.0*A*E) / (4.0*A*C - B * B);
	}
	else if (R1 != 0 && R2 == 0)
	{
		qeal H2 = -1.0*B / (2.0*C);
		qeal K2 = -1.0*E / (2.0*C);
		qeal W1 = pow((2.0*A + B * H2), 2) - 4.0*R1*R1*(A + B * H2 + C * H2*H2);
		qeal W2 = 2.0*(2.0*A + B * H2)*(B*K2 + D) - 4.0*R1*R1*(B*K2 + 2 * C*H2*K2 + D + E * H2);
		qeal W3 = pow((B*K2 + D), 2) - 4.0*R1*R1*(C*K2*K2 + E * K2 + F);
		t1 = (-W2 - sqrt(W2*W2 - 4.0*W1*W3)) / (2.0*W1);
		t2 = H2 * t1 + K2;
	}
	else
	{
		qeal L1 = 2.0*A*R2 - B * R1;
		qeal L2 = 2.0*C*R1 - B * R2;
		qeal L3 = E * R1 - D * R2;
		if (L1 == 0 && L2 != 0)
		{
			t2 = -1.0*L3 / L2;
			qeal W1 = 4.0*A*A - 4.0*R1*R1*A;
			qeal W2 = 4.0*A*(B*t2 + D) - 4.0*R1*R1*(B*t2 + D);
			qeal W3 = pow((B*t2 + D), 2) - (C*t2*t2 + E * t2 + F);
			t1 = (-W2 - sqrt(W2*W2 - 4.0*W1*W3)) / (2.0*W1);
		}
		else if (L1 != 0 && L2 == 0)
		{
			t1 = 1.0*L3 / L1;
			qeal W1 = 4.0*C*C - 4.0*R2*R2*C;
			qeal W2 = 4.0*C*(B*t1 + E) - 4.0*R2*R2*(B*t1 + E);
			qeal W3 = pow((B*t1 + E), 2) - (A*t1*t1 + D * t1 + F);
			t2 = (-W2 - sqrt(W2*W2 - 4.0*W1*W3)) / (2.0*W1);
		}
		else
		{
			qeal H3 = L2 / L1;
			qeal K3 = L3 / L1;
			qeal W1 = pow((2.0*C + B * H3), 2) - 4.0*R2*R2*(A*H3*H3 + B * H3 + C);
			qeal W2 = 2.0*(2.0*C + B * H3)*(B*K3 + E) - 4.0*R2*R2*(2.0*A*H3*K3 + B * K3 + D * H3 + E);
			qeal W3 = pow((B*K3 + E), 2) - 4.0*R2*R2*(A*K3*K3 + D * K3 + F);

			t2 = (-W2 - sqrt(W2*W2 - 4.0*W1*W3)) / (2.0*W1);
			t1 = H3 * t2 + K3;
		}
	}

	if ((t1 + t2) < 1.0 && t1 >= 0 && t1 <= 1.0 && t2 >= 0 && t2 <= 1.0)
	{
		Vector3 c = t1 * c11 + t2 * c12 + (1.0 - t1 - t2)*c13;
		qeal r = t1 * r11 + t2 * r12 + (1.0 - t1 - t2)*r13;
		dist = (c - cq).norm() - (r + rq);
		return;
	}
	else
	{
		qeal min_t1, min_t2;
		qeal  min_d;
		getNearestSphereOnMedialCone(sc, sr, c11, r11, c13, r13, t1, min_d);
		t2 = 0;
		min_t1 = t1;
		min_t2 = t2;
		getNearestSphereOnMedialCone(sc, sr, c12, r12, c13, r13, t2, dist);
		if (dist < min_d)
		{
			min_d = dist;
			min_t1 = 0;
			min_t2 = t2;
		}
		getNearestSphereOnMedialCone(sc, sr, c11, r11, c12, r12, t1, dist);
		if (dist < min_d)
		{
			min_d = dist;
			min_t1 = t1;
			min_t2 = 1.0 - t1;
		}
		t1 = min_t1;
		t2 = min_t2;
		dist = min_d;
	}
}


#endif