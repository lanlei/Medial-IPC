#pragma once
#ifndef FINITE_ELEMENT_METHOD_REDUCED_FRAME_H
#define FINITE_ELEMENT_METHOD_REDUCED_FRAME_H
#include "MatrixCore.h"

namespace FiniteElementMethod
{
	typedef enum {
		STATIC = 0,
		LINEAR = 1,
		Quadratic = 2,
		Translation = 3
	} ReducedFrameType;

	class ReducedFrame
	{
	public:
		ReducedFrame(int frameId, qeal* p, qeal* X, qeal* preX, qeal* Xtilde, qeal* Vel, qeal* preVel, qeal* Acc, qeal* preAcc, int offset)
		{
			_type = STATIC;
			_frameId = frameId;
			_offset = offset;
			_dim = 0;
			_p = p;

			_X = X + offset;
			_preX = preX + offset;
			_Xtilde = Xtilde + offset;
			_Vel = Vel + offset;
			_preVel = preVel + offset;
			_Acc = Acc + offset;
			_preAcc = preAcc + offset;

			_x = p[0];
			_y = p[1];
			_z = p[2];

			_lastX = _x, _lastY = _y, _lastZ = _y;
		}

		ReducedFrame(qeal* p, int offset = -1)
		{
			_type = STATIC;
			_frameId = offset;
			_offset = offset;
			_dim = 0;
			_p = p;

			_X = nullptr;
			_preX = nullptr;
			_Vel = nullptr;
			_preVel = nullptr;
			_Acc = nullptr;
			_preAcc = nullptr;

			_x = p[0];
			_y = p[1];
			_z = p[2];

			_lastX = _x, _lastY = _y, _lastZ = _z;
		}


		ReducedFrame(qeal* p, qeal* X, qeal* Xtilde, int offset = -1)
		{
			_type = STATIC;
			_frameId = offset;
			_offset = offset;
			_dim = 0;
			_p = p;

			_X = X + offset;
			_preX = nullptr;
			_Xtilde = Xtilde + offset;
			_Vel = nullptr;
			_preVel = nullptr;
			_Acc = nullptr;
			_preAcc = nullptr;

			_x = p[0];
			_y = p[1];
			_z = p[2];

			_lastX = _x, _lastY = _y, _lastZ = _z;
		}



		ReducedFrameType getFrameType() { return _type; }
		int getFrameId() { return _frameId; }
		int getOffset() { return _offset; }
		int getDim() { return _dim; }
		qeal* getDofs() { return _q; }

		qeal* getP() { return _p; }
		void getOriginalP(qeal* oriP) 
		{
			oriP[0] = _x; oriP[1] = _y; oriP[2] = _z;
		}

		void getLastPos(qeal* lastP)
		{
			lastP[0] = _lastX; lastP[1] = _lastY; lastP[2] = _lastZ;
		}
		qeal* getCurrentX() { return _X; }
		qeal* getCurrentVel() { return _Vel; }
		qeal* getCurrentAcc() { return _Acc; }
		qeal* getCurrentXtilde() { return _Xtilde; }
		qeal* getPreviousX() { return _preX; }
		qeal* getPreviousVel() { return _preVel; }
		qeal* getPreviousAcc() { return _preAcc; }

		virtual void projectFullspaceP(qeal* P, qeal w = 1.0, qeal* oriP = nullptr);
		virtual void projectFullspacePreP(qeal* P, qeal w = 1.0, qeal* oriP = nullptr);
		virtual void projectFullspaceX(qeal* X, qeal w = 1.0, qeal* oriP = nullptr);
		virtual void projectFullspacePreX(qeal* X, qeal w = 1.0, qeal* oriP = nullptr);
		virtual void projectFullspaceXtilde(qeal* Xtilde, qeal* oriP = nullptr, qeal w = 1.0);

		virtual MatrixX getUMatrix(qeal* oriP = nullptr, const qeal w = 1.0);
		virtual void transform(){}

	//protected:
	public:
		ReducedFrameType _type;
		int _frameId;
		qeal* _q;
		qeal* _p;

		qeal* _X;
		qeal* _preX;
		qeal* _Xtilde;
		
		qeal* _Vel;
		qeal* _preVel;

		qeal* _Acc;
		qeal* _preAcc;

		int _offset;
		int _dim;
		qeal _x, _y, _z;
		qeal _lastX, _lastY, _lastZ;
	};

	class LinearReducedFrame: public ReducedFrame
	{
	public:
		LinearReducedFrame(int frameId, qeal* p, qeal* X, qeal* preX, qeal* Xtilde, qeal* Vel, qeal* preVel, qeal* Acc, qeal* preAcc, int offset = -1) :ReducedFrame(frameId, p, X, preX, Xtilde, Vel, preVel, Acc, preAcc, offset)
		{
			_type = LINEAR;
			_dim = 12;
		}

		virtual void projectFullspaceX(qeal* X, qeal w = 1.0, qeal* oriP = nullptr);
		virtual void projectFullspacePreX(qeal* X, qeal w = 1.0, qeal* oriP = nullptr);
		virtual void projectFullspaceXtilde(qeal* Xtilde, qeal* oriP = nullptr, qeal w = 1.0);

		virtual MatrixX getUMatrix(qeal* oriP = nullptr, const qeal w = 1.0);
		virtual void transform();
	};

	class QuadraticReducedFrame : public ReducedFrame
	{
	public:
		QuadraticReducedFrame(int frameId, qeal* p, qeal* X, qeal* preX, qeal* Xtilde, qeal* Vel, qeal* preVel, qeal* Acc, qeal* preAcc, int offset = -1) :ReducedFrame(frameId, p, X, preX, Xtilde, Vel, preVel, Acc, preAcc, offset)
		{
			_type = Quadratic;
			_dim = 30;

			_xx = p[0] * p[0];
			_yy = p[1] * p[1];
			_zz = p[2] * p[2];
			_xy = p[0] * p[1];
			_yz = p[1] * p[2];
			_xz = p[0] * p[2];
		}
		virtual void projectFullspaceX(qeal* X, qeal w = 1.0, qeal* oriP = nullptr);
		virtual void projectFullspacePreX(qeal* X, qeal w = 1.0, qeal* oriP = nullptr);
		virtual void projectFullspaceXtilde(qeal* Xtilde, qeal* oriP = nullptr, qeal w = 1.0);

		virtual MatrixX getUMatrix(qeal* oriP = nullptr, const qeal w = 1.0);
		virtual void transform();
	protected:
		qeal _xx, _yy, _zz;
		qeal _xy, _yz, _xz;
	};


	// for debug
	class TranslationReducedFrame : public ReducedFrame
	{
	public:
		TranslationReducedFrame(int frameId, qeal* p, qeal* X, qeal* preX, qeal* Xtilde, qeal* Vel, qeal* preVel, qeal* Acc, qeal* preAcc, int offset = -1) :ReducedFrame(frameId, p, X, preX, Xtilde, Vel, preVel, Acc, preAcc, offset)
		{
			_type = Translation;
			_dim = 3;
		}

		virtual void projectFullspaceX(qeal* X, qeal w = 1.0, qeal* oriP = nullptr);
		virtual void projectFullspacePreX(qeal* X, qeal w = 1.0, qeal* oriP = nullptr);
		virtual void projectFullspaceXtilde(qeal* Xtilde, qeal* oriP = nullptr, qeal w = 1.0);

		virtual MatrixX getUMatrix(qeal* oriP = nullptr, const qeal w = 1.0);
		virtual void transform();
	};


}

#endif
