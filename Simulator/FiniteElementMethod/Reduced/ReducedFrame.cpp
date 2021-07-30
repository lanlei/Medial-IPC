#include "ReducedFrame.h"
#include <iostream>

void FiniteElementMethod::ReducedFrame::projectFullspaceP(qeal * P, qeal w, qeal * oriP)
{
	qeal d[3];
	projectFullspaceX(d, w, oriP);
	if (oriP != nullptr)
	{
		P[0] = oriP[0] + d[0];
		P[1] = oriP[1] + d[1];
		P[2] = oriP[2] + d[2];
	}
	else
	{
		P[0] = _x + d[0];
		P[1] = _y + d[1];
		P[2] = _z + d[2];
	}
}

void FiniteElementMethod::ReducedFrame::projectFullspacePreP(qeal * P, qeal w, qeal * oriP)
{
	qeal d[3];
	projectFullspacePreX(d, w, oriP);
	if (oriP != nullptr)
	{
		P[0] = oriP[0] + d[0];
		P[1] = oriP[1] + d[1];
		P[2] = oriP[2] + d[2];
	}
	else
	{
		P[0] = _x + d[0];
		P[1] = _y + d[1];
		P[2] = _z + d[2];
	}
}

void FiniteElementMethod::ReducedFrame::projectFullspaceX(qeal * X, qeal w, qeal * oriP)
{
	X[0] = 0.0;
	X[1] = 0.0;
	X[2] = 0.0;
}

void FiniteElementMethod::ReducedFrame::projectFullspacePreX(qeal * X, qeal w, qeal * oriP)
{
	X[0] = 0.0;
	X[1] = 0.0;
	X[2] = 0.0;
}

void FiniteElementMethod::ReducedFrame::projectFullspaceXtilde(qeal * Xtilde, qeal * oriP, qeal w)
{
	Xtilde[0] = 0.0;
	Xtilde[1] = 0.0;
	Xtilde[2] = 0.0;
}

MatrixX FiniteElementMethod::ReducedFrame::getUMatrix(qeal * oriP, const qeal w)
{	
	MatrixX u(3, 1);
	u.setZero();
	return u;
}

void  FiniteElementMethod::LinearReducedFrame::projectFullspaceX(qeal* X, qeal w, qeal* oriP)
{
	qeal dx, dy, dz;
	if (oriP != nullptr)
	{
		dx = oriP[0] * _X[0] + oriP[1] * _X[3] + oriP[2] * _X[6] + _X[9];
		dy = oriP[0] * _X[1] + oriP[1] * _X[4] + oriP[2] * _X[7] + _X[10];
		dz = oriP[0] * _X[2] + oriP[1] * _X[5] + oriP[2] * _X[8] + _X[11];

		X[0] = w * dx;
		X[1] = w * dy;
		X[2] = w * dz;
		return;
	}

	dx = _x * _X[0] + _y * _X[3] + _z * _X[6] + _X[9];
	dy = _x * _X[1] + _y * _X[4] + _z * _X[7] + _X[10];
	dz = _x * _X[2] + _y * _X[5] + _z * _X[8] + _X[11];

	X[0] = w * dx;
	X[1] = w * dy;
	X[2] = w * dz;
}

void FiniteElementMethod::LinearReducedFrame::projectFullspacePreX(qeal * X, qeal w, qeal * oriP)
{
	qeal dx, dy, dz;
	if (oriP != nullptr)
	{
		dx = oriP[0] * _preX[0] + oriP[1] * _preX[3] + oriP[2] * _preX[6] + _preX[9];
		dy = oriP[0] * _preX[1] + oriP[1] * _preX[4] + oriP[2] * _preX[7] + _preX[10];
		dz = oriP[0] * _preX[2] + oriP[1] * _preX[5] + oriP[2] * _preX[8] + _preX[11];

		X[0] = w * dx;
		X[1] = w * dy;
		X[2] = w * dz;
		return;
	}

	dx = _x * _preX[0] + _y * _preX[3] + _z * _preX[6] + _preX[9];
	dy = _x * _preX[1] + _y * _preX[4] + _z * _preX[7] + _preX[10];
	dz = _x * _preX[2] + _y * _preX[5] + _z * _preX[8] + _preX[11];

	X[0] = w * dx;
	X[1] = w * dy;
	X[2] = w * dz;
}

void  FiniteElementMethod::LinearReducedFrame::projectFullspaceXtilde(qeal* Xtilde, qeal* oriP, qeal w)
{
	qeal dx, dy, dz;
	if (oriP != nullptr)
	{
		dx = oriP[0] * _Xtilde[0] + oriP[1] * _Xtilde[3] + oriP[2] * _Xtilde[6] + _Xtilde[9];
		dy = oriP[0] * _Xtilde[1] + oriP[1] * _Xtilde[4] + oriP[2] * _Xtilde[7] + _Xtilde[10];
		dz = oriP[0] * _Xtilde[2] + oriP[1] * _Xtilde[5] + oriP[2] * _Xtilde[8] + _Xtilde[11];

		Xtilde[0] = w * dx;
		Xtilde[1] = w * dy;
		Xtilde[2] = w * dz;
		return;
	}

	dx = _x * _Xtilde[0] + _y * _Xtilde[3] + _z * _Xtilde[6] + _Xtilde[9];
	dy = _x * _Xtilde[1] + _y * _Xtilde[4] + _z * _Xtilde[7] + _Xtilde[10];
	dz = _x * _Xtilde[2] + _y * _Xtilde[5] + _z * _Xtilde[8] + _Xtilde[11];

	Xtilde[0] = w * dx;
	Xtilde[1] = w * dy;
	Xtilde[2] = w * dz;
}

MatrixX  FiniteElementMethod::LinearReducedFrame::getUMatrix(qeal* oriP, const qeal w)
{
	MatrixX u(3, 12);
	u.setZero();
	if (_offset < 0 || _X == nullptr)
		return u;

	qeal x, y, z, t;
	if (oriP == nullptr)
	{
		x = _x * w;
		y = _y * w;
		z = _z * w;
	}
	else
	{
		x = oriP[0] * w;
		y = oriP[1] * w;
		z = oriP[2] * w;
	}

	t = w;

	u.data()[0] = x; u.data()[9] = y; u.data()[18] = z;
	u.data()[4] = x; u.data()[13] = y; u.data()[22] = z;
	u.data()[8] = x; u.data()[17] = y; u.data()[26] = z;

	u.data()[27] = t; u.data()[31] = t; u.data()[35] = t;

	return u;
}

void  FiniteElementMethod::LinearReducedFrame::transform()
{
	if (_offset < 0 || _X == nullptr) return;

	_lastX = _p[0];
	_lastY = _p[1];
	_lastZ = _p[2];

	qeal dx = _x * _X[0] + _y * _X[3] + _z * _X[6] + _X[9];
	qeal dy = _x * _X[1] + _y * _X[4] + _z * _X[7] + _X[10];
	qeal dz = _x * _X[2] + _y * _X[5] + _z * _X[8] + _X[11];

	_p[0] = _x + dx;
	_p[1] = _y + dy;
	_p[2] = _z + dz;
}

void  FiniteElementMethod::QuadraticReducedFrame::projectFullspaceX(qeal* X, qeal w, qeal* oriP)
{
	qeal dx, dy, dz;
	if (oriP != nullptr)
	{
		qeal x, y, z, xx, yy, zz, xy, yz, xz;
		x = oriP[0]; y = oriP[1]; z = oriP[2];
		xx = x * x; yy = y * y; zz = z * z;
		xy = x * y; yz = y * z; xz = x * z;

		dx = x * _X[0] + y * _X[3] + z * _X[6] + xx * _X[9] + yy * _X[12] + zz * _X[15] + xy * _X[18] + yz * _X[21] + xz * _X[24] + _X[27];
		dy = x * _X[1] + y * _X[4] + z * _X[7] + xx * _X[10] + yy * _X[13] + zz * _X[16] + xy * _X[19] + yz * _X[22] + xz * _X[25] + _X[28];
		dz = x * _X[2] + y * _X[5] + z * _X[8] + xx * _X[11] + yy * _X[14] + zz * _X[17] + xy * _X[20] + yz * _X[23] + xz * _X[26] + _X[29];

		X[0] = w * dx;
		X[1] = w * dy;
		X[2] = w * dz;
		return;
	}

	dx = _x * _X[0] + _y * _X[3] + _z * _X[6] + _xx * _X[9] + _yy * _X[12] + _zz * _X[15] + _xy * _X[18] + _yz * _X[21] + _xz * _X[24] + _X[27];
	dy = _x * _X[1] + _y * _X[4] + _z * _X[7] + _xx * _X[10] + _yy * _X[13] + _zz * _X[16] + _xy * _X[19] + _yz * _X[22] + _xz * _X[25] + _X[28];
	dz = _x * _X[2] + _y * _X[5] + _z * _X[8] + _xx * _X[11] + _yy * _X[14] + _zz * _X[17] + _xy * _X[20] + _yz * _X[23] + _xz * _X[26] + _X[29];

	X[0] = w * dx;
	X[1] = w * dy;
	X[2] = w * dz;
}

void FiniteElementMethod::QuadraticReducedFrame::projectFullspacePreX(qeal * X, qeal w, qeal * oriP)
{
	qeal dx, dy, dz;
	if (oriP != nullptr)
	{
		qeal x, y, z, xx, yy, zz, xy, yz, xz;
		x = oriP[0]; y = oriP[1]; z = oriP[2];
		xx = x * x; yy = y * y; zz = z * z;
		xy = x * y; yz = y * z; xz = x * z;

		dx = x * _preX[0] + y * _preX[3] + z * _preX[6] + xx * _preX[9] + yy * _preX[12] + zz * _preX[15] + xy * _preX[18] + yz * _preX[21] + xz * _preX[24] + _preX[27];
		dy = x * _preX[1] + y * _preX[4] + z * _preX[7] + xx * _preX[10] + yy * _preX[13] + zz * _preX[16] + xy * _preX[19] + yz * _preX[22] + xz * _preX[25] + _preX[28];
		dz = x * _preX[2] + y * _preX[5] + z * _preX[8] + xx * _preX[11] + yy * _preX[14] + zz * _preX[17] + xy * _preX[20] + yz * _preX[23] + xz * _preX[26] + _preX[29];

		X[0] = w * dx;
		X[1] = w * dy;
		X[2] = w * dz;
		return;
	}


	dx = _x * _preX[0] + _y * _preX[3] + _z * _preX[6] + _xx * _preX[9] + _yy * _preX[12] + _zz * _preX[15] + _xy * _preX[18] + _yz * _preX[21] + _xz * _preX[24] + _preX[27];
	dy = _x * _preX[1] + _y * _preX[4] + _z * _preX[7] + _xx * _preX[10] + _yy * _preX[13] + _zz * _preX[16] + _xy * _preX[19] + _yz * _preX[22] + _xz * _preX[25] + _preX[28];
	dz = _x * _preX[2] + _y * _preX[5] + _z * _preX[8] + _xx * _preX[11] + _yy * _preX[14] + _zz * _preX[17] + _xy * _preX[20] + _yz * _preX[23] + _xz * _preX[26] + _preX[29];

	X[0] = w * dx;
	X[1] = w * dy;
	X[2] = w * dz;
}

void  FiniteElementMethod::QuadraticReducedFrame::projectFullspaceXtilde(qeal* Xtilde, qeal* oriP, qeal w)
{
	qeal dx, dy, dz;
	if (oriP != nullptr)
	{
		qeal x, y, z, xx, yy, zz, xy, yz, xz;
		x = oriP[0]; y = oriP[1]; z = oriP[2];
		xx = x * x; yy = y * y; zz = z * z;
		xy = x * y; yz = y * z; xz = x * z;

		dx = x * _Xtilde[0] + y * _X[3] + z * _Xtilde[6] + xx * _Xtilde[9] + yy * _Xtilde[12] + zz * _Xtilde[15] + xy * _Xtilde[18] + yz * _Xtilde[21] + xz * _Xtilde[24] + _Xtilde[27];
		dy = x * _Xtilde[1] + y * _X[4] + z * _Xtilde[7] + xx * _Xtilde[10] + yy * _Xtilde[13] + zz * _Xtilde[16] + xy * _Xtilde[19] + yz * _Xtilde[22] + xz * _Xtilde[25] + _Xtilde[28];
		dz = x * _Xtilde[2] + y * _X[5] + z * _Xtilde[8] + xx * _Xtilde[11] + yy * _Xtilde[14] + zz * _Xtilde[17] + xy * _Xtilde[20] + yz * _Xtilde[23] + xz * _Xtilde[26] + _Xtilde[29];

		Xtilde[0] = w * dx;
		Xtilde[1] = w * dy;
		Xtilde[2] = w * dz;
		return;
	}

	dx = _x * _Xtilde[0] + _y * _Xtilde[3] + _z * _Xtilde[6] + _xx * _Xtilde[9] + _yy * _Xtilde[12] + _zz * _Xtilde[15] + _xy * _Xtilde[18] + _yz * _Xtilde[21] + _xz * _Xtilde[24] + _Xtilde[27];
	dy = _x * _Xtilde[1] + _y * _Xtilde[4] + _z * _Xtilde[7] + _xx * _Xtilde[10] + _yy * _Xtilde[13] + _zz * _Xtilde[16] + _xy * _Xtilde[19] + _yz * _Xtilde[22] + _xz * _Xtilde[25] + _Xtilde[28];
	dz = _x * _Xtilde[2] + _y * _Xtilde[5] + _z * _Xtilde[8] + _xx * _Xtilde[11] + _yy * _Xtilde[14] + _zz * _Xtilde[17] + _xy * _Xtilde[20] + _yz * _Xtilde[23] + _xz * _Xtilde[26] + _Xtilde[29];

	Xtilde[0] = w * dx;
	Xtilde[1] = w * dy;
	Xtilde[2] = w * dz;
}

MatrixX  FiniteElementMethod::QuadraticReducedFrame::getUMatrix(qeal* oriP, const qeal w)
{
	MatrixX u(3, 30);
	u.setZero();
	if (_offset < 0 || _X == nullptr)
		return u;

	qeal x, y, z, xx, yy, zz, xy, yz, xz, t;
	if (oriP == nullptr)
	{
		x = _x * w;
		y = _y * w;
		z = _z * w;
		xx = _xx* w;
		yy = _yy * w;
		zz = _zz * w;
		xy = _xy * w;
		yz = _yz * w;
		xz = _xz * w;
	}
	else
	{
		x = oriP[0] * w;
		y = oriP[1] * w;
		z = oriP[2] * w;
		xx = oriP[0] * oriP[0] * w;
		yy = oriP[1] * oriP[1] * w;
		zz = oriP[2] * oriP[2] * w;
		xy = oriP[0] * oriP[1] * w;
		yz = oriP[1] * oriP[2] * w;
		xz = oriP[0] * oriP[2] * w;
	}
	t = w;

	u.data()[0] = x; u.data()[9] = y; 	u.data()[18] = z;
	u.data()[4] = x; u.data()[13] = y; u.data()[22] = z;
	u.data()[8] = x; u.data()[17] = y; u.data()[26] = z;

	u.data()[27] = xx; u.data()[36] = yy; u.data()[45] = zz;
	u.data()[31] = xx; u.data()[40] = yy; u.data()[49] = zz;
	u.data()[35] = xx; u.data()[44] = yy; u.data()[53] = zz;

	u.data()[54] = xy; u.data()[63] = yz; u.data()[72] = xz;
	u.data()[58] = xy; u.data()[67] = yz; u.data()[76] = xz;
	u.data()[62] = xy; u.data()[71] = yz; u.data()[80] = xz;

	u.data()[81] = t; u.data()[85] = t; u.data()[89] = t;

	return u;
}

void  FiniteElementMethod::QuadraticReducedFrame::transform()
{
	if (_offset < 0 || _X == nullptr) return;

	_lastX = _p[0];
	_lastY = _p[1];
	_lastZ = _p[2];

	qeal dx = _x * _X[0] + _y * _X[3] + _z * _X[6] + _xx * _X[9] + _yy * _X[12] + _zz * _X[15] + _xy * _X[18] + _yz * _X[21] + _xz * _X[24] + _X[27];
	qeal dy = _x * _X[1] + _y * _X[4] + _z * _X[7] + _xx * _X[10] + _yy * _X[13] + _zz * _X[16] + _xy * _X[19] + _yz * _X[22] + _xz * _X[25] + _X[28];
	qeal dz = _x * _X[2] + _y * _X[5] + _z * _X[8] + _xx * _X[11] + _yy * _X[14] + _zz * _X[17] + _xy * _X[20] + _yz * _X[23] + _xz * _X[26] + _X[29];

	_p[0] = _x + dx;
	_p[1] = _y + dy;
	_p[2] = _z + dz;
}

void FiniteElementMethod::TranslationReducedFrame::projectFullspaceX(qeal * X, qeal w, qeal * oriP)
{
	qeal dx, dy, dz;
	dx = _X[0];
	dy = _X[1];
	dz = _X[2];

	X[0] = w * dx;
	X[1] = w * dy;
	X[2] = w * dz;
}

void FiniteElementMethod::TranslationReducedFrame::projectFullspacePreX(qeal * X, qeal w, qeal * oriP)
{
	qeal dx, dy, dz;
	dx = _preX[0];
	dy = _preX[1];
	dz = _preX[2];

	X[0] = w * dx;
	X[1] = w * dy;
	X[2] = w * dz;
}

void FiniteElementMethod::TranslationReducedFrame::projectFullspaceXtilde(qeal * Xtilde, qeal * oriP, qeal w)
{
	qeal dx, dy, dz;
	dx = _Xtilde[0];
	dy = _Xtilde[1];
	dz = _Xtilde[2];

	Xtilde[0] = w * dx;
	Xtilde[1] = w * dy;
	Xtilde[2] = w * dz;
}

MatrixX FiniteElementMethod::TranslationReducedFrame::getUMatrix(qeal * oriP, const qeal w)
{
	MatrixX u(3, 3);
	u.setZero();
	if (_offset < 0 || _X == nullptr)
		return u;
	u.setIdentity();
	u *= w;
	return u;
}

void FiniteElementMethod::TranslationReducedFrame::transform()
{
	if (_offset < 0 || _X == nullptr) return;

	_lastX = _p[0];
	_lastY = _p[1];
	_lastZ = _p[2];

	qeal dx = _X[0];
	qeal dy = _X[1];
	qeal dz =_X[2];

	_p[0] = _x + dx;
	_p[1] = _y + dy;
	_p[2] = _z + dz;
}
