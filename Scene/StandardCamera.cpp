#include "StandardCamera.h"

StandardCamera::StandardCamera()
{
	_halfWidth = 10;
	_halfHight = 10;
	_zNear = 0.01;
	_zFar = 100;

	_locked = false;
	_showFov = false;
	setType(qglviewer::Camera::PERSPECTIVE);
	_contraint = new qglviewer::WorldConstraint();
}

StandardCamera::~StandardCamera()
{
}

void StandardCamera::draw(bool drawFarPlane, qreal scale) const
{
	if (!_showFov)
		return;
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLineWidth(4.0);
	glColor4f(0.0, 0.0, 1.0, 0.5);
	qglviewer::Camera::draw(drawFarPlane, scale);
	glEnable(GL_LIGHTING);
	glDisable(GL_BLEND);
}

qreal StandardCamera::zNear() const
{
	return _zNear;
}

qreal StandardCamera::zFar() const
{
	return _zFar;
}

void StandardCamera::getOrthoWidthHeight(GLdouble & halfWidth, GLdouble & halfHeight) const
{
	halfWidth = _halfWidth;
	halfHeight = _halfHight;
}

void StandardCamera::setOrtho(const GLdouble halfWidth, const GLdouble halfHeight, const GLdouble zNear, const GLdouble zFar)
{
	_halfWidth = halfWidth;
	_halfWidth = halfHeight;
	_zNear = zNear;
	_zFar = zFar;
}

QMatrix4x4 StandardCamera::getViewQMatrix()
{
	GLdouble m[16];
	getModelViewMatrix(m);
	QMatrix4x4 matrix = QMatrix4x4(m[0], m[4], m[8], m[12],
		m[1], m[5], m[9], m[13],
		m[2], m[6], m[10], m[14],
		m[3], m[7], m[11], m[15]);
	return matrix;
}

QMatrix4x4 StandardCamera::getProjectionQMatrix()
{
	GLdouble m[16];
	getProjectionMatrix(m);
	QMatrix4x4 matrix = QMatrix4x4(m[0], m[4], m[8], m[12],
		m[1], m[5], m[9], m[13],
		m[2], m[6], m[10], m[14],
		m[3], m[7], m[11], m[15]);
	return matrix;
}

QMatrix4x4 StandardCamera::getProjectionViewQMatrix()
{
	QMatrix4x4 v, p;
	v = getViewQMatrix();
	p = getProjectionQMatrix();
	return p * v;
}

void StandardCamera::enableOrthoMode(bool enable)
{
	if (enable)
		setType(Camera::ORTHOGRAPHIC);
	else setType(Camera::PERSPECTIVE);
}

void StandardCamera::enableLockCamera(bool enable)
{
	_locked = enable;
	if (_locked)
	{
		_contraint->setRotationConstraintType(qglviewer::WorldConstraint::FORBIDDEN);
		_contraint->setTranslationConstraintType(qglviewer::WorldConstraint::FORBIDDEN);
	}
	else
	{
		_contraint->setRotationConstraintType(qglviewer::WorldConstraint::FREE);
		_contraint->setTranslationConstraintType(qglviewer::WorldConstraint::FREE);
	}
	qglviewer::ManipulatedCameraFrame * f = frame();
	f->setConstraint(_contraint);
}

void StandardCamera::enableShowCameraFOV(bool enable)
{
	_showFov = enable;
}


