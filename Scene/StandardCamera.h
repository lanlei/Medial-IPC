#pragma once
#ifndef STANDARD_CAMERA_H
#define STANDARD_CAMERA_H
#include <QGLViewer/manipulatedCameraFrame.h>
#include <QGLViewer/camera.h>
#include <QMatrix4x4>
#include "MatrixCore.h"

class StandardCamera : public qglviewer::Camera
{
public:
	StandardCamera();
	~StandardCamera();

	virtual void 	draw(bool drawFarPlane = true, qreal scale = 1.0) const;
	virtual qreal zNear() const;
	virtual qreal zFar() const;

	virtual void getOrthoWidthHeight(GLdouble &halfWidth,
		GLdouble &halfHeight) const;

	void setOrtho(const GLdouble halfWidth,
		const GLdouble halfHeight, const GLdouble zNear,
		const GLdouble zFar);

	QMatrix4x4 getViewQMatrix();
	QMatrix4x4 getProjectionQMatrix();
	QMatrix4x4 getProjectionViewQMatrix();

	void enableOrthoMode(bool enable = false);
	bool queryOrthoMode() { return type() == qglviewer::Camera::ORTHOGRAPHIC ? 1 : 0; }

	void enableLockCamera(bool enable = false);
	bool queryLockCamera() { return _locked; }

	void enableShowCameraFOV(bool enable = false);
	bool queryShowFOV() { return _showFov; }
protected:
	qeal _halfWidth;
	qeal _halfHight;
	qeal _zNear;
	qeal _zFar;

	qglviewer::WorldConstraint* _contraint;	
	bool _locked;
	bool _showFov;
};

#endif // !STANDARD_CAMERA_H