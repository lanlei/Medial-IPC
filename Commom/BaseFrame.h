#ifndef BASE_FRAME_H
#define BASE_FRAME_H
#include <QGLViewer/qglviewer.h>
#include <QGLViewer/vec.h>
#include <QGLViewer/frame.h>
#include <QGLViewer/manipulatedFrame.h>
#include <QGLViewer/constraint.h>
#include <vector>
#include <set>
#include <memory>

enum ConstraintType
{
	Selected = 0,
	Fixed = 1,
	Handled = 2,
	Free = 3
};

class  BaseFrame : public qglviewer::ManipulatedFrame
{
public:
	BaseFrame(BaseFrame* parent = nullptr) : cType(Free)
	{
		if (parent != nullptr)
			this->setReferenceFrame(parent);
	}
	unsigned int id;
	ConstraintType cType;
};

typedef BaseFrame vertexFrame;
typedef BaseFrame edgeFrame;
typedef BaseFrame faceFrame;
typedef BaseFrame localFrame;
typedef BaseFrame OrentationFrame;

struct FrameIndexComp
{
	bool operator() (const BaseFrame &a, BaseFrame &b)
	{
		return a.id < b.id;
	}
};

typedef std::set<BaseFrame, FrameIndexComp> FrameSet;

class BaseFrameConstraint : public qglviewer::Constraint
{
public:
	BaseFrameConstraint() {}
	~BaseFrameConstraint() {}
	void freeFrame();
	void setFrame(BaseFrame* f);

	virtual void constrainTranslation(qglviewer::Vec& translation, qglviewer::Frame* const frame);
	virtual void constrainRotation(qglviewer::Quaternion& rotation, qglviewer::Frame* const frame);

protected:
	QList<BaseFrame*>	select_frame;
};

#endif

#pragma once
