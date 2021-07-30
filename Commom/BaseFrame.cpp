#include "BaseFrame.h"
void BaseFrameConstraint::freeFrame()
{
	select_frame.clear();
}
void BaseFrameConstraint::setFrame(BaseFrame * f)
{
	select_frame.append(f);
}
void BaseFrameConstraint::constrainTranslation(qglviewer::Vec& translation, qglviewer::Frame* const frame)
{
	for (QList<localFrame*>::iterator it = select_frame.begin(), end = select_frame.end(); it != end; ++it)
	{
		(*it)->translate(translation);
	}
}

void BaseFrameConstraint::constrainRotation(qglviewer::Quaternion& rotation, qglviewer::Frame* const frame)
{
	const qglviewer::Vec worldAxis = frame->inverseTransformOf(rotation.axis());
	const qglviewer::Vec pos = frame->position();

	const float angle = rotation.angle();

	for (QList<localFrame*>::iterator it = select_frame.begin(), end = select_frame.end(); it != end; ++it)
	{
		qglviewer::Quaternion qVertex((*it)->transformOf(worldAxis), angle);

		(*it)->rotate(qVertex);
		qglviewer::Quaternion qWorld(worldAxis, angle);
		(*it)->setPosition(pos + qWorld.rotate((*it)->position() - pos));
	}
}