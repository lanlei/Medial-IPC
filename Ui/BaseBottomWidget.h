#pragma on
#ifndef BASE_BOTTOM_WIDGET_H
#define BASE_BOTTOM_WIDGET_H

#include <QtWidgets/qapplication.h>
#include <QtWidgets>
#include "ui_BottomFormWidget.h"
#include "Scene\BaseScene.h"

class BaseBottomWidget : public QWidget, public Ui::BottomWidget
{
	Q_OBJECT
public:
	BaseBottomWidget();
	~BaseBottomWidget();

	void resizeEvent(QResizeEvent* size);

	void setBaseScene(BaseScene* scene);
signals:
	void changeViewerNum(int);
	void recordAnimation(bool);
	void playAnimation(bool);
	void resetAnimation();
public slots:
	void bindMultiViewerGroup();
	void bindAnimationGroup();
	void handleMultiViewerButtonGroup(QAbstractButton *);
	void handleInitSimulator();
	void handleRecordAnimation();
	void handlePlayAnimation();
	void handleResetAnimation();
	void resetInitSimulatorStatus();
protected:
	BaseScene* _scene;
	QButtonGroup* multiViewerButtonGroup;
	QIcon* playAnimationIcon;
	QIcon* pauseAnimationIcon;
	QIcon* recordAnimationIcon;
	QIcon* resetAnimationIcon;
	QIcon* initSimulatorIcon;
};


#endif
