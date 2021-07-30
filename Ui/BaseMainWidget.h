#pragma once
#ifndef BASE_MAIN_WIDGET_H
#define BASE_MAIN_WIDGET_H
#include <QtWidgets/QMainWindow>
#include <QtWidgets/qapplication.h>
#include <QtWidgets>
#include "BaseBottomWidget.h"
#include "BaseRightWidget.h"
#include "MainLayout.h"
#include "Scene\BaseScene.h"
#include "BaseToolBar.h"
#include "Commom\tinyxml\tinyxml.h"

#include "Simulator\BaseSimulator.h"

class BaseMainWidget : public QWidget
{
	Q_OBJECT
public:
	BaseMainWidget(){
		_scene = new BaseScene();
		_bottom = new BaseBottomWidget();
		_right = new BaseRightWidget();

		_layout = new MainLayout();
		_layout->addWidget(_scene, MainLayout::CENTER);
		_layout->addWidget(_bottom, MainLayout::BOTTOM);
		//_layout->addWidget(_right, MainLayout::RIGHT);
		setLayout(_layout);
		_topToolBar = new BaseToolBar();
		_bottom->setBaseScene(_scene);
		_right->setBaseScene(_scene);
		_topToolBar->setObjectName(QStringLiteral("mainToolBar"));
		bindWidget();
	}
	~BaseMainWidget()
	{
			
	}

	BaseBottomWidget* getBottomWidget() { return _bottom; }
	void setBottomWidget(BaseBottomWidget* ptr)
	{
		if (!_bottom)
			delete _bottom;
		_bottom = ptr;
	}

	BaseRightWidget* getRightWidget() { return _right; }

	void setRightWidget(BaseRightWidget* ptr)
	{
		if (!_right)
			delete _right;
		_right = ptr;
	}
	
	BaseToolBar* getTopToolBar() { return _topToolBar; }
	void setTopToolBar(BaseToolBar* bar) { _topToolBar = bar; bindTopToolBar();}

	BaseScene* getCenteralScene() {return _scene; }

	void setBaseScene(BaseScene* ptr)
	{
		if (_scene)
		{
			_layout->remove(MainLayout::CENTER);
			_layout->removeWidget(_scene);
			free(_scene);
		}			
		_scene = ptr;
		_layout->addWidget(_scene, MainLayout::CENTER);

		bindWidget();
	}

	void connectUI(QMainWindow* w)
	{
		if (!w)
			return;
		_w = w;
		_w->setCentralWidget(this);
		_w->show();
		_w->addToolBar(Qt::TopToolBarArea, getTopToolBar());
	}

	void bindWidget();
	void bindTopToolBar();
	void bindBottomAndScene();
	void bindRightAndBottomWidget();
	void bindRightAndMainWidget();

	bool importFromConfigFile(const std::string filename);

	virtual bool readSimulatorConfigFile(TiXmlElement* subItem);
signals:
	void addFileSuccess(bool);
public slots:
	// tool bar
	void handleOpenFileAction();
	void handleSaveFileAction();
	void handleRenderLineModeAction();
	void handleSelectSurfacePointsAction();
	void handleSelectTetPointsAction();

	// bottom connect scene
	void handlePlayAnimationBtn(bool);
	void handleRecorAnimationBtn(bool);
	void handleResetAnimationBtn();
	void handleMultiViewerButtonGroup(int);

protected:
	QMainWindow* _w;
	MainLayout* _layout = nullptr;
	BaseBottomWidget* _bottom = nullptr;
	BaseRightWidget* _right = nullptr;
	BaseToolBar* _topToolBar = nullptr;
	BaseScene* _scene = nullptr;
};

#endif