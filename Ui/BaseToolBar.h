#pragma once
#ifndef BASE_TOOL_BAR_H
#define BASE_TOOL_BAR_H
#include <QtWidgets/qapplication.h>
#include <QtWidgets>
class BaseToolBar :public QToolBar
{
public:
	BaseToolBar();

public:
	//
	QIcon* openFileIcon;
	QAction *openFileAction;

	QIcon* saveFileIcon;
	QAction *saveFileAction;

	QIcon* renderLineModeIcon;
	QAction *renderLineModeAction;

	QIcon* selectSurfacePointsIcon;
	QAction *selectSurfacePointsIAction;

	QIcon* selectTetPointsIcon;
	QAction *selectTetPointsAction;
};

#endif // !BASE_TOOL_BAR_H

