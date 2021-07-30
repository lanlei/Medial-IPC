#include "BaseToolBar.h"

BaseToolBar::BaseToolBar()
{
	openFileIcon = new QIcon("./icon/open_file.png");
	openFileAction = new QAction(this);
	openFileAction->setIcon(*openFileIcon);

	saveFileIcon = new QIcon("./icon/save_file.png");
	saveFileAction = new QAction(this);
	saveFileAction->setIcon(*saveFileIcon);

	renderLineModeIcon = new QIcon("./icon/render_line_mode.png");
	renderLineModeAction = new QAction(this);
	renderLineModeAction->setCheckable(true);
	renderLineModeAction->setIcon(*renderLineModeIcon);

	selectSurfacePointsIcon = new QIcon("./icon/select_surface_point.png");
	selectSurfacePointsIAction = new QAction(this);
	selectSurfacePointsIAction->setCheckable(true);
	selectSurfacePointsIAction->setIcon(*selectSurfacePointsIcon);

	selectTetPointsIcon = new QIcon("./icon/select_tet_point.png");
	selectTetPointsAction = new QAction(this);
	selectTetPointsAction->setCheckable(true);
	selectTetPointsAction->setIcon(*selectTetPointsIcon);

	addAction(openFileAction);
	addAction(saveFileAction);
	saveFileAction->setEnabled(false);
	addAction(renderLineModeAction);
	addAction(selectSurfacePointsIAction);
	addAction(selectTetPointsAction);
}
