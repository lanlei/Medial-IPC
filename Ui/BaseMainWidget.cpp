     #include "BaseMainWidget.h"
#include "Simulator\SimulatorFactor.h"


void BaseMainWidget::bindWidget()
{
	bindTopToolBar();
	bindRightAndMainWidget();
	bindBottomAndScene();
	bindRightAndBottomWidget();

}

void BaseMainWidget::bindTopToolBar()
{
	connect(_topToolBar->openFileAction, SIGNAL(triggered()), this, SLOT(handleOpenFileAction()));
	connect(_topToolBar->saveFileAction, SIGNAL(triggered()), this, SLOT(handleSaveFileAction()));
	connect(_topToolBar->renderLineModeAction, SIGNAL(triggered()), this, SLOT(handleRenderLineModeAction()));
	connect(_topToolBar->selectSurfacePointsIAction, SIGNAL(triggered()), this, SLOT(handleSelectSurfacePointsAction()));
	connect(_topToolBar->selectTetPointsAction, SIGNAL(triggered()), this, SLOT(handleSelectTetPointsAction()));
}

void BaseMainWidget::bindRightAndMainWidget()
{
	connect(this, SIGNAL(addFileSuccess(bool)), _right, SLOT(refreshModelList()));
}

void BaseMainWidget::bindBottomAndScene()
{
	connect(_bottom, SIGNAL(recordAnimation(bool)), this, SLOT(handleRecorAnimationBtn(bool)));
	connect(_bottom, SIGNAL(playAnimation(bool)), this, SLOT(handlePlayAnimationBtn(bool)));
	connect(_bottom, SIGNAL(resetAnimation()), this, SLOT(handleResetAnimationBtn()));
	connect(_bottom, SIGNAL(changeViewerNum(int)), this, SLOT(handleMultiViewerButtonGroup(int)));
	connect(_scene, SIGNAL(needInitSim()), _bottom, SLOT(resetInitSimulatorStatus()));
}

void BaseMainWidget::bindRightAndBottomWidget()
{
	connect(_right, SIGNAL(needInitSimulator()), _bottom, SLOT(resetInitSimulatorStatus()));
	_right->setVisible(false);
}

void BaseMainWidget::handleOpenFileAction()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Select a Surface to import"), "./example/", tr(""));
	if (filename.isEmpty())
		return;
	std::string dir, name, format;
	getFilenameInfo(filename.toStdString(), dir, name, format);
	if (format == std::string("xml"))
	{
		bool flag = importFromConfigFile(filename.toStdString());
		emit addFileSuccess(flag);
	}
	else
	{
		// read binary file
	}
}

void BaseMainWidget::handleSaveFileAction()
{
	// to do
}

void BaseMainWidget::handleRenderLineModeAction()
{
	if (_scene == nullptr)
		return;
	_scene->enableRenderLineMode(_topToolBar->renderLineModeAction->isChecked());
}

void BaseMainWidget::handleSelectSurfacePointsAction()
{
	if (_scene == nullptr)
		return;
	_scene->enableSelectSurfacePointsMode(_topToolBar->selectSurfacePointsIAction->isChecked());
	if (_topToolBar->selectSurfacePointsIAction->isChecked())
		_topToolBar->selectTetPointsAction->setChecked(false);
}

void BaseMainWidget::handleSelectTetPointsAction()
{
	if (_scene == nullptr)
		return;
	if (_scene == nullptr)
		return;
	_scene->enableSelectTetNodesMode(_topToolBar->selectTetPointsAction->isChecked());
	if (_topToolBar->selectTetPointsAction->isChecked())
		_topToolBar->selectSurfacePointsIAction->setChecked(false);
}

void BaseMainWidget::handlePlayAnimationBtn(bool flag)
{
	if (flag)
		_scene->startAnimation();
	else _scene->stopAnimation();
}

void BaseMainWidget::handleRecorAnimationBtn(bool flag)
{
	_scene->getSceneStaus()->recordAnimation = flag;
}

void BaseMainWidget::handleResetAnimationBtn()
{
	_scene->getCurrentSimulator()->reset();
}

void BaseMainWidget::handleMultiViewerButtonGroup(int id)
{
	// to do
}

bool BaseMainWidget::importFromConfigFile(const std::string filename)
{
	TiXmlDocument doc(filename.c_str());
	doc.LoadFile();
	if (doc.Error() && doc.ErrorId() == TiXmlBase::TIXML_ERROR_OPENING_FILE) {
		std::cout << "Error: can't read config file !" << std::endl;
		return false;
	}
	TiXmlElement* headerItem = doc.FirstChildElement();
	if (!headerItem)
		return false;
	std::string checkItemName = headerItem->Value();
	if (checkItemName != std::string("SimFramework"))
		return false;

	// readSimulator
	if (!readSimulatorConfigFile(headerItem)) return false;

	// read OtherInfo
	TiXmlElement* subItem = headerItem->FirstChildElement();

	while (subItem)
	{
		std::string subItemName = subItem->Value();
		if (subItemName == std::string("simulator"))
		{
			if (!readSimulatorFromConfigFile(getCenteralScene()->getCurrentSimulator(), filename, subItem))
			{
				std::cout << "Error: can't read config file !!!" << std::endl;
				return false;
			}
			_bottom->handleInitSimulator();

		}
		else if (subItemName == std::string("scene"))
		{		
			getCenteralScene()->setFromConfigFile(filename, subItem);
		}
		subItem = subItem->NextSiblingElement();
	}

	return true;
}

bool BaseMainWidget::readSimulatorConfigFile(TiXmlElement* subItem)
{
	TiXmlAttribute * attri = subItem->FirstAttribute();
	std::string simName;
	std::string runPlaformName;
	while (attri)
	{
		std::string attriName = attri->Name();
		if (attriName == std::string("name"))
		{
			simName = attri->Value();
		}
		else if (attriName == std::string("run"))
		{
			runPlaformName = attri->Value();
		}
		attri = attri->Next();
	}

	RunPlatform run = RunPlatform::CPU;
	if (runPlaformName == "OPEMMP")
		run = RunPlatform::OPEMMP;
	else if (runPlaformName == "CUDA")
		run = RunPlatform::CUDA;

	getCenteralScene()->bindSimulator(createSimulator(simName, run));
	return true;
}
