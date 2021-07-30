#include "BaseRightWidget.h"

BaseRightWidget::BaseRightWidget():_selectModel(-1)
{
	setupUi(this);
	colorPanel = new ColorSelectionPanel();
	materialPanel = new BaseMaterialEditWidget();
}

BaseRightWidget::~BaseRightWidget()
{
}

void BaseRightWidget::resizeEvent(QResizeEvent * size)
{
	int h = height();
	int w = width();
	tabWidget->resize(QSize(w,h));
}

void BaseRightWidget::setBaseScene(BaseScene * scene)
{
	_scene = scene;
	if (_scene == nullptr)
		return;

	connect(_scene, SIGNAL(initSceneSignal()), this, SLOT(bindSceneTab()));
	connect(_scene, SIGNAL(initSceneSignal()), this, SLOT(bindModelTab()));
}

void BaseRightWidget::bindSceneTab()
{
	bindCameraGroup();
	bindLightingGroup();
	bindFloorGroup();
	refreshPhysicalEnvironmentGroup();
	connect(this, SIGNAL(openColorPanel(QColor*)), this, SLOT(getColorFromPanel(QColor*)));
	connect(_scene, SIGNAL(changeSimulator()), this, SLOT(refreshPhysicalEnvironmentGroup()));
}

void BaseRightWidget::bindCameraGroup()
{	
	ortho_mode_checkBox->disconnect(SIGNAL(stateChanged(int)));
	ortho_mode_checkBox->setChecked(((StandardCamera*)_scene->camera())->queryOrthoMode());
	connect(ortho_mode_checkBox, SIGNAL(stateChanged(int)), this, SLOT(useCameraOrthoMode()));

	lock_camera_checkBox->disconnect(SIGNAL(stateChanged(int)));
	lock_camera_checkBox->setChecked(((StandardCamera*)_scene->camera())->queryLockCamera());
	connect(lock_camera_checkBox, SIGNAL(stateChanged(int)), this, SLOT(lockCamera()));
	
	show_FOV_checkBox->disconnect(SIGNAL(stateChanged(int)));
	show_FOV_checkBox->setChecked(((StandardCamera*)_scene->camera())->queryShowFOV());
	connect(show_FOV_checkBox, SIGNAL(stateChanged(int)), this, SLOT(showCameraFOV()));
}

void BaseRightWidget::bindLightingGroup()
{
	QVector3D lightPos = _scene->getLightSetting()->getMainLightPosition();
	light_pos_x->disconnect(SIGNAL(textChanged(const QString &)));
	light_pos_y->disconnect(SIGNAL(textChanged(const QString &)));
	light_pos_z->disconnect(SIGNAL(textChanged(const QString &)));
	light_pos_x->setText(QString::number(lightPos[0]));
	light_pos_y->setText(QString::number(lightPos[1]));
	light_pos_z->setText(QString::number(lightPos[2]));
	connect(light_pos_x, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightPos()));
	connect(light_pos_y, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightPos()));
	connect(light_pos_z, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightPos()));

	light_attenuation_constant->setText(QString::number(_scene->getLightSetting()->getPointLight()->constant));
	light_attenuation_linear->setText(QString::number(_scene->getLightSetting()->getPointLight()->linear));
	light_attenuation_quadratic->setText(QString::number(_scene->getLightSetting()->getPointLight()->quadratic));
	connect(light_attenuation_constant, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightAttenuation()));
	connect(light_attenuation_linear, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightAttenuation()));
	connect(light_attenuation_quadratic, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightAttenuation()));

	connect(light_ambient_color_button, SIGNAL(clicked(bool)), this, SLOT(pressMainLightAmbientBtn()));
	connect(light_diffuse_color_button, SIGNAL(clicked(bool)), this, SLOT(pressMainLightDiffuseBtn()));
	connect(light_specular_color_button, SIGNAL(clicked(bool)), this, SLOT(pressMainLightSpecularBtn()));

	setMainLightAmbientEditLine();
	light_amibent_color_r->disconnect(SIGNAL(textChanged(const QString &)));
	light_amibent_color_g->disconnect(SIGNAL(textChanged(const QString &)));
	light_amibent_color_b->disconnect(SIGNAL(textChanged(const QString &)));
	connect(light_amibent_color_r, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightAmbientRGB()));
	connect(light_amibent_color_g, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightAmbientRGB()));
	connect(light_amibent_color_b, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightAmbientRGB()));

	setMainLightDiffuseEditLine();
	light_diffuse_color_r->disconnect(SIGNAL(textChanged(const QString &)));
	light_diffuse_color_g->disconnect(SIGNAL(textChanged(const QString &)));
	light_diffuse_color_b->disconnect(SIGNAL(textChanged(const QString &)));
	connect(light_diffuse_color_r, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightDiffuseRGB()));
	connect(light_diffuse_color_g, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightDiffuseRGB()));
	connect(light_diffuse_color_b, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightDiffuseRGB()));

	setMainLightSpecularEditLine();
	light_specular_color_r->disconnect(SIGNAL(textChanged(const QString &)));
	light_specular_color_g->disconnect(SIGNAL(textChanged(const QString &)));
	light_specular_color_b->disconnect(SIGNAL(textChanged(const QString &)));
	connect(light_specular_color_r, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightSpecularRGB()));
	connect(light_specular_color_g, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightSpecularRGB()));
	connect(light_specular_color_b, SIGNAL(textChanged(const QString &)), this, SLOT(setMainLightSpecularRGB()));

	shadow_checkBox->setChecked(_scene->getLightSetting()->useShadow());
	shadow_checkBox->disconnect(SIGNAL(stateChanged(int)));
	connect(shadow_checkBox, SIGNAL(stateChanged(int)), this, SLOT(enableShadow()));

}

void BaseRightWidget::bindFloorGroup()
{
	connect(floor_metrial_button, SIGNAL(clicked(bool)), this, SLOT(handleFloorMaterialBtn()));
	
	floor_hide_checkBox->disconnect(SIGNAL(stateChanged(int)));
	floor_hide_checkBox->setChecked(_scene->getFloor()->isHide());
	connect(floor_hide_checkBox, SIGNAL(stateChanged(int)), this, SLOT(hideFloor()));
	
	floor_xz_scale->disconnect(SIGNAL(textChanged(const QString &)));
	floor_xz_scale->setText(QString::number(_scene->getFloor()->getScaleXZ()));
	connect(floor_xz_scale, SIGNAL(textChanged(const QString &)), this, SLOT(scaleFloorXZ()));

	floor_y_plane->disconnect(SIGNAL(textChanged(const QString &)));
	floor_y_plane->setText(QString::number(_scene->getFloor()->getYPlane()));
	connect(floor_y_plane, SIGNAL(textChanged(const QString &)), this, SLOT(setFloorY()));
}

void BaseRightWidget::bindModelTab()
{
	bindModelList();
	bindMeshGroup();
}

void BaseRightWidget::bindModelList()
{
	connect(mesh_listWidget, SIGNAL(currentRowChanged(int)), this, SLOT(selectModel(int)));
	connect(mesh_listWidget, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(clickModel(QListWidgetItem*)));
}

void BaseRightWidget::bindMeshGroup()
{
	refreshMeshGroup();
	refreshVolumetricMeshGroup();
	refreshMedialMeshGroup();
}

/* 
slot
*/
void BaseRightWidget::getColorFromPanel(QColor* color)
{
	//
	colorPanel->bindColor(color);
	colorPanel->show();
}

// Scene Tab: Camera group
void BaseRightWidget::useCameraOrthoMode()
{
	((StandardCamera*)_scene->camera())->enableOrthoMode(ortho_mode_checkBox->isChecked());
}

void BaseRightWidget::lockCamera()
{
	((StandardCamera*)_scene->camera())->enableLockCamera(lock_camera_checkBox->isChecked());
}

void BaseRightWidget::showCameraFOV()
{
	((StandardCamera*)_scene->camera())->enableShowCameraFOV(show_FOV_checkBox->isChecked());
}

// Scene Tab: Lighting group
void BaseRightWidget::setMainLightPos()
{
	qeal x = light_pos_x->text().toDouble();
	qeal y = light_pos_y->text().toDouble();
	qeal z = light_pos_z->text().toDouble();

	_scene->getLightSetting()->getPointLight()->lightPos = QVector3D(x, y, z);
}

void BaseRightWidget::setMainLightAttenuation()
{
	_scene->getLightSetting()->getPointLight()->constant = light_attenuation_constant->text().toDouble();
	_scene->getLightSetting()->getPointLight()->linear = light_attenuation_linear->text().toDouble();
	_scene->getLightSetting()->getPointLight()->quadratic = light_attenuation_quadratic->text().toDouble();
}

void BaseRightWidget::pressMainLightAmbientBtn()
{
	colorPanel->disconnect(SIGNAL(panelClose()));
	connect(colorPanel, SIGNAL(panelClose()), this, SLOT(setMainLightAmbientEditLine()));
	emit openColorPanel(&(_scene->getLightSetting()->getPointLight()->ambientColor));
}

void BaseRightWidget::pressMainLightDiffuseBtn()
{
	colorPanel->disconnect(SIGNAL(panelClose()));
	connect(colorPanel, SIGNAL(panelClose()), this, SLOT(setMainLightDiffuseEditLine()));
	emit openColorPanel(&(_scene->getLightSetting()->getPointLight()->diffuseColor));
}

void BaseRightWidget::pressMainLightSpecularBtn()
{
	colorPanel->disconnect(SIGNAL(panelClose()));
	connect(colorPanel, SIGNAL(panelClose()), this, SLOT(setMainLightSpecularEditLine()));
	emit openColorPanel(&(_scene->getLightSetting()->getPointLight()->specularColor));
}

void BaseRightWidget::setMainLightAmbientRGB()
{
	int r = light_amibent_color_r->text().toInt();
	int g = light_amibent_color_g->text().toInt();
	int b = light_amibent_color_b->text().toInt();
	_scene->getLightSetting()->getPointLight()->ambientColor = QColor(r, g, b);
}

void BaseRightWidget::setMainLightAmbientEditLine()
{
	QColor color = _scene->getLightSetting()->getPointLight()->ambientColor;
	int r = color.red();
	int g = color.green();
	int b = color.blue();
	light_amibent_color_r->setText(QString::number(r));
	light_amibent_color_g->setText(QString::number(g));
	light_amibent_color_b->setText(QString::number(b));
}

void BaseRightWidget::setMainLightDiffuseRGB()
{
	int r = light_diffuse_color_r->text().toInt();
	int g = light_diffuse_color_g->text().toInt();
	int b = light_diffuse_color_b->text().toInt();
	_scene->getLightSetting()->getPointLight()->diffuseColor = QColor(r, g, b);
}

void BaseRightWidget::setMainLightDiffuseEditLine()
{
	QColor color = _scene->getLightSetting()->getPointLight()->diffuseColor;
	int r = color.red();
	int g = color.green();
	int b = color.blue();
	light_diffuse_color_r->setText(QString::number(r));
	light_diffuse_color_g->setText(QString::number(g));
	light_diffuse_color_b->setText(QString::number(b));
}

void BaseRightWidget::setMainLightSpecularRGB()
{
	int r = light_specular_color_r->text().toInt();
	int g = light_specular_color_g->text().toInt();
	int b = light_specular_color_b->text().toInt();
	_scene->getLightSetting()->getPointLight()->specularColor = QColor(r, g, b);
}

void BaseRightWidget::setMainLightSpecularEditLine()
{
	QColor color = _scene->getLightSetting()->getPointLight()->specularColor;
	int r = color.red();
	int g = color.green();
	int b = color.blue();
	light_specular_color_r->setText(QString::number(r));
	light_specular_color_g->setText(QString::number(g));
	light_specular_color_b->setText(QString::number(b));
}

void BaseRightWidget::enableShadow()
{
	_scene->getLightSetting()->enableShadow(shadow_checkBox->isChecked());
}

// Scene Tab: Floor group
void BaseRightWidget::handleFloorMaterialBtn()
{
	materialPanel->openWidget(_scene->getFloor()->renderMaterials[0], _scene->getFloor()->renderMaterials.size());
}

void BaseRightWidget::hideFloor()
{
	_scene->getFloor()->enableHide(floor_hide_checkBox->isChecked());
}

void BaseRightWidget::scaleFloorXZ()
{
	_scene->getFloor()->scaleXZ(floor_xz_scale->text().toInt());
}

void BaseRightWidget::setFloorY()
{
	_scene->getFloor()->setYPlane(floor_y_plane->text().toDouble());
}

// Scene Tab: Physical Environment
void BaseRightWidget::refreshPhysicalEnvironmentGroup()
{
	BaseSimulator* sim = _scene->getCurrentSimulator();
	bool valid = true;
	if (sim == nullptr)
		valid = false;
	gravity_lineEdit->setEnabled(valid);
	gravity_checkBox->setEnabled(valid);
	add_extra_force_button->setEnabled(valid);
	time_step_lineEdit->setEnabled(valid);
	if (!valid)
		return;

	gravity_lineEdit->disconnect(SIGNAL(textChanged(const QString &)));
	gravity_lineEdit->setText(QString::number(sim->getGraviry()));
	connect(gravity_lineEdit, SIGNAL(textChanged(const QString &)), this, SLOT(setGrivaty()));

	gravity_checkBox->disconnect(SIGNAL(stateChanged(int)));
	gravity_checkBox->setChecked(sim->isUseGravity());
	connect(gravity_checkBox, SIGNAL(stateChanged(int)), this, SLOT(enableGrivaty()));

	time_step_lineEdit->disconnect(SIGNAL(textChanged(const QString &)));
	time_step_lineEdit->setText(QString::number(sim->getTimeStep()));
	connect(time_step_lineEdit, SIGNAL(textChanged(const QString &)), this, SLOT(setTimeStep()));

	connect(add_extra_force_button, SIGNAL(clicked(bool)), this, SLOT(handleAddExtraForceBtn()));
}

void BaseRightWidget::enableGrivaty()
{
	BaseSimulator* sim = _scene->getCurrentSimulator();
	if (sim == nullptr)
		return;
	sim->enableGravity(gravity_checkBox->isChecked());
}

void BaseRightWidget::setGrivaty()
{
	BaseSimulator* sim = _scene->getCurrentSimulator();
	if (sim == nullptr)
		return;
	sim->setGraviry(gravity_lineEdit->text().toDouble());
	emit needInitSimulator();
}

void BaseRightWidget::setTimeStep()
{
	BaseSimulator* sim = _scene->getCurrentSimulator();
	if (sim == nullptr)
		return;
	sim->setTimeStep(time_step_lineEdit->text().toDouble());
	emit needInitSimulator();
}

void BaseRightWidget::handleAddExtraForceBtn()
{
	//to do
	QString filename = QFileDialog::getOpenFileName(this, tr("Select a Surface to import"), "./example/", tr(""));
	if (filename.isEmpty())
		return;
	debugModifiledFile(filename.toStdString());
}

// Model Tab: Model List
void BaseRightWidget::refreshModelList()
{	
	mesh_listWidget->clear();
	if (_scene == nullptr)
		return;

	BaseSimulator* sim = _scene->getCurrentSimulator();
	if (sim == nullptr)
		return;

	for (int i = 0; i < sim->models.size(); i++)
	{
		mesh_listWidget->addItem(sim->models[i]->nickName.c_str());
		if (i == 0)
			mesh_listWidget->setCurrentRow(0);
	}

	for (int i = 0; i < sim->staticModels.size(); i++)
	{
		mesh_listWidget->addItem(sim->staticModels[i]->nickName.c_str());
		if (i == 0 && sim->models.size() == 0)
			mesh_listWidget->setCurrentRow(0);
	}
}

void BaseRightWidget::selectModel(int index)
{
	_selectModel = index;
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;

	refreshMeshGroup();
	refreshVolumetricMeshGroup();
	refreshMedialMeshGroup();
}

void BaseRightWidget::clickModel(QListWidgetItem * item)
{
	if (_selectModel < 0 || _selectModel >= _scene->getCurrentSimulator()->models.size())
		return;
	BaseModel* m = _scene->getCurrentSimulator()->models[_selectModel];

	m->clickModelList();
}

void BaseRightWidget::refreshMeshGroup()
{
	bool valid = true;
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		valid = false;

	BaseModel* m = nullptr;
	if(valid)
	{
		if (_selectModel >= _scene->getCurrentSimulator()->models.size())
		{
			int id = _selectModel - _scene->getCurrentSimulator()->models.size();
			m = _scene->getCurrentSimulator()->staticModels[id];
		}
		else m = _scene->getCurrentSimulator()->models[_selectModel];
		if (m == nullptr)
			valid = false;
	}
	
	mesh_material_button->setEnabled(valid);
	mesh_hide_checkBox->setEnabled(valid);
	mesh_translation_x->setEnabled(valid);
	mesh_translation_y->setEnabled(valid);
	mesh_translation_z->setEnabled(valid);
	mesh_rotation_x->setEnabled(valid);
	mesh_rotation_y->setEnabled(valid);
	mesh_rotation_z->setEnabled(valid);
	mesh_rotation_sita->setEnabled(valid);
	mesh_scale_lineEdit->setEnabled(valid);

	mesh_material_button->disconnect(SIGNAL(clicked(bool)));
	mesh_hide_checkBox->disconnect(SIGNAL(stateChanged(int)));
	mesh_translation_x->disconnect(SIGNAL(textChanged(const QString &)));
	mesh_translation_y->disconnect(SIGNAL(textChanged(const QString &)));
	mesh_translation_z->disconnect(SIGNAL(textChanged(const QString &)));
	mesh_rotation_x->disconnect(SIGNAL(textChanged(const QString &)));
	mesh_rotation_y->disconnect(SIGNAL(textChanged(const QString &)));
	mesh_rotation_z->disconnect(SIGNAL(textChanged(const QString &)));
	mesh_rotation_sita->disconnect(SIGNAL(textChanged(const QString &)));
	mesh_scale_lineEdit->disconnect(SIGNAL(textChanged(const QString &)));
	if (!valid)
		return;

	mesh_hide_checkBox->setChecked(m->getSurfaceHandle()->getMesh()->isHide());
	qeal tx, ty, tz, rx, ry, rz, rsita, scale;
	m->getTranslation(tx, ty, tz);
	m->getRotation(rx, ry, rz, rsita);
	m->getScale(scale);

	mesh_translation_x->setText(QString::number(tx));
	mesh_translation_y->setText(QString::number(ty));
	mesh_translation_z->setText(QString::number(tz));
	mesh_rotation_x->setText(QString::number(rx));
	mesh_rotation_y->setText(QString::number(ry));
	mesh_rotation_z->setText(QString::number(rz));
	mesh_rotation_sita->setText(QString::number(rsita));
	mesh_scale_lineEdit->setText(QString::number(scale));

	connect(mesh_material_button, SIGNAL(clicked(bool)), this, SLOT(handleMeshMaterialBtn()));
	connect(mesh_hide_checkBox, SIGNAL(stateChanged(int)), this, SLOT(hideMesh()));
	connect(mesh_translation_x, SIGNAL(textChanged(const QString &)), this, SLOT(setMeshTranslation()));
	connect(mesh_translation_y, SIGNAL(textChanged(const QString &)), this, SLOT(setMeshTranslation()));
	connect(mesh_translation_z, SIGNAL(textChanged(const QString &)), this, SLOT(setMeshTranslation()));
	connect(mesh_rotation_x, SIGNAL(textChanged(const QString &)), this, SLOT(setMeshRotation()));
	connect(mesh_rotation_y, SIGNAL(textChanged(const QString &)), this, SLOT(setMeshRotation()));
	connect(mesh_rotation_z, SIGNAL(textChanged(const QString &)), this, SLOT(setMeshRotation()));
	connect(mesh_rotation_sita, SIGNAL(textChanged(const QString &)), this, SLOT(setMeshRotation()));
	connect(mesh_scale_lineEdit, SIGNAL(textChanged(const QString &)), this, SLOT(setMeshScale()));
}

void BaseRightWidget::handleMeshMaterialBtn()
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;

	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr) return;
	materialPanel->openWidget(m->renderMaterials[0], m->renderMaterials.size());
}

void BaseRightWidget::hideMesh()
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;

	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr) return;

	m->getSurfaceHandle()->getMesh()->enableHide(mesh_hide_checkBox->isChecked());
}

void BaseRightWidget::setMeshScale()
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;
	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr) return;
	m->scaleModel(mesh_scale_lineEdit->text().toDouble());
	emit needInitSimulator();
}

void BaseRightWidget::setMeshTranslation()
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;
	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr) return;
	qeal tx = mesh_translation_x->text().toDouble();
	qeal ty = mesh_translation_y->text().toDouble();
	qeal tz = mesh_translation_z->text().toDouble();
	m->translateModel(tx, ty, tz);
	emit needInitSimulator();
}

void BaseRightWidget::setMeshRotation()
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;
	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr) return;

	qeal rx = mesh_rotation_x->text().toDouble();
	qeal ry = mesh_rotation_y->text().toDouble();
	qeal rz = mesh_rotation_z->text().toDouble();
	qeal rsita = mesh_rotation_sita->text().toDouble();
	m->rotateModel(rx, ry, rz, rsita);
	emit needInitSimulator();
}

void BaseRightWidget::refreshVolumetricMeshGroup()
{
	element_set_listWidget->clear();
	element_set_listWidget->disconnect(SIGNAL(currentRowChanged(int)));
	element_set_listWidget->disconnect(SIGNAL(itemDoubleClicked(QListWidgetItem*)));

	bool valid = true;
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;

	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr)
	{
		tet_mesh_hide_checkBox->setEnabled(false);
		tet_mesh_transparent_checkBox->setEnabled(false);

		tet_mesh_hide_checkBox->disconnect(SIGNAL(stateChanged(int)));
		tet_mesh_transparent_checkBox->disconnect(SIGNAL(stateChanged(int)));
		return;
	}

	if (!m->isTetMeshValid()) valid = false;

	tet_mesh_hide_checkBox->setEnabled(valid);
	tet_mesh_transparent_checkBox->setEnabled(valid);

	tet_mesh_hide_checkBox->disconnect(SIGNAL(stateChanged(int)));
	tet_mesh_transparent_checkBox->disconnect(SIGNAL(stateChanged(int)));

	if (!valid)
		return;

	tet_mesh_hide_checkBox->setChecked(m->getTetMeshHandle()->getTetMesh()->isHide());
	tet_mesh_transparent_checkBox->setChecked(m->getTetMeshHandle()->getTetMesh()->isTransparent());

	connect(tet_mesh_hide_checkBox, SIGNAL(stateChanged(int)), this, SLOT(hideTetMesh()));
	connect(tet_mesh_transparent_checkBox, SIGNAL(stateChanged(int)), this, SLOT(setTetMeshTransparent()));

	//
	connect(element_set_listWidget, SIGNAL(currentRowChanged(int)), this, SLOT(selectTetElementSet(int)));
	connect(element_set_listWidget, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(clickTetElementSetList(QListWidgetItem*)));
	std::vector<BaseTetElementSet> elementSet = m->getTetMeshHandle()->getTetMeshElementSet();
	element_set_listWidget->addItem("all_elements");
	element_set_listWidget->setCurrentRow(0);

	if (elementSet.size() > 1)
	{
		for (int i = 0; i < elementSet.size(); i++)
		{
			element_set_listWidget->addItem(elementSet[i].getName().c_str());
		}
	}	
}

void BaseRightWidget::hideTetMesh()
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;

	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr)
		return;
	if (!m->getTetMeshHandle()->getTetMesh()->isTetMeshValid())
		return;
	m->getTetMeshHandle()->getTetMesh()->enableHide(tet_mesh_hide_checkBox->isChecked());
}

void BaseRightWidget::setTetMeshTransparent()
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;

	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr)
		return;
	if (!m->getTetMeshHandle()->getTetMesh()->isTetMeshValid())
		return;
	m->getTetMeshHandle()->getTetMesh()->enableTransparent(tet_mesh_transparent_checkBox->isChecked());
}

void BaseRightWidget::selectTetElementSet(int id)
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;

	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (!m || !m->isTetMeshValid()) return;

	m->getTetMeshHandle()->setRenderElementSetId(id);
}

void BaseRightWidget::clickTetElementSetList(QListWidgetItem * item)
{
}

void BaseRightWidget::refreshMedialMeshGroup()
{
	medial_mesh_set_listWidget->clear();
	medial_mesh_set_listWidget->disconnect(SIGNAL(currentRowChanged(int)));
	medial_mesh_set_listWidget->disconnect(SIGNAL(itemDoubleClicked(QListWidgetItem*)));
	bool valid = true;
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;

	BaseModel* m = nullptr;

	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr)
	{
		medial_mesh_hide_checkBox->setEnabled(false);
		medial_mesh_transparent_checkBox->setEnabled(false);

		medial_mesh_hide_checkBox->disconnect(SIGNAL(stateChanged(int)));
		medial_mesh_transparent_checkBox->disconnect(SIGNAL(stateChanged(int)));
		return;
	}

	if(!m->isMedialMeshValid()) valid = false;

	medial_mesh_hide_checkBox->setEnabled(valid);
	medial_mesh_transparent_checkBox->setEnabled(valid);

	medial_mesh_hide_checkBox->disconnect(SIGNAL(stateChanged(int)));
	medial_mesh_transparent_checkBox->disconnect(SIGNAL(stateChanged(int)));

	if (!valid)
		return;

	medial_mesh_hide_checkBox->setChecked(m->getMedialMeshHandle()->getMesh()->isHide());
	medial_mesh_transparent_checkBox->setChecked(m->getMedialMeshHandle()->getMesh()->isTransparent());

	connect(medial_mesh_hide_checkBox, SIGNAL(stateChanged(int)), this, SLOT(hideMedialMesh()));
	connect(medial_mesh_transparent_checkBox, SIGNAL(stateChanged(int)), this, SLOT(setMedialMeshTransparent()));

	//
	connect(medial_mesh_set_listWidget, SIGNAL(currentRowChanged(int)), this, SLOT(selectMedialMeshSet(int)));
	connect(medial_mesh_set_listWidget, SIGNAL(itemDoubleClicked(QListWidgetItem*)), this, SLOT(clickMedialMeshSetList(QListWidgetItem*)));
	std::vector<BaseMedialPrimitiveSet> _primitiveSet;
	std::vector<BaseMedialPointSet> _pointSet;
	m->getMedialMeshHandle()->getMedialMeshSet(_primitiveSet, _pointSet);
	medial_mesh_set_listWidget->addItem("all_Primitives");
	medial_mesh_set_listWidget->setCurrentRow(0);
	medial_mesh_set_listWidget->addItem("all_Spheres");

	if (_primitiveSet.size() > 1)
	{
		for (int i = 0; i < _primitiveSet.size(); i++)
		{
			medial_mesh_set_listWidget->addItem(_primitiveSet[i].getName().c_str());
			medial_mesh_set_listWidget->addItem(_pointSet[i].getName().c_str());
		}
	}
}

void BaseRightWidget::hideMedialMesh()
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;

	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr)
		return;
	if (!m->getMedialMeshHandle()->getMesh()->isMedialMeshValid())
		return;
	m->getMedialMeshHandle()->getMesh()->enableHide(medial_mesh_hide_checkBox->isChecked());
}

void BaseRightWidget::setMedialMeshTransparent()
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;

	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (m == nullptr)
		return;
	if (!m->getMedialMeshHandle()->getMesh()->isMedialMeshValid())
		return;
	m->getMedialMeshHandle()->getMesh()->enableTransparent(medial_mesh_transparent_checkBox->isChecked());
}

void BaseRightWidget::selectMedialMeshSet(int id)
{
	if (_selectModel < 0 || _selectModel >= (_scene->getCurrentSimulator()->models.size() + _scene->getCurrentSimulator()->staticModels.size()))
		return;
	BaseModel* m = nullptr;

	if (_selectModel >= _scene->getCurrentSimulator()->models.size())
	{
		int id = _selectModel - _scene->getCurrentSimulator()->models.size();
		m = _scene->getCurrentSimulator()->staticModels[id];
	}
	else m = _scene->getCurrentSimulator()->models[_selectModel];

	if (!m || !m->isMedialMeshValid()) return;

	m->getMedialMeshHandle()->setRenderMedialMeshSetId(id);
}

void BaseRightWidget::clickMedialMeshSetList(QListWidgetItem * item)
{
}

void BaseRightWidget::debugModifiledFile(const std::string filename)
{
	std::ifstream fin(filename.c_str());
	int vn, cn, sn;
	std::vector<qeal> vlist;
	std::vector<int> vtlist;
	std::vector<int> clist;
	std::vector<int> slist;

	qeal sx, miny, maxy, minz, maxz;
	qeal rx;
	miny = DBL_MAX;
	maxy = -DBL_MAX;
	minz = DBL_MAX;
	maxz = -DBL_MAX;

	fin >> vn >> cn >> sn;
	for (int i = 0; i < vn; i++)
	{
		qeal x, y, z, r;
		int d;
		char ch;
		fin >> ch >> x >> y >> z >> r >> d;
		vlist.push_back(x);
		vlist.push_back(y);
		vlist.push_back(z);
		vlist.push_back(r);
		vtlist.push_back(d);

		if (y < miny) miny = y;
		if (y > maxy) maxy = y;

		if (z < minz) minz = z;
		if (z > maxz) maxz = z;
	}
	sx = vlist[0];
	rx = vlist[3];
	for (int i = 0; i < cn; i++)
	{
		int v0, v1;
		char ch;
		fin >> ch >> v0 >> v1;
		clist.push_back(v0);
		clist.push_back(v1);
	}
	for (int i = 0; i < sn; i++)
	{
		int v0, v1, v2;
		char ch;
		fin >> ch >> v0 >> v1 >> v2;
		slist.push_back(v0);
		slist.push_back(v1);
		slist.push_back(v2);
	}
	fin.close();

	Vector3 s0(sx, miny, minz);
	Vector3 s1(sx, maxy, minz);
	Vector3 s2(sx, maxy, maxz);
	Vector3 s3(sx, miny, maxz);


	int step = 10;
	qeal yStep = (maxy - miny) / (step - 1);
	qeal zStep = (maxz - minz) / (step - 1);

	std::vector<std::vector<Vector3>> nslist(step);
	std::vector<std::vector<int>> nsIndexlist(step);
	std::vector<std::vector<int>> nsIndexTdlist(step);
	int index = 0;
	int tds = 9346;
	for (int r = 0; r < step; r++)
	{
		for (int c = 0; c < step; c++)
		{
			Vector3 s(sx, miny + c * yStep, minz + r * zStep);
			nslist[r].push_back(s);
			nsIndexlist[r].push_back(index);
			nsIndexTdlist[r].push_back(tds);
			index++;
			tds++;
		}
	}

	std::vector<Vector3i> slabList;

	for (int i = 0; i < nsIndexlist.size(); i++)
	{
		for (int j = 0; j < nsIndexlist[i].size(); j++)
		{
			std::cout << nsIndexlist[i][j] << " ";
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}

	for (int i = 1; i < step; i++)
	{
		for (int j = i - 1; j < i; j++)
		{
			for (int ci = 0; ci < step - 1; ci++)
			{
				int vi0, vi1, vj0, vj1;
				vi0 = nsIndexlist[i][ci];
				vi1 = nsIndexlist[i][ci + 1];
				vj0 = nsIndexlist[j][ci];
				vj1 = nsIndexlist[j][ci + 1];

				Vector3i s0 = Vector3i(vj0, vj1, vi0);
				Vector3i s1 = Vector3i(vj1, vi1, vi0);
				slabList.push_back(s0);
				slabList.push_back(s1);
			}
		}
	}

	std::string fi = filename;
	fi += "t";

	std::ofstream fout(fi.c_str());
	fout << index << " " << 0 << " " << slabList.size() << std::endl;
	for (int i = 0; i < nslist.size(); i++)
	{
		for (int j = 0; j < nslist[i].size(); j++)
		{
			fout << "v " << nslist[i][j].data()[0] << " " << nslist[i][j].data()[1] << " " << nslist[i][j].data()[2] << " " << rx << " " << nsIndexTdlist[i][j] << std::endl;
		}
	}
	for (int i = 0; i < slabList.size(); i++)
	{
		fout << "s " << slabList[i].data()[0] << " " << slabList[i].data()[1] << " " << slabList[i].data()[2] << std::endl;
	}
	fout.close();

	fi = filename;
	fi += ".a.node";
	fout.open(fi.c_str());
	fout << index << " 3 0 0"<< std::endl;
	index = 0;
	for (int i = 0; i < nslist.size(); i++)
	{
		for (int j = 0; j < nslist[i].size(); j++)
		{
			fout << index  <<" " << nslist[i][j].data()[0] << " " << nslist[i][j].data()[1] << " " << nslist[i][j].data()[2] << std::endl;
			index++;
		}
	}
	fout.close();
}