#include "BaseMaterialEditWidget.h"

BaseMaterialEditWidget::BaseMaterialEditWidget():_materials(nullptr), _materialNum(0)
{
	setupUi(this);
	setFixedSize(QSize(495, 320));
	setWindowTitle(QString("MaterialEditor"));
	ShinnessSlider->setMaximum(1280);
	ShinnessSlider->setMinimum(0);

	colorPanel = new ColorSelectionPanel();

	connect(this, SIGNAL(openColorPanel(QColor*)), this, SLOT(handleColorPanel(QColor*)));

	connect(material_listWidget, SIGNAL(currentRowChanged(int)), this, SLOT(selectMaterial(int)));

	connect(amibent_color_r, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromAmbientEditLine()));
	connect(amibent_color_g, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromAmbientEditLine()));
	connect(amibent_color_b, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromAmbientEditLine()));

	connect(diffuse_color_r, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromDiffuseEditLine()));
	connect(diffuse_color_g, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromDiffuseEditLine()));
	connect(diffuse_color_b, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromDiffuseEditLine()));

	connect(specular_color_r, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromSpecularEditLine()));
	connect(specular_color_g, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromSpecularEditLine()));
	connect(specular_color_b, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromSpecularEditLine()));

	connect(ambient_color_button, SIGNAL(clicked(bool)), this, SLOT(handleAmbientBtn()));
	connect(diffuse_color_button, SIGNAL(clicked(bool)), this, SLOT(handleDiffuseBtn()));
	connect(specular_color_button, SIGNAL(clicked(bool)), this, SLOT(handleSpecularBtn()));

	connect(transparent_checkBox, SIGNAL(stateChanged(int)), this, SLOT(handleTransparentCheckBox()));

	connect(shinness_lineEdit, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromShinnessEditLine()));
	connect(ShinnessSlider, SIGNAL(valueChanged(int)), this, SLOT(updateFromShinnessSlider(int)));

	ambient_texture_dir_lineEdit->setReadOnly(true);
	diffuse_texture_dir_lineEdit->setReadOnly(true);
	specular_texture_dir_lineEdit->setReadOnly(true);
	bump_texture_dir_lineEdit->setReadOnly(true);

	connect(enable_amibent_texture_checkBox, SIGNAL(stateChanged(int)), this, SLOT(enableAmbientTexture(int)));
	connect(enable_diffuse_texture_checkBox, SIGNAL(stateChanged(int)), this, SLOT(enableDiffuseTexture(int)));
	connect(enable_specular_texture_checkBox, SIGNAL(stateChanged(int)), this, SLOT(enableSpecularTexture(int)));
	connect(enable_bump_texture_checkBox, SIGNAL(stateChanged(int)), this, SLOT(enableBumpTexture(int)));

	connect(ambient_texture_button, SIGNAL(clicked(bool)), this, SLOT(handleTextureAmbientBtn()));
	connect(diffuse_texture_button, SIGNAL(clicked(bool)), this, SLOT(handleTextureDiffuseBtn()));
	connect(specular_texture_button, SIGNAL(clicked(bool)), this, SLOT(handleTextureSpecularBtn()));
	connect(bump_texture_button, SIGNAL(clicked(bool)), this, SLOT(handleTextureBumpBtn()));
}

void BaseMaterialEditWidget::openWidget(BaseRenderMaterial* materials, int size)
{
	_materials = materials;
	_materialNum = size;
	show();
	if (_materials == nullptr || _materialNum == 0)
		return;	
	
	for (int i = 0; i < _materialNum; i++)
	{
		QString name = QString("matrial ") + QString::number(i);
		material_listWidget->addItem(name);
		if (i == 0)
			material_listWidget->setCurrentRow(0);
	}	
}

void BaseMaterialEditWidget::selectMaterial(int index)
{
	if (_materials == nullptr)
		return;
	_selectedMaterial = _materials + index;
	refreshTextureBar();
	refreshAmbient();
	refreshDiffuse();
	refreshSpecular();
	refreshShinness();

	enable_amibent_texture_checkBox->setChecked(_selectedMaterial->useTextureMap[0]);
	enable_diffuse_texture_checkBox->setChecked(_selectedMaterial->useTextureMap[1]);
	enable_specular_texture_checkBox->setChecked(_selectedMaterial->useTextureMap[2]);
	enable_bump_texture_checkBox->setChecked(_selectedMaterial->useTextureMap[3]);

	enableAmbientTexture(enable_amibent_texture_checkBox->isChecked());
	enableDiffuseTexture(enable_diffuse_texture_checkBox->isChecked());
	enableSpecularTexture(enable_specular_texture_checkBox->isChecked());
	enableBumpTexture(enable_bump_texture_checkBox->isChecked());
}

void BaseMaterialEditWidget::handleColorPanel(QColor *color)
{
	colorPanel->bindColor(color);
	colorPanel->show();
}

void BaseMaterialEditWidget::refreshAmbient()
{
	if (_selectedMaterial == nullptr)
		return;
	QColor color = _selectedMaterial->ambient;

	amibent_color_r->setText(QString::number(color.red()));
	amibent_color_g->setText(QString::number(color.green()));
	amibent_color_b->setText(QString::number(color.blue()));
}

void BaseMaterialEditWidget::updateFromAmbientEditLine()
{
	if (_selectedMaterial == nullptr)
		return;
	qeal red = amibent_color_r->text().toDouble();
	qeal green = amibent_color_g->text().toDouble();
	qeal blue = amibent_color_b->text().toDouble();
	_selectedMaterial->ambient = QColor(red, green, blue);
}

void BaseMaterialEditWidget::refreshDiffuse()
{
	if (_selectedMaterial == nullptr)
		return;
	QColor color = _selectedMaterial->diffuse;
	diffuse_color_r->setText(QString::number(color.red()));
	diffuse_color_g->setText(QString::number(color.green()));
	diffuse_color_b->setText(QString::number(color.blue()));
}

void BaseMaterialEditWidget::updateFromDiffuseEditLine()
{
	if (_selectedMaterial == nullptr)
		return;
	qeal red = diffuse_color_r->text().toDouble();
	qeal green = diffuse_color_g->text().toDouble();
	qeal blue = diffuse_color_b->text().toDouble();
	_selectedMaterial->diffuse = QColor(red, green, blue);
}

void BaseMaterialEditWidget::refreshSpecular()
{
	if (_selectedMaterial == nullptr)
		return;
	QColor color = _selectedMaterial->specular;
	specular_color_r->setText(QString::number(color.red()));
	specular_color_g->setText(QString::number(color.green()));
	specular_color_b->setText(QString::number(color.blue()));
}

void BaseMaterialEditWidget::updateFromSpecularEditLine()
{
	if (_selectedMaterial == nullptr)
		return;
	qeal red = specular_color_r->text().toDouble();
	qeal green = specular_color_g->text().toDouble();
	qeal blue = specular_color_b->text().toDouble();
	_selectedMaterial->specular = QColor(red, green, blue);
}

void BaseMaterialEditWidget::handleAmbientBtn()
{
	if (_selectedMaterial == nullptr)
		return;
	colorPanel->disconnect(SIGNAL(panelClose()));
	connect(colorPanel, SIGNAL(panelClose()), this, SLOT(refreshAmbient()));
	emit openColorPanel(&(_selectedMaterial->ambient));
}

void BaseMaterialEditWidget::handleDiffuseBtn()
{
	if (_selectedMaterial == nullptr)
		return;
	colorPanel->disconnect(SIGNAL(panelClose()));
	connect(colorPanel, SIGNAL(panelClose()), this, SLOT(refreshDiffuse()));
	emit openColorPanel(&(_selectedMaterial->diffuse));
}

void BaseMaterialEditWidget::handleSpecularBtn()
{
	if (_selectedMaterial == nullptr)
		return;
	colorPanel->disconnect(SIGNAL(panelClose()));
	connect(colorPanel, SIGNAL(panelClose()), this, SLOT(refreshSpecular()));
	emit openColorPanel(&(_selectedMaterial->specular));
}

void BaseMaterialEditWidget::handleTransparentCheckBox()
{
	if (_selectedMaterial == nullptr)
		return;
	_selectedMaterial->enableTransparent(transparent_checkBox->isChecked());
}

void BaseMaterialEditWidget::refreshShinness()
{
	if (_selectedMaterial == nullptr)
		return;
	qeal shinness = _selectedMaterial->shinness;
	shinness_lineEdit->setText(QString::number(shinness));
	shinness *= 10;
	ShinnessSlider->setValue(shinness);
}

void BaseMaterialEditWidget::updateFromShinnessEditLine()
{
	if (_selectedMaterial == nullptr)
		return;
	qeal shinness = shinness_lineEdit->text().toDouble();
	ShinnessSlider->setValue(shinness * 10);
	_selectedMaterial->shinness = shinness;
}

void BaseMaterialEditWidget::updateFromShinnessSlider(int value)
{
	if (_selectedMaterial == nullptr)
		return;
	qeal shinness = qeal(value) / 10.0;
	shinness_lineEdit->setText(QString::number(shinness));
	_selectedMaterial->shinness = shinness;
}

void BaseMaterialEditWidget::refreshTextureBar()
{
	if (_selectedMaterial == nullptr)
		return;
	ambient_texture_dir_lineEdit->setText(QString(_selectedMaterial->textureMapFilename[0].c_str()));
	diffuse_texture_dir_lineEdit->setText(QString(_selectedMaterial->textureMapFilename[1].c_str()));
	specular_texture_dir_lineEdit->setText(QString(_selectedMaterial->textureMapFilename[2].c_str()));
	bump_texture_dir_lineEdit->setText(QString(_selectedMaterial->textureMapFilename[3].c_str()));
}

void BaseMaterialEditWidget::enableAmbientTexture(int flag)
{
	ambient_texture_button->setEnabled(flag);
	ambient_texture_dir_lineEdit->setEnabled(flag);
	ambient_color_button->setEnabled(!flag);
	amibent_color_r->setEnabled(!flag);
	amibent_color_g->setEnabled(!flag);
	amibent_color_b->setEnabled(!flag);

	if (!flag)
	{
		updateFromAmbientEditLine();
		_selectedMaterial->useTextureMap[AmbientMapIndex] = 0;
	}
	else
		_selectedMaterial->useTextureMap[AmbientMapIndex] = 1;
}

void BaseMaterialEditWidget::enableDiffuseTexture(int flag)
{
	diffuse_texture_button->setEnabled(flag);
	diffuse_texture_dir_lineEdit->setEnabled(flag);
	diffuse_color_button->setEnabled(!flag);
	diffuse_color_r->setEnabled(!flag);
	diffuse_color_g->setEnabled(!flag);
	diffuse_color_b->setEnabled(!flag);

	if (!flag)
	{
		updateFromDiffuseEditLine();
		_selectedMaterial->useTextureMap[DiffuseMapIndex] = 0;
	}
	else
		_selectedMaterial->useTextureMap[DiffuseMapIndex] = 1;
}

void BaseMaterialEditWidget::enableSpecularTexture(int flag)
{
	specular_texture_button->setEnabled(flag);
	specular_texture_dir_lineEdit->setEnabled(flag);
	specular_color_button->setEnabled(!flag);
	specular_color_r->setEnabled(!flag);
	specular_color_g->setEnabled(!flag);
	specular_color_b->setEnabled(!flag);

	if (!flag)
	{
		updateFromSpecularEditLine();
		_selectedMaterial->useTextureMap[SpecularMapIndex] = 0;
	}
	else
		_selectedMaterial->useTextureMap[SpecularMapIndex] = 1;
}

void BaseMaterialEditWidget::enableBumpTexture(int flag)
{
	bump_texture_button->setEnabled(flag);
	bump_texture_dir_lineEdit->setEnabled(flag);
	if (!flag)
		_selectedMaterial->useTextureMap[BumpMapIndex] = 0;
	else
		_selectedMaterial->useTextureMap[BumpMapIndex] = 1;
}

void BaseMaterialEditWidget::handleTextureAmbientBtn()
{
	if (_selectedMaterial == nullptr)
		return;
	QString filename = QFileDialog::getOpenFileName(this, tr(""), "./texture/", tr(""));
	if (filename.isEmpty())
		return;
	_selectedMaterial->readTextureMap(filename, AmbientMapIndex);
	refreshTextureBar();
}

void BaseMaterialEditWidget::handleTextureDiffuseBtn()
{
	if (_selectedMaterial == nullptr)
		return;
	QString filename = QFileDialog::getOpenFileName(this, tr(""), "./texture/", tr(""));
	if (filename.isEmpty())
		return;
	_selectedMaterial->readTextureMap(filename, DiffuseMapIndex);
	refreshTextureBar();
}

void BaseMaterialEditWidget::handleTextureSpecularBtn()
{
	if (_selectedMaterial == nullptr)
		return;
	QString filename = QFileDialog::getOpenFileName(this, tr(""), "./texture/", tr(""));
	if (filename.isEmpty())
		return;
	_selectedMaterial->readTextureMap(filename, SpecularMapIndex);
	refreshTextureBar();
}

void BaseMaterialEditWidget::handleTextureBumpBtn()
{
	if (_selectedMaterial == nullptr)
		return;
	QString filename = QFileDialog::getOpenFileName(this, tr(""), "./texture/", tr(""));
	if (filename.isEmpty())
		return;
	_selectedMaterial->readTextureMap(filename, BumpMapIndex);
	refreshTextureBar();
}