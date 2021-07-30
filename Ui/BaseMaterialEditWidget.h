#pragma once
#ifndef  BASE_METRIAL_EDIT_WIDGET_H
#define BASE_METRIAL_EDIT_WIDGET_H
#include<QColorDialog>
#include <QtWidgets/qapplication.h>
#include <QtWidgets>
#include "ui_MaterialEditWidget.h"
#include "ColorSelectionPanel.h"
#include "Model\TriangleMesh\BaseRenderMaterial.h"

class BaseMaterialEditWidget : public QWidget, public Ui::MaterialEditWidget
{
	Q_OBJECT
public:
	BaseMaterialEditWidget();
	~BaseMaterialEditWidget(){}

	void openWidget(BaseRenderMaterial* materials, int size);
signals:
	void openColorPanel(QColor*);
public slots:
	void selectMaterial(int);
	void handleColorPanel(QColor*);

	void refreshAmbient();
	void updateFromAmbientEditLine();

	void refreshDiffuse();
	void updateFromDiffuseEditLine();

	void refreshSpecular();
	void updateFromSpecularEditLine();

	void handleAmbientBtn();
	void handleDiffuseBtn();
	void handleSpecularBtn();

	void handleTransparentCheckBox();

	void refreshShinness();
	void updateFromShinnessEditLine();
	void updateFromShinnessSlider(int);

	void refreshTextureBar();

	void enableAmbientTexture(int);
	void enableDiffuseTexture(int);
	void enableSpecularTexture(int);
	void enableBumpTexture(int);

	void handleTextureAmbientBtn();
	void handleTextureDiffuseBtn();
	void handleTextureSpecularBtn();
	void handleTextureBumpBtn();
public:

	ColorSelectionPanel* colorPanel;
private:

	BaseRenderMaterial* _materials;
	BaseRenderMaterial* _selectedMaterial;
	int _materialNum;
};

#endif
