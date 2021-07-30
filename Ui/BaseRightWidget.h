#pragma warning(disable:4068)
#ifndef BASE_RIGHT_WIDGET_H
#define BASE_RIGHT_WIDGET_H

#include <QtWidgets/qapplication.h>
#include <QtWidgets>
#include <QResizeEvent>
#include "Scene\BaseScene.h"
#include "ui_RightFormWidget.h"
#include "ColorSelectionPanel.h"
#include "BaseMaterialEditWidget.h"

class BaseRightWidget : public QWidget, public Ui::RightWidget
{
	Q_OBJECT
public:
	BaseRightWidget();
	~BaseRightWidget();
	void resizeEvent(QResizeEvent* size);

	void setBaseScene(BaseScene* scene);

signals:
	void openColorPanel(QColor*);
	void needInitSimulator();
public slots:
	void getColorFromPanel(QColor*);

	void bindSceneTab();
	void bindCameraGroup();
	void bindLightingGroup();
	void bindFloorGroup();

	void bindModelTab();
	void bindModelList();
	void bindMeshGroup();

	// scene Tab
	//ColorSelection
	//camera group
	void useCameraOrthoMode();
	void lockCamera();
	void showCameraFOV();

	//lighting group
	void setMainLightPos();
	void setMainLightAttenuation();
	void pressMainLightAmbientBtn();
	void pressMainLightDiffuseBtn();
	void pressMainLightSpecularBtn();

	void setMainLightAmbientRGB();
	void setMainLightAmbientEditLine();
	void setMainLightDiffuseRGB();
	void setMainLightDiffuseEditLine();
	void setMainLightSpecularRGB();
	void setMainLightSpecularEditLine();
	void enableShadow();

	//Floor group
	void handleFloorMaterialBtn();
	void hideFloor();
	void scaleFloorXZ();
	void setFloorY();

	// Physical Environment
	void refreshPhysicalEnvironmentGroup();
	void enableGrivaty();
	void setGrivaty();
	void setTimeStep();
	void handleAddExtraForceBtn();

	// Model Tab
	// Model List
	void refreshModelList();
	void selectModel(int);
	void clickModel(QListWidgetItem * item);
	// Mesh Gruop
	void refreshMeshGroup();
	void handleMeshMaterialBtn();
	void hideMesh();
	void setMeshScale();
	void setMeshTranslation();
	void setMeshRotation();

	// Volumetric Mesh Gruop
	void refreshVolumetricMeshGroup();
	void hideTetMesh();
	void setTetMeshTransparent();

	void selectTetElementSet(int);
	void clickTetElementSetList(QListWidgetItem * item);
	// Medial Mesh Gruop
	void refreshMedialMeshGroup();
	void hideMedialMesh();
	void setMedialMeshTransparent();

	void selectMedialMeshSet(int);
	void clickMedialMeshSetList(QListWidgetItem * item);
protected:
	BaseScene* _scene;
	ColorSelectionPanel* colorPanel;
	BaseMaterialEditWidget* materialPanel;

	int _selectModel;
	

	// debug
	void debugModifiledFile(const std::string filename);
};




#endif // !BASE_RIGHT_WIDGET_H
