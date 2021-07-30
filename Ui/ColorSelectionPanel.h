#pragma once
#ifndef  COLOR_SELECTION_PANEL_H
#define COLOR_SELECTION_PANEL_H

#include <QtWidgets/qapplication.h>
#include <QtWidgets>
#include "ui_ColorSelectionWidget.h"
#include <QCloseEvent>

#include "colorpanelbtn.h"
#include "colorbutton.h"
#include "colorpanelhsb.h"
#include <QSlider>
#include <QColor>

class ColorSelectionPanel : public QWidget, public Ui::ColorSelectionWidget
{
	Q_OBJECT
public:
	ColorSelectionPanel();
	~ColorSelectionPanel(){}

	void refreshWidget();
signals:
	void colorChange(QColor* c);
	void panelClose();
public slots:
	void updateFromColorPanelBtn(const QColor &color); 
	void updateFromColorPanelHSB(const QColor &color, double hue, double sat);

	void updateFromRedEditLine();
	void updateFromGreenEditLine();
	void updateFromBlueEditLine();
	void updateFromAlphaEditLine();

	void updateFromRedSlider(int);
	void updateFromGreenSlider(int);
	void updateFromBlueSlider(int);
	void updateFromAlphaSlider(int);
protected:
	void closeEvent(QCloseEvent *event);
public:
	QGridLayout *gridLayout;

	ColorPanelBtn* colorPanelBtn;
	ColorPanelHSB* colorPanelHSB;
	ColorButton* colorButton;
	QSlider* redSlider;
	QSlider* greenSlider;
	QSlider* blueSlider;
	QSlider* alphaSlider;

	void bindColor(QColor* color);
	QColor getColor() { return _color; }
private:
	QColor _color;
	QColor* _boundColor;
};



#endif // ! COLOR_SELECTION_WIDGET_H
