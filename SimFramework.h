#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_SimFramework.h"
#include <QToolButton>
#include "Ui\BaseMainWidget.h"
class SimFramework : public QMainWindow
{
	Q_OBJECT

public:
	SimFramework(QWidget *parent = Q_NULLPTR);
private:
	Ui::SimFrameworkClass _ui;
};
