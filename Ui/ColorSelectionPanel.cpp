#include "ColorSelectionPanel.h"

ColorSelectionPanel::ColorSelectionPanel():
	_boundColor(nullptr)
{
	setupUi(this);
	setFixedSize(QSize(680, 350));

	gridLayout = new QGridLayout(this);
	gridLayout->setObjectName(QStringLiteral("gridLayout"));
	setWindowTitle(QString("ColorSelection"));

	colorPanelBtn = new ColorPanelBtn(this);
	colorPanelHSB = new ColorPanelHSB(this);
	colorButton = new ColorButton(this);
	colorButton->setFixedWidth(180);
	colorButton->setCanMove(false);
	colorButton->setEnabled(false);
	redSlider = new QSlider(Qt::Horizontal, this);
	greenSlider = new QSlider(Qt::Horizontal, this);
	blueSlider = new QSlider(Qt::Horizontal, this);
	alphaSlider = new QSlider(Qt::Horizontal, this);
	redSlider->setMinimum(0);
	redSlider->setMaximum(255);
	greenSlider->setMinimum(0);
	greenSlider->setMaximum(255);
	blueSlider->setMinimum(0);
	blueSlider->setMaximum(255);
	alphaSlider->setMinimum(0);
	alphaSlider->setMaximum(255);

	int rows = 4;
	int cols = 10;

	gridLayout->addWidget(colorPanelHSB, 0, 0, 1, 9);
	gridLayout->addWidget(colorPanelBtn, 0, 9, 5, 1);

	gridLayout->addWidget(colorButton, 1, 0, 4, 5);

	gridLayout->addWidget(r_color_lineEdit, 1, 5, 1, 2);
	gridLayout->addWidget(redSlider, 1, 7, 1, 1);
	gridLayout->addWidget(r_color_label, 1, 8, 1, 1);

	gridLayout->addWidget(g_color_lineEdit, 2, 5, 1, 2);
	gridLayout->addWidget(greenSlider, 2, 7, 1, 1);
	gridLayout->addWidget(g_color_label, 2, 8, 1, 1);

	gridLayout->addWidget(b_color_lineEdit, 3, 5, 1, 2);
	gridLayout->addWidget(blueSlider, 3, 7, 1, 1);
	gridLayout->addWidget(b_color_label, 3, 8, 1, 1);

	gridLayout->addWidget(alpha_color_lineEdit, 4, 5, 1, 2);
	gridLayout->addWidget(alphaSlider, 4, 7, 1, 1);
	gridLayout->addWidget(alpha_color_label, 4, 8, 1, 1);

	gridLayout->setHorizontalSpacing(10);
	gridLayout->setVerticalSpacing(10);
	gridLayout->setContentsMargins(10, 10, 10, 10);
	setLayout(gridLayout);

	connect(colorPanelBtn, SIGNAL(colorChanged(const QColor &)), this, SLOT(updateFromColorPanelBtn(const QColor &)));
	connect(colorPanelHSB, SIGNAL(colorChanged(const QColor &, double, double)), this, SLOT(updateFromColorPanelHSB(const QColor&, double, double)));

	connect(r_color_lineEdit, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromRedEditLine()));
	connect(g_color_lineEdit, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromGreenEditLine()));
	connect(b_color_lineEdit, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromBlueEditLine()));
	connect(alpha_color_lineEdit, SIGNAL(textChanged(const QString &)), this, SLOT(updateFromAlphaEditLine()));

	connect(redSlider, SIGNAL(valueChanged(int)), this, SLOT(updateFromRedSlider(int)));
	connect(greenSlider, SIGNAL(valueChanged(int)), this, SLOT(updateFromGreenSlider(int)));
	connect(blueSlider, SIGNAL(valueChanged(int)), this, SLOT(updateFromBlueSlider(int)));
	connect(alphaSlider, SIGNAL(valueChanged(int)), this, SLOT(updateFromAlphaSlider(int)));
}

void ColorSelectionPanel::updateFromColorPanelBtn(const QColor &color)
{
	_color = color;
	refreshWidget();
}

void ColorSelectionPanel::updateFromColorPanelHSB(const QColor & color, double hue, double sat)
{
	_color = color;
	refreshWidget();
}

void ColorSelectionPanel::updateFromRedEditLine()
{
	_color.setRed(r_color_lineEdit->text().toInt());
	refreshWidget();
}

void ColorSelectionPanel::updateFromGreenEditLine()
{
	_color.setGreen(g_color_lineEdit->text().toInt());
	refreshWidget();
}


void ColorSelectionPanel::updateFromBlueEditLine()
{
	_color.setBlue(b_color_lineEdit->text().toInt());
	refreshWidget();
}

void ColorSelectionPanel::updateFromAlphaEditLine()
{
	_color.setAlpha(alpha_color_lineEdit->text().toInt());
	refreshWidget();
}

void ColorSelectionPanel::updateFromRedSlider(int r)
{
	_color.setRed(r);
	refreshWidget();
}

void ColorSelectionPanel::updateFromGreenSlider(int g)
{
	_color.setGreen(g);
	refreshWidget();
}

void ColorSelectionPanel::updateFromBlueSlider(int b)
{
	_color.setBlue(b);
	refreshWidget();
}

void ColorSelectionPanel::updateFromAlphaSlider(int a)
{
	_color.setAlpha(a);
	refreshWidget();
}

void ColorSelectionPanel::closeEvent(QCloseEvent * event)
{
	emit panelClose();
}

void ColorSelectionPanel::refreshWidget()
{
	colorButton->setNormalColor(_color);
	r_color_lineEdit->setText(QString::number(_color.red()));
	g_color_lineEdit->setText(QString::number(_color.green()));
	b_color_lineEdit->setText(QString::number(_color.blue()));
	alpha_color_lineEdit->setText(QString::number(_color.alpha()));
	redSlider->setValue(_color.red());
	greenSlider->setValue(_color.green());
	blueSlider->setValue(_color.blue());
	alphaSlider->setValue(_color.alpha());

	if (_boundColor)
	{
		*_boundColor = _color;
		emit colorChange(_boundColor);
	}
}

void ColorSelectionPanel::bindColor(QColor* color)
{
	if (color)
	{
		_boundColor = color;
		_color = *_boundColor;
	}
	else _color = QColor(126, 126, 126);

	refreshWidget();
}
