#include "DataConversion.h"

QVector4D getOpenglRGBA(const QColor color)
{
	return QVector4D(color.red() / 255.0, color.green() / 255.0, color.blue() / 255.0, color.alpha() / 255.0);
}

QColor getQColorRGBA(const QVector4D color)
{
	if (color[0] <= 1.0 && color[1] <= 1.0 && color[2] <= 1.0 && color[3] <= 1.0)
	{		
		return QColor(int(color[0] * 255.0), int(color[1] * 255.0), int(color[2] * 255.0), int(color[3] * 255.0));
	}
	else
	{
		return QColor(int(color[0]), int(color[1]), int(color[2]), int(color[3]));
	}		
}

QVector3D getOpenglRGB(const QColor color)
{
	return QVector3D(color.red() / 255.0, color.green() / 255.0, color.blue() / 255.0);
}

QColor getQColorRGB(const QVector3D color)
{
	if (color[0] <= 1.0 && color[1] <= 1.0 && color[2] <= 1.0)
	{
		return QColor(int(color[0] * 255.0), int(color[1] * 255.0), int(color[2] * 255.0), 1.0);
	}
	else
	{
		return QColor(int(color[0]), int(color[1]), int(color[2]), 1.0);
	}
}
