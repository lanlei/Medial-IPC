#pragma once
#ifndef DATA_CONVERSION_H
#define DATA_CONVERSION_H
#include<QVector3D>
#include<QVector4D>
#include<QColor>

QVector4D getOpenglRGBA(const QColor color);
QColor getQColorRGBA(const QVector4D color);
QVector3D getOpenglRGB(const QColor color);
QColor getQColorRGB(const QVector3D color);

//


#endif

////

