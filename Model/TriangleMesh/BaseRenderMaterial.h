#pragma once
#ifndef BASE_RENDER_MATERIAL_H
#define BASE_RENDER_MATERIAL_H

#include <QOpenGLTexture>
#include <vector>
#include <array>
#include <QColor>
#include <QOpenGLShaderProgram>
#include "MatrixCore.h"

#define AmbientMapIndex 0
#define DiffuseMapIndex 1
#define SpecularMapIndex 2
#define BumpMapIndex 3

class BaseRenderMaterial
{
public:
	BaseRenderMaterial();
	QOpenGLTexture* getTextureMap(int mapIndex) { return _textureMap[mapIndex]; }
	bool readTextureMap(const QString filename, int mapIndex);
	void transferToShader(QOpenGLShaderProgram * program);

	QColor ambient;
	QColor diffuse;
	QColor specular;
	qeal shinness;

	std::array<std::string, 4> textureMapFilename;
	std::array<int, 4> useTextureMap;

	void enableTransparent(bool enable) {_transparent = enable;}
	bool isTransparent() {return _transparent;}
protected:
	std::array<QOpenGLTexture*, 4> _textureMap;
	bool _transparent;
};
#endif