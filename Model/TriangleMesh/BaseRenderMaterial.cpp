#include "BaseRenderMaterial.h"
#include <iostream>
#include <string>
#include "Commom\DataConversion.h"
BaseRenderMaterial::BaseRenderMaterial()
{
//	ambient = QColor(64, 64, 166);
	ambient = QColor(240, 128, 128);	
	diffuse = QColor(76, 76, 76);
	specular = QColor(76, 76, 76);
	shinness = 6.0;
	
	for (int i = 0; i < 4; i++)
	{
		textureMapFilename[i] = std::string("");
		useTextureMap[i] = 0;
		_textureMap[i] = nullptr;
	}
	_transparent = false;
}

bool BaseRenderMaterial::readTextureMap(const QString filename, int mapIndex)
{
	std::string img_filename = filename.toStdString();
	QImage image = QImage(img_filename.c_str());
	
	if (image.isNull())
	{
		std::cout << "The texture " << img_filename << " don't exist !" << std::endl;
		return false;
	}

	if (_textureMap[mapIndex] != nullptr)
		delete _textureMap[mapIndex];
	_textureMap[mapIndex] = new QOpenGLTexture(QImage(img_filename.c_str()).mirrored());

	if (_textureMap[mapIndex]->isCreated())
	{
		_textureMap[mapIndex]->setMinificationFilter(QOpenGLTexture::Nearest);
		_textureMap[mapIndex]->setMagnificationFilter(QOpenGLTexture::Linear);
		_textureMap[mapIndex]->setWrapMode(QOpenGLTexture::Repeat);
		textureMapFilename[mapIndex] = img_filename;
		useTextureMap[mapIndex] = 1;
	}else useTextureMap[mapIndex] = 0;
	return true;
}

void BaseRenderMaterial::transferToShader(QOpenGLShaderProgram * program)
{
	program->setUniformValue("material.ambient", getOpenglRGBA(ambient));
	program->setUniformValue("material.diffuse", getOpenglRGBA(diffuse));
	program->setUniformValue("material.specular", getOpenglRGBA(specular));
	program->setUniformValue("material.shininess", getOpenglRGBA(shinness));

	if (useTextureMap[AmbientMapIndex] == 1 && _textureMap[AmbientMapIndex]->isCreated())
	{
		program->setUniformValue("material.useAmbientMap", useTextureMap[AmbientMapIndex]);
		getTextureMap(AmbientMapIndex)->bind(1);
	}else program->setUniformValue("material.useAmbientMap", 0);
		
	if (useTextureMap[DiffuseMapIndex] == 1 && _textureMap[DiffuseMapIndex]->isCreated())
	{
		program->setUniformValue("material.useDiffuseMap", useTextureMap[DiffuseMapIndex]);
		getTextureMap(DiffuseMapIndex)->bind(2);
	}else program->setUniformValue("material.useDiffuseMap", 0);
		
	if (useTextureMap[SpecularMapIndex] == 1 && _textureMap[SpecularMapIndex]->isCreated())
	{
		program->setUniformValue("material.useSpecularMap", useTextureMap[SpecularMapIndex]);
		getTextureMap(SpecularMapIndex)->bind(3);
	}else program->setUniformValue("material.useSpecularMap", 0);
		
	if (useTextureMap[BumpMapIndex] == 1 && _textureMap[BumpMapIndex]->isCreated())
	{
		program->setUniformValue("material.useBumpMap", useTextureMap[BumpMapIndex]);
		getTextureMap(BumpMapIndex)->bind(4);
	}else program->setUniformValue("material.useBumpMap", 0);
		
}
