#include "BaseLighting.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "Commom\DataConversion.h"

void BaseLight::transferToShader(QOpenGLShaderProgram * program, unsigned int lightId)
{
	QString idStr;
	idStr.setNum(lightId);
	program->setUniformValue(("directionalLight[" + idStr.toStdString() + "].position").c_str(), lightPos);
	program->setUniformValue(("directionalLight[" + idStr.toStdString() + "].ambientColor").c_str(), getOpenglRGBA(ambientColor));
	program->setUniformValue(("directionalLight[" + idStr.toStdString() + "].diffuseColor").c_str(), getOpenglRGBA(diffuseColor));
	program->setUniformValue(("directionalLight[" + idStr.toStdString() + "].specularColor").c_str(), getOpenglRGBA(specularColor));
}

void DirectionalLight::transferToShader(QOpenGLShaderProgram * program, unsigned int lightId)
{
	QString idStr;
	idStr.setNum(lightId);
	program->setUniformValue(("directionalLight["+ idStr.toStdString()+"].position").c_str(), lightPos);
	program->setUniformValue(("directionalLight[" + idStr.toStdString() + "].ambientColor").c_str(), getOpenglRGBA(ambientColor));
	program->setUniformValue(("directionalLight[" + idStr.toStdString() + "].diffuseColor").c_str(), getOpenglRGBA(diffuseColor));
	program->setUniformValue(("directionalLight[" + idStr.toStdString() + "].specularColor").c_str(), getOpenglRGBA(specularColor));
	program->setUniformValue(("directionalLight[" + idStr.toStdString() + "].direction").c_str(), lightdir);

}

void PointLight::transferToShader(QOpenGLShaderProgram * program, unsigned int lightId)
{
	QString idStr;
	idStr.setNum(lightId);
	program->setUniformValue(("pointLight[" + idStr.toStdString() + "].position").c_str(), lightPos);
	program->setUniformValue(("pointLight[" + idStr.toStdString() + "].ambientColor").c_str(), getOpenglRGBA(ambientColor));
	program->setUniformValue(("pointLight[" + idStr.toStdString() + "].diffuseColor").c_str(), getOpenglRGBA(diffuseColor));
	program->setUniformValue(("pointLight[" + idStr.toStdString() + "].specularColor").c_str(), getOpenglRGBA(specularColor));
	program->setUniformValue(("pointLight[" + idStr.toStdString() + "].constant").c_str(), constant);
	program->setUniformValue(("pointLight[" + idStr.toStdString() + "].linear").c_str(), linear);
	program->setUniformValue(("pointLight[" + idStr.toStdString() + "].quadratic").c_str(), quadratic);
}


BaseLight* convertBaseLightingPtr(PointLight* p)
{
	return dynamic_cast<BaseLight*>(p);
}

DirectionalLight* convertDirectionalLightingPtr(BaseLight* p)
{
	if (p->getLightType() == LightType::DIRECTIONAL_LIGHT)
		return dynamic_cast<DirectionalLight*>(p);
	return nullptr;
}

BaseLight* convertBaseLightingPtr(DirectionalLight* p)
{
	return dynamic_cast<BaseLight*>(p);
}

PointLight* convertPointLightingPtr(BaseLight* p)
{
	if (p->getLightType() == LightType::POINT_LIGHT)
		return dynamic_cast<PointLight*>(p);
	return nullptr;
}

int SceneLighting::addDirectionalLight(DirectionalLight light, bool replaceLast)
{
	if (_dLNum < MAX_LIGHT_NUM)
	{
		int L = _dLNum;
		_directionalLights[L] = light;
		_dLNum++;
		return L;
	}
	else if(replaceLast){
		_directionalLights[_dLNum - 1] = light;
		return _dLNum - 1;
	}
	else return -1;
}

int SceneLighting::removeDirectionalLight(const unsigned int lightId)
{
	if (lightId >= _dLNum)
		return -1;

	if (lightId == 0)
	{
		_directionalLights[0] = _directionalLights[1];
		_directionalLights[1] = _directionalLights[2];
		_directionalLights[2] = DirectionalLight();
		_dLNum--;
	}
	else if (lightId == 1)
	{
		_directionalLights[1] = _directionalLights[2];
		_directionalLights[2] = DirectionalLight();
		_dLNum--;
	}
	else if (lightId == 2)
	{
		_directionalLights[2] = DirectionalLight();
		_dLNum--;
	}

	return 	_dLNum;
}

int SceneLighting::addPointLight(PointLight light, bool replaceLast)
{
	if (_pLNum < MAX_LIGHT_NUM)
	{
		int L = _dLNum;
		_pointLights[L] = light;
		_dLNum++;
		return L;
	}
	else if (replaceLast) {
		_pointLights[_pLNum - 1] = light;
		return _pLNum - 1;
	}
	else return -1;
}

int SceneLighting::removePointLight(const unsigned int lightId)
{
	if (lightId >= _pLNum)
		return -1;

	if (lightId == 1)
	{
		_pointLights[1] = _pointLights[2];
		_pointLights[2] = PointLight();
		_pLNum--;
	}
	else if (lightId == 2)
	{
		_pointLights[2] = PointLight();
		_pLNum--;
	}
	return 	_pLNum;
}

DirectionalLight * SceneLighting::getDirectionalLight(const unsigned int lightId)
{
	if(lightId >= _dLNum)
		return nullptr;
	return &(_directionalLights[lightId]);
}

PointLight * SceneLighting::getPointLight(const unsigned int lightId)
{
	if (lightId >= _pLNum)
		return nullptr;
	return &(_pointLights[lightId]);
}

void SceneLighting::transferToShader(QOpenGLShaderProgram * program)
{
	program->setUniformValue("directionalLightNum", _dLNum);
	program->setUniformValue("pointLightNum", _pLNum);

	for (int i = 0; i < _dLNum; i++)
		_directionalLights[i].transferToShader(program, i);
	for (int i = 0; i < _pLNum; i++)
		_pointLights[i].transferToShader(program, i);

	if(useShadow())
		program->setUniformValue("enableShadow", 1);
	else program->setUniformValue("enableShadow", 0);
}
