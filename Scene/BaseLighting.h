#pragma once
#ifndef BASE_LIGHTING_H
#define BASE_LIGHTING_H

#include <QVector3D>
#include <array>
#include <QOpenGLShaderProgram>
#include <QColor>
#include "MatrixCore.h"

enum LightType
{
	NONE = 0,
	DIRECTIONAL_LIGHT = 1,
	POINT_LIGHT = 2,
};

class PointLight;
class DirectionalLight;

class BaseLight
{
public:
	BaseLight(LightType type = NONE, QVector3D pos = QVector3D(0, 5.0, 5.0), QColor ac = QColor(153, 153, 153), QColor dc = QColor(128, 128, 128), QColor sc = QColor(128, 128, 128)):
		_type(type),
		lightPos(pos),
		ambientColor(ac),
		diffuseColor(dc),
		specularColor(sc)
	{}
	virtual ~BaseLight(){}

	LightType getLightType() {return _type;}

	virtual void transferToShader(QOpenGLShaderProgram* program, unsigned int lightId);
	
	QVector3D lightPos;
	QColor ambientColor;
	QColor diffuseColor;
	QColor specularColor;

protected:
	LightType _type;
};

class DirectionalLight : public BaseLight
{
public:
	DirectionalLight(QVector3D dir = QVector3D(0, -1, 1).normalized(), QVector3D pos = QVector3D(0, 5.0, 5.0), QColor ac = QColor(153, 153, 153), QColor dc = QColor(128, 128, 128), QColor sc = QColor(128, 128, 128)):
		BaseLight(LightType::DIRECTIONAL_LIGHT, pos, ac, dc, sc),
		lightdir(dir.normalized())
	{}
	void transferToShader(QOpenGLShaderProgram* program, unsigned int lightId);

	QVector3D lightdir;
};

BaseLight* convertBaseLightingPtr(DirectionalLight* p);
PointLight* convertPointLightingPtr(BaseLight* p);

class PointLight : public BaseLight
{
public:
	PointLight(QVector3D pos = QVector3D(0, 5.0, 5.0), QColor ac = QColor(153, 153, 153), QColor dc = QColor(128, 128, 128), QColor sc = QColor(128, 128, 128)) :
		BaseLight(LightType::POINT_LIGHT, pos, ac, dc, sc),
		constant(1.0),
		linear(0.00),
		quadratic(0.0) {}

	void transferToShader(QOpenGLShaderProgram* program, unsigned int lightId);

	//Attenuation
	float constant;
	float linear;
	float quadratic;
};

BaseLight* convertBaseLightingPtr(PointLight* p);
DirectionalLight* convertDirectionalLightingPtr(BaseLight* p);

#define MAX_LIGHT_NUM 3

class SceneLighting
{
public:
	SceneLighting()
	{
		_dLNum = 0;
		_pLNum = 1;
		shadow = true;
	}
	~SceneLighting(){}

	int addDirectionalLight(DirectionalLight light, bool replaceLast = false);
	int removeDirectionalLight(const unsigned int lightId);

	int addPointLight(PointLight light, bool replaceLast = false);
	int removePointLight(const unsigned int lightId);

	DirectionalLight* getDirectionalLight(const unsigned int lightId = 0);
	PointLight* getPointLight(const unsigned int lightId = 0);

	void transferToShader(QOpenGLShaderProgram* program);

	QVector3D getMainLightPosition() {
		return _pointLights
			[0].lightPos;
	}

	void enableShadow(bool enable) { shadow = enable; };
	bool useShadow() { return shadow; }
protected:
	// _pointLights[0] is main lighting
	int _pLNum;
	std::array<PointLight, MAX_LIGHT_NUM> _pointLights;
	int _dLNum;
	std::array<DirectionalLight, MAX_LIGHT_NUM> _directionalLights;
	bool shadow;
};


#endif // !BASE_LIGHTING_H
