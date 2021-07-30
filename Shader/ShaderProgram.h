#pragma once
#ifndef SHADER_PROGRAM_H
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions>
#include <QOpenGLBuffer>
#include <iostream>
using namespace std;
class ShaderProgram : public QOpenGLShaderProgram
{
public:
	ShaderProgram();
	~ShaderProgram();
	bool initProgram(char* vshader, char* fshader, char* gshader = nullptr, bool useGeom = false);
};
#endif // !SHADER_PROGRAM_H
