#include "ShaderProgram.h"

ShaderProgram::ShaderProgram()
{

}

ShaderProgram::~ShaderProgram()
{

}

bool ShaderProgram::initProgram(char* vshader, char* fshader, char* gshader, bool useGeom)
{
	// compile vertex shader
	if (!addShaderFromSourceFile(QOpenGLShader::Vertex, vshader))
	{
		cout << log().toStdString() << endl;
		return false;
	}

	// compile fragment shader
	if (!addShaderFromSourceFile(QOpenGLShader::Fragment, fshader))
	{
		cout << log().toStdString() << endl;
		return false;
	}

	// compile geom shader
	if (useGeom)
	{
		if (!addShaderFromSourceFile(QOpenGLShader::Geometry, gshader))
		{
			cout << log().toStdString() << endl;
			return false;
		}
	}


	// Link shader pipeline
	if (!link())
	{
		cout << log().toStdString() << endl;
		return false;
	}
	return true;
}