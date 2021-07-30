#pragma once
#ifndef SCREEN_FONT_TYPE_H
#define  SCREEN_FONT_TYPE_H
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLTexture>
#include <QGLViewer/camera.h>
#include <ft2build.h>
#include <ft2build.h>
#include <freetype/freetype.h>
#include <freetype/ftglyph.h>
#include <freetype/ftoutln.h>
#include <freetype/fttrigon.h>
#include FT_FREETYPE_H
#include <vector>
#include <map>
#include <string>
#include "MatrixCore.h"

namespace ScreenTextPainter
{
	using namespace std;
	using namespace Eigen;
	
	struct Character
	{
		GLuint textureId;
		Vector2i size;
		Vector2i bearing;
		FT_Pos advance;
	};

	struct FontCoords
	{
		QVector2D pos;
		QVector2D tex;
	};

	class TextPainter /*: protected QOpenGLFunctions*/
	{
	public:
		TextPainter();
		~TextPainter();
		void generateFont(QOpenGLFunctions* f, string filename = "Resources/times.ttf", float height = 36);
		void renderText(QOpenGLShaderProgram* program, string text, GLfloat x, GLfloat y, GLfloat scale, QVector3D color);

		map <GLchar , Character > charList;
		FT_Library ft;
		FT_Face face;
		//render
		QOpenGLBuffer arrayBuf;

		vector<FontCoords> fontData;
		vector<QOpenGLTexture*> textures;
		vector<GLushort> index;
	protected:
		QOpenGLFunctions* _opengl;
	};


};



#endif