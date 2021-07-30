#include "BaseScene.h"
#include<direct.h>

BaseScene::BaseScene():_sim(nullptr)
{
	QGLFormat glFormat;
	glFormat.setVersion(3, 2);
	glFormat.setProfile(QGLFormat::CoreProfile);

	QSurfaceFormat glSurfaceFormat;
	glSurfaceFormat.setSamples(32);
	setFormat(glSurfaceFormat);

	_selectedTetPointId = -1;
	_mouseForce.setZero();

}

BaseScene::BaseScene(QWidget* parent) :
	QGLViewer(parent)
{

}

BaseScene::~BaseScene()
{
	if (_camera == nullptr)
		delete _camera;
}

bool BaseScene::setFromConfigFile(const std::string filename, TiXmlElement * item)
{
	return false;
}

void BaseScene::bindSimulator(BaseSimulator * sim)
{
	if (_sim)
		free(_sim);
	_sim = sim;
	emit changeSimulator();
}

void BaseScene::draw()
{
	renderDepthToBuffer();

	glViewport(0, 0, width(), height());
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	setBackgroundColor(QColor(255, 255, 255));

	_phongProgram.bind();
	_phongProgram.setUniformValue("projectMatrix", ((StandardCamera*)camera())->getProjectionQMatrix());
	_phongProgram.setUniformValue("viewMatrix", ((StandardCamera*)camera())->getViewQMatrix());
	QVector3D lightPos = _lights.getMainLightPosition();

	_lightCamera.setPosition(qglviewer::Vec(lightPos[0], lightPos[1], lightPos[2]));
	_lightCamera.lookAt(qglviewer::Vec(0.0, 0.0, 0.0));

	_phongProgram.setUniformValue("lightSpaceMatrix", _lightCamera.getProjectionViewQMatrix());
	QVector3D viewPos = QVector3D(camera()->position().x, camera()->position().y, camera()->position().z);
	_phongProgram.setUniformValue("viewPos", viewPos);

	_lights.transferToShader(&_phongProgram);

	_phongProgram.setUniformValue("shadowMap", 0);
	_phongProgram.setUniformValue("ambientMap", 1);
	_phongProgram.setUniformValue("diffuseMap", 2);
	_phongProgram.setUniformValue("specularMap", 3);
	_phongProgram.setUniformValue("bumpMap", 4);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, _depthMapId);
	
	renderScene(&_phongProgram);

	glBindTexture(GL_TEXTURE_2D, 0);
	_phongProgram.release();

	renderScene();
	update();
}

void BaseScene::animate()
{
	if (!_sim)
		return;

	_status.animationTimer.start();

	_sim->animate(_status.animationLifeTime);

	qeal time = _status.animationTimer.nsecsElapsed() / 1e6;
	_status.animationCost += time;
	_status.animationLifeTime++;

	if (_status.recordAnimation)
		emit saveAnimation(true);
}

void BaseScene::postDraw()
{	
	renderSelectionRect();
	renderSelectionPoints();
	renderForceLine();
	renderCornerAxis();
	renderInfoOnScene();
	update();
}

void BaseScene::init()
{
	initializeOpenGLFunctions();
	glEnable(GL_MULTISAMPLE);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	initShaderProgram();
	initCamera();
	initShadowFramebuffer();
	initAnimation();

	_floor = new Floor();
	emit initSceneSignal();

	select_tet_group = 1;
}

void BaseScene::mousePressEvent(QMouseEvent* e)
{
	_rectangle = QRect(e->pos(), e->pos());
	if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::ShiftModifier))
		_selectionMode = SHIFT_SELECT;
	else if ((e->button() == Qt::RightButton) && (e->modifiers() == Qt::ShiftModifier))
		_selectionMode = SHIFT_REMOVE;
	else if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::AltModifier))
		_selectionMode = ALT_SELECT;
	else if ((e->button() == Qt::RightButton) && (e->modifiers() == Qt::AltModifier))
		_selectionMode = ALT_REMOVE;
	else if ((e->button() == Qt::LeftButton) && (e->modifiers() == Qt::ControlModifier))
	{
		_selectionMode = CTRL_SELECT;
		setSelectRegionWidth(30);
		setSelectRegionHeight(30);
		select(e->pos());
		_forceLine.setP1(e->pos());
	}
	else if ((e->button() == Qt::RightButton) && (e->modifiers() == Qt::ControlModifier))
		_selectionMode = CTRL_REMOVE;
	else
		QGLViewer::mousePressEvent(e);
	update();
}

void BaseScene::mouseMoveEvent(QMouseEvent* e)
{
	if (_selectionMode != NONE_SELECT && _status.enableSelectSurfacePoints)
	{
		_rectangle.setBottomRight(e->pos());
		update();
	}
	else if (_selectionMode != NONE_SELECT && _status.enableSelectTetNodes)
	{
		_rectangle.setBottomRight(e->pos());
		update();
	}
	else if (_selectionMode == CTRL_SELECT)
	{
		_forceLine.setP2(e->pos());
	}
	else QGLViewer::mouseMoveEvent(e);
}

void BaseScene::mouseReleaseEvent(QMouseEvent *e)
{
	if (_selectionMode != NONE_SELECT && _status.enableSelectSurfacePoints)
	{
		// to do
		_rectangle = _rectangle.normalized();
		setSelectRegionWidth(_rectangle.width());
		setSelectRegionHeight(_rectangle.height());
		select(_rectangle.center());
		_selectionMode = NONE_SELECT;
		_rectangle = QRect();
		update();
	}
	else if (_selectionMode != NONE_SELECT && _status.enableSelectTetNodes)
	{
		_rectangle = _rectangle.normalized();
		setSelectRegionWidth(_rectangle.width());
		setSelectRegionHeight(_rectangle.height());
		select(_rectangle.center());
		_selectionMode = NONE_SELECT;
		_rectangle = QRect();
		update();
	}
	else if (_selectionMode == CTRL_SELECT)
	{
		_forceLine = QLine();
		_selectionMode = NONE_SELECT;
		_selectedTetPointId = -1;
		_mouseForce.setZero();
	}
	else QGLViewer::mouseReleaseEvent(e);
}

void BaseScene::keyPressEvent(QKeyEvent * e)
{
	if (e->key() == Qt::Key_D)
		startAnimation();
}

void BaseScene::wheelEvent(QWheelEvent *e)
{
	QGLViewer::wheelEvent(e);
}

void BaseScene::drawWithNames()
{
	if (!_sim)
		return;
	if (_sim->getModelsNum() == 0)return;

	if (_status.enableSelectSurfacePoints)
	{
		for (int i = 0; i < _sim->totalPointsNum; i++)
		{
			glPushName(i);
			glBegin(GL_POINTS);
			qeal x, y, z;
			_sim->getSurfacePoint(i, x, y, z);
			glVertex3f(x, y, z);
			glEnd();
			glPopName();
		}
	}
	else //if (_status.enableSelectTetNodes)
	{
		for (int i = 0; i < _sim->totalTetPointsNum; i++)
		{
			glPushName(i);
			glBegin(GL_POINTS);
			qeal x, y, z;
			_sim->getTetPoint(i, x, y, z);
			glVertex3f(x, y, z);
			glEnd();
			glPopName();
		}
	}
}

void BaseScene::endSelection(const QPoint& point)
{
	update();
	if (_sim == nullptr)
		return;
	GLint nbHits = glRenderMode(GL_RENDER);
	if (nbHits <= 0)
		return;
	//std::cout << "Select : "<< nbHits << std::endl;
	for (int i = 0; i < nbHits; ++i)
	{
		int id = (selectBuffer())[4 * i + 3];
		if (_status.enableSelectSurfacePoints)
		{
			if (id >= _sim->totalPointsNum) continue;
			if (_selectionMode == ALT_SELECT)
			{
				_sim->selectSurfacePoints.insert(id);
			}
			else if (_selectionMode == ALT_REMOVE)
			{
				_sim->selectSurfacePoints.erase(id);
			}
		}
		else if (_status.enableSelectTetNodes)
		{			
			if (id >= _sim->totalTetPointsNum) continue;
			if (_selectionMode == ALT_SELECT)
			{
				_sim->selectTetPoints.insert(id);
				std::cout << id << std::endl;
			}
			else if (_selectionMode == ALT_REMOVE)
			{
				_sim->selectTetPoints.erase(id);
			}
		}
		else if(_selectionMode == CTRL_SELECT)
		{
			_selectedTetPointId = -1;
			_mouseForce.setZero();

			std::vector<Vector3> selectList;
			std::vector<uint32_t> select_index_list;
			Vector3 cent = Vector3(0, 0, 0);
			for (uint32_t i = 0; i < nbHits; i++)
			{
				int buf_id = selectBuffer()[4 * i + 3];
				if (buf_id < _sim->totalTetPointsNum)
				{
					qeal x, y, z;
					_sim->getTetPoint(buf_id, x, y, z);
					Vector3 pos = Vector3(x, y, z);
					selectList.push_back(pos);
					select_index_list.push_back(buf_id);
					cent += pos;
				}
			}

			cent /= selectList.size();
			qeal min_dist = QEAL_MAX;
			for (uint32_t i = 0; i < selectList.size(); i++)
			{
				qeal len = (selectList[i] - cent).norm();
				if (len < min_dist)
				{
					min_dist = len;
					_selectedTetPointId = select_index_list[i];
				}
			}
		}
	}

}

void BaseScene::initCamera()
{
	_camera = new StandardCamera();
	qglviewer::Camera *c = camera();
	setCamera(_camera);
	delete c;
	
	qglviewer::Vec world_origin = qglviewer::Vec(0.f, 0.f, 0.f);
	setSceneCenter(world_origin);

	camera()->setType(qglviewer::Camera::PERSPECTIVE);
	camera()->setPosition(qglviewer::Vec(0.0, 3.0, 10.0));
	camera()->lookAt(qglviewer::Vec(0.0, 0.0, 0.0));

	_lightCamera.setType(qglviewer::Camera::ORTHOGRAPHIC);
}

void BaseScene::initShaderProgram()
{
	_phongProgram.initProgram("Shader/Program/phong_shader_program.vs", "Shader/Program/phong_shader_program.fs");
	_shadowProgram.initProgram("Shader/Program/depth_map.vs", "Shader/Program/depth_map.fs");
	_textProgram.initProgram("Shader/Program/text_shader_program.vs", "Shader/Program/text_shader_program.fs");
	_textPainter.generateFont(QOpenGLContext::currentContext()->functions());
}

void BaseScene::initShadowFramebuffer()
{
	glGenFramebuffers(1, &_depthFramebuffer);
	glGenTextures(1, &_depthMapId);
	glBindTexture(GL_TEXTURE_2D, _depthMapId);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, 1024, 1024, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

	GLfloat borderColor[] = { 1.0, 1.0, 1.0, 1.0 };
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

	glBindFramebuffer(GL_FRAMEBUFFER, _depthFramebuffer);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, _depthMapId, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void BaseScene::initAnimation()
{
	stopAnimation();
	setAnimationPeriod(0);
	QDateTime current_date_time = QDateTime::currentDateTime();
	QString snapshotFileName = QString(_status.animationFilePath.c_str()) + QString("//") + current_date_time.toString("yyyy-MM-dd-hh-mm");
	int code = mkdir(snapshotFileName.toStdString().c_str());
	snapshotFileName += QString("/frame");
	setSnapshotFileName(snapshotFileName);
	setSnapshotFormat(QString("PNG"));
	setSnapshotQuality(100);
	connect(this, SIGNAL(saveAnimation(bool)), SLOT(saveSnapshot(bool)));
}

void BaseScene::renderDepthToBuffer()
{	
	if (!_lights.useShadow())
		return;
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_CULL_FACE);
	_shadowProgram.bind();
	QVector3D lightPos = _lights.getMainLightPosition();
	_lightCamera.setPosition(qglviewer::Vec(lightPos[0], lightPos[1], lightPos[2]));
	_lightCamera.lookAt(qglviewer::Vec(0.0, 0.0, 0.0));

	_shadowProgram.setUniformValue("lightSpaceMatrix", _lightCamera.getProjectionViewQMatrix());
	glBindFramebuffer(GL_FRAMEBUFFER, _depthFramebuffer);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, 1024, 1024);
	glCullFace(GL_FRONT);
	renderScene(&_shadowProgram);
	glCullFace(GL_BACK);
	_shadowProgram.release();
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	//glDisable (GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
}

void BaseScene::renderScene(QOpenGLShaderProgram * program)
{
	_floor->render(program, QOpenGLContext::currentContext()->functions());
	if (_sim != nullptr)
	{
		_sim->render(program, QOpenGLContext::currentContext()->functions());
		if (_status.renderLineMode)
		{
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			_sim->render(program, QOpenGLContext::currentContext()->functions(), true);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}
	}		
}

void BaseScene::renderScene()
{
	if (_sim != nullptr)
	{
		_sim->renderExtraElementOnCPU(QOpenGLContext::currentContext()->functions());
	}
}

void BaseScene::renderCornerAxis()
{

	int viewport[4];
	int scissor[4];

	// The viewport and the scissor are changed to fit the lower left
	// corner. Original values are saved.
	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetIntegerv(GL_SCISSOR_BOX, scissor);

	// Axis viewport size, in pixels
	const int size = 150;
	glViewport(0, 0, size, size);
	glScissor(0, 0, size, size);

	// The Z-buffer is cleared to make the axis appear over the
	// original image.
	glClear(GL_DEPTH_BUFFER_BIT);

	// Tune for best line rendering
	glDisable(GL_LIGHTING);
	glLineWidth(3.0);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -1, 1);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixd(camera()->orientation().inverse().matrix());

	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(1.0, 0.0, 0.0);

	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 1.0, 0.0);

	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 1.0);
	glEnd();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glEnable(GL_LIGHTING);

	// The viewport and the scissor are restored.
	glScissor(scissor[0], scissor[1], scissor[2], scissor[3]);
	glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
}

void BaseScene::renderInfoOnScene()
{
	double sw = camera()->screenWidth();
	double sh = camera()->screenHeight();
	QMatrix4x4 ortho_matrix = QMatrix4x4(2.0 / (sw), 0, 0, -1,
		0, 2.0 / (sh), 0, -1,
		0, 0, 2, -1,
		0, 0, 0, 1);

	_textProgram.bind();
	_textProgram.setUniformValue("mvp_matrix", ortho_matrix);

	const int w_step = 150;
	const int h_step = 20;

	int fps_w = sw - w_step;
	int fps_h = sh - h_step;

	QVector3D fontColor = QVector3D(0.3, 0.8, 0.3);
	qeal avg_cost = _status.animationCost / _status.animationLifeTime;
	if (avg_cost != 0.0 && _status.animationLifeTime > 0)
		_status.avg_fps = 1000.0 / (avg_cost);
	else _status.avg_fps = 1000.0;
	if (_status.avg_fps > 1000.0) _status.avg_fps = 1000.0;
	QString fps_str = QString::number(_status.avg_fps, 'f', 2);
	_textPainter.renderText(&_textProgram, "FPS: " + fps_str.toStdString() + " Hz", fps_w, fps_h, 0.5, fontColor);

	QString frame_str = QString::number(_status.animationLifeTime);
	_textPainter.renderText(&_textProgram, "Frames:   " + frame_str.toStdString(), fps_w, fps_h - h_step, 0.5, fontColor);

	_textProgram.release();
}

void BaseScene::renderSelectionRect()
{
	glBlendFunc(GL_ONE, GL_ONE);
	startScreenCoordinatesSystem();
	glDisable(GL_LIGHTING);
	glEnable(GL_BLEND);

	glColor4f(0.0, 0.0, 0.3f, 0.3f);
	glBegin(GL_QUADS);
	glVertex2i(_rectangle.left(), _rectangle.top());
	glVertex2i(_rectangle.right(), _rectangle.top());
	glVertex2i(_rectangle.right(), _rectangle.bottom());
	glVertex2i(_rectangle.left(), _rectangle.bottom());
	glEnd();

	glLineWidth(2.0);
	glColor4f(0.4f, 0.4f, 0.5f, 0.5f);
	glBegin(GL_LINE_LOOP);
	glVertex2i(_rectangle.left(), _rectangle.top());
	glVertex2i(_rectangle.right(), _rectangle.top());
	glVertex2i(_rectangle.right(), _rectangle.bottom());
	glVertex2i(_rectangle.left(), _rectangle.bottom());
	glEnd();

	glDisable(GL_BLEND);
	glEnable(GL_LIGHTING);
	stopScreenCoordinatesSystem();
	update();
}

void BaseScene::renderSelectionPoints()
{
	if (_sim == nullptr)
		return;
	glPointSize(3.0);
	if (_status.enableSelectSurfacePoints)
	{
		for (auto it = _sim->selectSurfacePoints.begin(); it != _sim->selectSurfacePoints.end(); it++)
		{
			int i = *it;
			glColor3f(0.8, 0.2, 0.2);
			glBegin(GL_POINTS);
			qeal x, y, z;
			_sim->getSurfacePoint(i, x, y, z);
			glVertex3f(x, y, z);
			glEnd();
			glPopName();
		}
	}
	else if (_status.enableSelectTetNodes)
	{
		for (auto it = _sim->selectTetPoints.begin(); it != _sim->selectTetPoints.end(); it++)
		{
			int i = *it;
			glColor3f(0.2, 0.8, 0.2);
			glBegin(GL_POINTS);
			qeal x, y, z;
			_sim->getTetPoint(i, x, y, z);
			glVertex3f(x, y, z);
			glEnd();
			glPopName();
		}
	}
}

void BaseScene::renderForceLine()
{
	if (_selectedTetPointId == -1 || _selectionMode != CTRL_SELECT || _sim == nullptr)
		return;
	if (_mouseForce.norm() < 1e-8)
		return;
	qeal x, y, z;
	_sim->getTetPoint(_selectedTetPointId, x, y, z);
	Vector3 selectedPos = Vector3(x, y, z);
	qglviewer::Vec pro_pos = camera()->projectedCoordinatesOf(qglviewer::Vec(x, y, z));

	qglviewer::Vec norm = qglviewer::Vec(_mouseForce.data()[0], _mouseForce.data()[1], _mouseForce.data()[2]);
	norm.normalize();
	qglviewer::Vec epos = qglviewer::Vec(selectedPos.data()[0], selectedPos.data()[1], selectedPos.data()[2]) + norm * 0.4;
	if (epos.y < 0.00001) epos.y = 0.00001;
	qglviewer::Vec spos = qglviewer::Vec(selectedPos.data()[0], selectedPos.data()[1], selectedPos.data()[2]);

	glBegin(GL_LINES);
	glColor3f(240.0 / 255.0, 240.0 / 255.0, 240.0 / 255.0);
	glVertex3d(spos.x, spos.y, spos.z);
	glColor3f(240.0 / 255.0, 240.0 / 255.0, 240.0 / 255.0);
	glVertex3d(epos.x, epos.y, epos.z);
	glEnd();

	double r = 0.03;

	glPushMatrix();
	glTranslated(spos.x, spos.y, spos.z);
	glColor3f(0.8, 0.2, 0.2);
	gluSphere(gluNewQuadric(), r, 50, 50);
	glPopMatrix();
	update();

	glPushMatrix();
	glTranslated(epos.x, epos.y, epos.z);
	glColor3f(0.8, 0.2, 0.2);
	gluSphere(gluNewQuadric(), r, 50, 50);
	glPopMatrix();
}

void BaseScene::computeMouseForce()
{
	if (_selectedTetPointId == -1 || _sim == nullptr)
	{
		_mouseForce = Vector3(0, 0, 0);
		return;
	}

	if (_forceLine.x1() == _forceLine.x2() && _forceLine.y1() == _forceLine.y2())
	{
		_mouseForce = Vector3(0, 0, 0);
		return;
	}

	qeal x, y, z;
	_sim->getTetPoint(_selectedTetPointId, x, y, z);
	Vector3 selectedPos = Vector3(x, y, z);
	qglviewer::Vec pro_pos = camera()->projectedCoordinatesOf(qglviewer::Vec(x, y, z));

	qglviewer::Vec endPos = qglviewer::Vec(_forceLine.x2(), _forceLine.y2(), pro_pos.z);
	qglviewer::Vec spos = qglviewer::Vec(x, y, z);
	qglviewer::Vec epos = camera()->unprojectedCoordinatesOf(endPos);

	qglviewer::Vec dir = epos - spos;
	_mouseForce = Vector3(dir.x, dir.y, dir.z);

	qeal scale = 50.0;

	_mouseForce *= scale;
	
	getCurrentSimulator()->handleMouseForce(_selectedTetPointId, _mouseForce.data()[0], _mouseForce.data()[1], _mouseForce.data()[2]);
}

void BaseScene::removeTet()
{
	BaseTetElementSet ele_set = _sim->models[0]->getTetMeshHandle()->_elementSet[select_tet_group];
	std::vector<int> elelist;
	std::set<int> list;
	ele_set.getElements(list);
	std::set<int>::iterator cit = list.begin();
	for (; cit != list.end(); ++cit)
	{
		elelist.push_back(*cit);
	}
	
	std::set<int>::iterator it = select_tet_element.begin();
	for (; it != select_tet_element.end(); ++it)
	{		
		elelist[*it] = -1;
	}

	ele_set.clear();
	for (int i = 0; i < elelist.size(); i++)
	{
		if (elelist[i] != -1)
		{
			ele_set.insert(elelist[i]);
		}
	}

	_sim->models[0]->getTetMeshHandle()->_elementSet[select_tet_group] = ele_set;
	removeSelect();
}

void BaseScene::removeSelect()
{
	select_tet_element.clear();
}
