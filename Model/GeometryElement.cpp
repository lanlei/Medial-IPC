#include "GeometryElement.h"

Box::Box(qeal halfX, qeal halfY, qeal halfZ, QVector3D center, QColor uniformColor):
	_halfX(halfX),
	_halfY(halfY),
	_halfZ(halfZ),
	_center(center)
{
	renderGroupNum = 1;
	pointsNum = 8;
	facesNum = 12;
	renderGroupFaceNum.buffer = (int*)malloc(sizeof(int));
	*(renderGroupFaceNum.buffer) = facesNum;
	renderGroupFaceNum.offset = 0;
	renderGroupFaceNum.span = 1;
	
	QVector3D min = center - QVector3D(_halfX, _halfY, _halfZ);
	QVector3D max = center + QVector3D(_halfX, _halfY, _halfZ);

	setFromBounding(min, max);
	texCoords.buffer = (qeal*)malloc(2 * 8 * sizeof(qeal));
	texCoords.offset = 0;
	texCoords.span = 2 * 8;
	//

	qeal buffer[16] = { 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1 };
	std::copy(buffer, buffer + 16, texCoords.buffer);

	faceIndices.buffer = (int*)malloc(3 * 12 * sizeof(int));
	faceIndices.offset = 0;
	faceIndices.span = 36;

	int indices[36] = {
	0, 1, 2,
	0, 2, 3,
	3, 2, 6,
	3, 6, 7,
	7, 6, 5,
	7, 5, 4,
	4, 5, 1,
	4, 1, 0,
	2, 1, 5,
	2, 5, 6,
	7, 4, 0,
	3, 7, 0 };
	std::copy(indices, indices + 36, faceIndices.buffer);
	renderMaterials.resize(1);

	renderMaterials[0] = new BaseRenderMaterial();
	renderMaterials[0]->ambient = uniformColor;
	renderMaterials[0]->diffuse = uniformColor;
	renderMaterials[0]->specular = uniformColor;
	initVBO();
}

void Box::scale(qeal sx, qeal sy, qeal sz)
{
	_halfX *= sx;
	_halfY *= sy;
	_halfZ *= sz;
	QVector3D min = _center - QVector3D(_halfX, _halfY, _halfZ);
	QVector3D max = _center + QVector3D(_halfX, _halfY, _halfZ);
	setFromBounding(min, max);	
	updateVBO();
}

void Box::setFromBounding(const QVector3D min, const QVector3D max)
{
	//points.resize(24);
	if (points.buffer != nullptr)
		free(points.buffer);
	points.buffer = (qeal*)malloc(2 * 12 * sizeof(qeal));
	points.offset = 0;
	points.span = 24;

	points.buffer[0] = min[0];
	points.buffer[1] = min[1];
	points.buffer[2] = min[2];

	points.buffer[3] = min[0];
	points.buffer[4] = max[1];
	points.buffer[5] = min[2];

	points.buffer[6] = max[0];
	points.buffer[7] = max[1];
	points.buffer[8] = min[2];

	points.buffer[9] = max[0];
	points.buffer[10] = min[1];
	points.buffer[11] = min[2];

	points.buffer[12] = min[0];
	points.buffer[13] = min[1];
	points.buffer[14] = max[2];

	points.buffer[15] = min[0];
	points.buffer[16] = max[1];
	points.buffer[17] = max[2];

	points.buffer[18] = max[0];
	points.buffer[19] = max[1];
	points.buffer[20] = max[2];

	points.buffer[21] = max[0];
	points.buffer[22] = min[1];
	points.buffer[23] = max[2];

	_halfX = (max[0] - min[0]) / 2;
	_halfY = (max[1] - min[1]) / 2;
	_halfZ = (max[2] - min[2]) / 2;

	_center = (max + min) / 2.0;
	updateVBO();
}

void Box::setCenter(QVector3D center)
{
	QVector3D min = center - QVector3D(_halfX, _halfY, _halfZ);
	QVector3D max = center + QVector3D(_halfX, _halfY, _halfZ);
	setFromBounding(min, max);
}

Floor::Floor(int xz, QColor uniformColor)
{
	renderGroupNum = 1;
	pointsNum = 8;
	facesNum = 12;
	
	renderGroupFaceNum.buffer = (int*)malloc(sizeof(int));
	*(renderGroupFaceNum.buffer) = facesNum;
	renderGroupFaceNum.offset = 0;
	renderGroupFaceNum.span = 1;

	hasTextureCoords = true;
	scaleXZ(xz);
	_y = 0;

	faceIndices.buffer = (int*)malloc(3 * 12 * sizeof(int));
	faceIndices.offset = 0;
	faceIndices.span = 36;

	int indices[36] = {
	0, 1, 2,
	0, 2, 3,
	3, 2, 6,
	3, 6, 7,
	7, 6, 5,
	7, 5, 4,
	4, 5, 1,
	4, 1, 0,
	2, 1, 5,
	2, 5, 6,
	7, 4, 0,
	3, 7, 0 };
	std::copy(indices, indices + 36, faceIndices.buffer);
	renderMaterials.resize(1);
	renderMaterials[0] = new BaseRenderMaterial();
	renderMaterials[0]->ambient = uniformColor;
	renderMaterials[0]->readTextureMap(QString("./texture/floor_texture.jpg"), AmbientMapIndex);
	initVBO();
}

void Floor::scaleXZ(int xz)
{
	_scaleXZ = xz;

	texCoords.buffer = (qeal*)malloc(2 * 8 * sizeof(qeal));
	texCoords.offset = 0;
	texCoords.span = 2 * 8;

	points.buffer = (qeal*)malloc(2 * 12 * sizeof(qeal));
	points.offset = 0;
	points.span = 24;

	QVector3D min(-0.5 * _scaleXZ, -1.0, -0.5 * _scaleXZ);
	QVector3D max(0.5 * _scaleXZ, 0, 0.5 * _scaleXZ);
	points.buffer[0] = min[0];
	points.buffer[1] = min[1];
	points.buffer[2] = min[2];

	points.buffer[3] = min[0];
	points.buffer[4] = max[1];
	points.buffer[5] = min[2];

	points.buffer[6] = max[0];
	points.buffer[7] = max[1];
	points.buffer[8] = min[2];

	points.buffer[9] = max[0];
	points.buffer[10] = min[1];
	points.buffer[11] = min[2];

	points.buffer[12] = min[0];
	points.buffer[13] = min[1];
	points.buffer[14] = max[2];

	points.buffer[15] = min[0];
	points.buffer[16] = max[1];
	points.buffer[17] = max[2];

	points.buffer[18] = max[0];
	points.buffer[19] = max[1];
	points.buffer[20] = max[2];

	points.buffer[21] = max[0];
	points.buffer[22] = min[1];
	points.buffer[23] = max[2];
	//
	qeal buffer[16] = { 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1 };
	std::copy(buffer, buffer + 16, texCoords.buffer);
	for (int i = 0; i < texCoords.span; i++)
		texCoords.buffer[i] *= _scaleXZ;
}

void Floor::setYPlane(qeal Y)
{
	_y = Y;
	for (int i = 0; i < pointsNum; i++)
		points.buffer[3 * i + 1] = _y;
}

void SphereElement::operator=(const SphereElement& b)
{
	center = b.center;
	radius = b.radius;
	frame.setReferenceFrame(b.frame.referenceFrame());
}

void SphereElement::show()
{
	glPushMatrix();
	glTranslated(center.data()[0], center.data()[1], center.data()[2]);
	glColor3d(204.0 / 255.0, 204.0 / 255.0, 153.0 / 255.0);
	gluSphere(gluNewQuadric(), radius, 50, 50);
	glPopMatrix();
}

SplintElement::SplintElement(const Vector3 v0, const Vector3 v1, const Vector3 v2, Vector3 n)
{
	vt[0] = v0;
	vt[1] = v1;
	vt[2] = v2;
	if (n == Vector3(0, 0, 0))
		updateNormal();
}

void SplintElement::updateNormal(bool reverse)
{
	const Vector3 v0v1 = vt[1] - vt[0];
	const Vector3 v0v2 = vt[2] - vt[0];
	nt = v0v1.cross(v0v2);
	nt.normalize();
	if (reverse)
	{
		nt = nt * -1;
	}
}

void ConeElement::computeCone()
{
	Vector3 c0c1 = _sp1 - _sp0;
	// one sphere is included in another sphere
	if (c0c1.norm() - abs(_r1 - _r0) < 1e-8)
	{
		apex = _r1 > _r0 ? _sp1 : _sp0;
		axis = Vector3(0, 0, 1);
		base = _r0 > _r1 ? _r0 : _r1;
		top = _r1 < _r0 ? _r1 : _r0;
		height = 0.0;
		rot_axis = Vector3(0, 0, 1);
		rot_angle = 0.;
		return;
	}

	if (c0c1.norm() < 1e-8)
	{
		apex = _sp0;
		axis = Vector3(0, 0, 1);
		base = _r0;
		top = _r0;
		height = 0.;
		rot_axis = Vector3(0, 0, 1);
		rot_angle = 0.;
	}
	else
	{
		qeal dr0r1 = fabs(_r0 - _r1);
		if (dr0r1 < 1e-8)
		{
			apex = _sp0;
			axis = _sp1 - _sp0;
			axis.normalize();
			base = _r0;
			top = _r0;
			height = (_sp1 - _sp0).norm();
		}
		else
		{
			apex = (_r1 * _sp0 - _r0 * _sp1) / (_r1 - _r0);
			axis = (_r0 > _r1) ? (_sp1 - _sp0) : (_sp0 - _sp1);
			axis.normalize();

			qeal cangle;
			Vector3 apexc0 = apex - _sp0;
			qeal vc0len = apexc0.norm();
			Vector3 apexc1 = apex - _sp1;
			qeal vc1len = apexc1.norm();

			cangle = sqrt(1. - _r0 * _r0 / vc0len / vc0len);

			if (_r0 < _r1)
			{
				apex = _sp0 + apexc0 * cangle * cangle;
				base = _r1 * cangle;
				top = _r0 * cangle;
				height = abs(vc1len - vc0len) * cangle * cangle;
			}
			else
			{
				apex = _sp1 + apexc1 * cangle * cangle;
				base = _r0 * cangle;
				top = _r1 * cangle;
				height = abs(vc0len - vc1len) * cangle * cangle;
			}
		}

		Vector3 za(0, 0, 1);
		rot_angle = acos(axis.dot(za));
		if ((fabs(rot_angle) < 1e-12) || (fabs(rot_angle - M_PI) < 1e-12))
			rot_axis = Vector3(0, 0, 1);
		else
			rot_axis = axis.cross(za);
		rot_axis.normalize();
		rot_angle *= (180. / M_PI);
	}
}

void ConeElement::show()
{
	double br = base;

	Vector3 c = _sp0;
	qeal r = _r0;
	if (r < _r1)
	{
		c = _sp1;
		r = _r1;
	}

	double t = sqrt(r* r - br * br);
	Vector3 translation = c + t * axis;
	Vector3 cone_axis = axis;
	qglviewer::Vec y_axis(0, 0, 1);

	qglviewer::Quaternion quat(y_axis, qglviewer::Vec(cone_axis.data()));

	BaseFrame frame;
	frame.setPosition(translation[0], translation[1], translation[2]);
	frame.rotate(quat);

	glPushMatrix();
	glMultMatrixd(frame.worldMatrix());
	glColor3d(0.16, 0.16, 0.86);
	gluCylinder(gluNewQuadric(), base, top, height, 200, 1);
	glPopMatrix();
}

void TwoSplintElements::show()
{
	glPushMatrix();
	glBegin(GL_POLYGON);
	glColor3d(204.0 / 255.0, 204.0 / 255.0, 153.0 / 255.0);
	glNormal3f(st1.nt[0], st1.nt[1], st1.nt[2]);
	glVertex3f(st1.vt[0][0], st1.vt[0][1], st1.vt[0][2]);
	glVertex3f(st1.vt[1][0], st1.vt[1][1], st1.vt[1][2]);
	glVertex3f(st1.vt[2][0], st1.vt[2][1], st1.vt[2][2]);
	glEnd();

	glBegin(GL_POLYGON);
	glColor3d(204.0 / 255.0, 204.0 / 255.0, 153.0 / 255.0);
	glNormal3f(-st1.nt[0], -st1.nt[1], -st1.nt[2]);
	glVertex3f(st1.vt[0][0], st1.vt[0][1], st1.vt[0][2]);
	glVertex3f(st1.vt[2][0], st1.vt[2][1], st1.vt[2][2]);
	glVertex3f(st1.vt[1][0], st1.vt[1][1], st1.vt[1][2]);
	glEnd();

	glBegin(GL_POLYGON);
	glColor3d(204.0 / 255.0, 204.0 / 255.0, 153.0 / 255.0);
	glNormal3f(st2.nt[0], st2.nt[1], st2.nt[2]);
	glVertex3f(st2.vt[0][0], st2.vt[0][1], st2.vt[0][2]);
	glVertex3f(st2.vt[1][0], st2.vt[1][1], st2.vt[1][2]);
	glVertex3f(st2.vt[2][0], st2.vt[2][1], st2.vt[2][2]);
	glEnd();

	glBegin(GL_POLYGON);
	glColor3d(204.0 / 255.0, 204.0 / 255.0, 153.0 / 255.0);
	glNormal3f(-st2.nt[0], -st2.nt[1], -st2.nt[2]);
	glVertex3f(st2.vt[0][0], st2.vt[0][1], st2.vt[0][2]);
	glVertex3f(st2.vt[2][0], st2.vt[2][1], st2.vt[2][2]);
	glVertex3f(st2.vt[1][0], st2.vt[1][1], st2.vt[1][2]);
	glEnd();

	glPopMatrix();
}
