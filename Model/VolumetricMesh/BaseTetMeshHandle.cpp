#include "BaseTetMeshHandle.h"

BaseTetMeshHandle::BaseTetMeshHandle(BaseTetMesh * tetMesh):_tetMesh(tetMesh), renderElementSetId(0)
{
}

BaseTetMesh * BaseTetMeshHandle::getTetMesh() const
{
		return _tetMesh;
}

void BaseTetMeshHandle::init()
{
	if (!_tetMesh)
		return;

	_elementsParam.resize(_tetMesh->tetElementNum);
	for (size_t eleid = 0; eleid < _tetMesh->tetElementNum; eleid++)
		_elementsParam[eleid] = new BaseTetElementParam(eleid);

	_nodesParam.resize(_tetMesh->tetPointsNum);
	for (size_t nid = 0; nid < _tetMesh->tetPointsNum; nid++)
		_nodesParam[nid] = new BaseTetNodeParam();

	_materials.resize(1);
	_materials[0] = new ENuMaterial("defalut", 1.0, 100, 0.4);

	_elementMaterials.resize(_tetMesh->tetElementNum);
	std::fill(_elementMaterials.begin(), _elementMaterials.end(), 0);
	_nodeMaterials.resize(_tetMesh->tetPointsNum);
	std::fill(_nodeMaterials.begin(), _nodeMaterials.end(), 0);

	BaseTetElementSet elementSet("defalut");
	for (size_t eleid = 0; eleid < _tetMesh->tetElementNum; eleid++)
		elementSet.insert(eleid);
	_elementSet.push_back(elementSet);

	BaseTetNodeSet nodeSet("defalut");
	for (size_t nid = 0; nid < _tetMesh->tetPointsNum; nid++)
		nodeSet.insert(nid);
	_nodeSet.push_back(nodeSet);

	correctElementNodeOrder();
	computeVolume();
	computeElementFaceArea();
	computeTetMeshNormal();
	computeElementShapeFunctionDerivate();
	initElementABCD();
	markBoundaryElement();
}

void BaseTetMeshHandle::correctElementNodeOrder()
{
	for (size_t eleid = 0; eleid < _tetMesh->tetElementNum; eleid++)
	{
		BaseTetElementParam* ele = getTetElementParam(eleid);
		Vector4i indices = _tetMesh->getTetElement(eleid);
		int n0, n1, n2, n3;
		n0 = indices.data()[0];
		n1 = indices.data()[1];
		n2 = indices.data()[2];
		n3 = indices.data()[3];

		Vector3 direction_v;

		TriangleNormal(_tetMesh->getTetPoint(n0), _tetMesh->getTetPoint(n1), _tetMesh->getTetPoint(n2), direction_v);

		Vector3 n3n0 = _tetMesh->getTetPoint(n3) - _tetMesh->getTetPoint(n0);
		if (n3n0.dot(direction_v) < 0)
		{
			indices.data()[1] = n2;
			indices.data()[2] = n1;
		}

		for (int i = 0; i < 4; i++)
		{
			int n0 = indices.data()[ele->faceIndices[i * 3 + 0]];
			int n1 = indices.data()[ele->faceIndices[i * 3 + 1]];
			int n2 = indices.data()[ele->faceIndices[i * 3 + 2]];

			TriangleNormal(_tetMesh->getTetPoint(n0), _tetMesh->getTetPoint(n1), _tetMesh->getTetPoint(n2), ele->facesNormal[i]);

			int n3;
			for (int j = 0; j < 4; j++)
			{
				if (!(indices.data()[j] == n0 || indices.data()[j] == n1 || indices.data()[j] == n2))
				{
					n3 = indices.data()[j];
					break;
				}
			}

			Vector3 n3n0 = _tetMesh->getTetPoint(n3) - _tetMesh->getTetPoint(n0);
			if (n3n0.dot(ele->facesNormal[i]) > 0)
			{
				ele->facesNormal[i] *= -1;
				int temp = ele->faceIndices[i * 3 + 0];
				ele->faceIndices[i * 3 + 0] = ele->faceIndices[i * 3 + 2];
				ele->faceIndices[i * 3 + 2] = temp;
			}
		}

		_tetMesh->setTetElement(eleid, indices.data());
	}
}

void BaseTetMeshHandle::computeVolume()
{
	for (size_t i = 0; i < _tetMesh->tetElementNum; i++)
	{
		VectorX points = _tetMesh->getTetElementPoint(i);
		getTetElementParam(i)->computeVolume(points.data());
	}
	for (size_t i = 0; i < _tetMesh->tetPointsNum; i++)
	{
		qeal v = 0;
		for (size_t j = 0; j < _tetMesh->tetPointElementList[i].span; j++)
		{
			int eleid = _tetMesh->tetPointElementList[i].buffer[j];
			v += getTetElementParam(eleid)->volume / 4.0;
		}
		getTetNodeParam(i)->volume = v;
	}
}

void BaseTetMeshHandle::computeElementFaceArea()
{
	for (size_t i = 0; i < _tetMesh->tetElementNum; i++)
	{
		VectorX points = _tetMesh->getTetElementPoint(i);
		getTetElementParam(i)->computeFaceArea(points.data());
	}
}

void BaseTetMeshHandle::computeTetMeshNormal()
{
	// element normal
	for (size_t i = 0; i < _tetMesh->tetElementNum; i++)
	{
		VectorX points = _tetMesh->getTetElementPoint(i);
		getTetElementParam(i)->computeElementNormal(points.data());
	}
}

void BaseTetMeshHandle::initElementABCD()
{
	for (size_t eleid = 0; eleid < _tetMesh->tetElementNum; eleid++)
	{
		getTetElementParam(eleid)->computeABCD();
	}
}

void BaseTetMeshHandle::computeElementShapeFunctionDerivate()
{
	for (size_t eleid = 0; eleid < _tetMesh->tetElementNum; eleid++)
	{
		VectorX points = _tetMesh->getTetElementPoint(eleid);
		getTetElementParam(eleid)->computeShapeFunctionDerivate(points.data());
	}
}

void BaseTetMeshHandle::markBoundaryElement()
{
	_elementNeighbor.resize(_tetMesh->tetElementNum);
	_boundaryElementFlag.resize(_tetMesh->tetElementNum);
	for (int eleId = 0; eleId < _tetMesh->tetElementNum; eleId++)
	{
		_elementNeighbor[eleId].resize(4);
		BaseTetElementParam* eleParam = getTetElementParam(eleId);
		Vector4i ele = _tetMesh->getTetElement(eleId);
		for (int f = 0; f < 4; f++)
		{
			_elementNeighbor[eleId][f] = -1;
			int n0 = ele.data()[eleParam->faceIndices[f * 3 + 0]];
			int n1 = ele.data()[eleParam->faceIndices[f * 3 + 1]];
			int n2 = ele.data()[eleParam->faceIndices[f * 3 + 2]];

			for (int j = 0; j < _tetMesh->tetPointElementList[n0].span; j++)
			{
				int neleId = _tetMesh->tetPointElementList[n0].buffer[j];
				if (neleId == eleId)continue;
				int n1FindSame = false;
				for (int k = 0; k < _tetMesh->tetPointElementList[n1].span; k++)
				{
					int n1_neleId = _tetMesh->tetPointElementList[n1].buffer[k];
					if (n1_neleId == neleId)
					{
						n1FindSame = true;
						break;
					}
				}

				int n2FindSame = false;
				for (int k = 0; k < _tetMesh->tetPointElementList[n2].span; k++)
				{
					int n2_neleId = _tetMesh->tetPointElementList[n2].buffer[k];
					if (n2_neleId == neleId)
					{
						n2FindSame = true;
						break;
					}
				}

				if (n1FindSame && n2FindSame)
				{
					_elementNeighbor[eleId][f] = neleId;
					break;
				}
			}

			if (_elementNeighbor[eleId][f] == -1)
				_boundaryElementFlag[eleId] = true;
		}
	}
}

BaseTetElementParam::BaseTetElementParam(int eleid)
{
	this->eleid = eleid;
	faceIndices[0] = 3;
	faceIndices[1] = 1;
	faceIndices[2] = 0;

	faceIndices[3] = 3;
	faceIndices[4] = 2;
	faceIndices[5] = 1;

	faceIndices[6] = 3;
	faceIndices[7] = 0;
	faceIndices[8] = 2;

	faceIndices[9] = 2;
	faceIndices[10] = 1;
	faceIndices[11] = 0;
}

void BaseTetElementParam::computeElementDeformationByShapeMatrix(Matrix3 & Ds, qeal * displacement, int * indices)
{
	Ds.setZero();

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			Ds.data()[j * 3 + i] = displacement[indices[j] * 3 + i] - displacement[indices[3] * 3 + i];
		}
	}

	Ds += this->Dm;
}

void BaseTetElementParam::computeShapeFunctionDerivate(qeal* points)
{
	Vector4 phi;
	Matrix4 referenceShape;

	for (int i = 0; i < 4; i++)
	{
		Vector3 p = Vector3(points[3 * i], points[3 * i + 1], points[3 * i + 2]);
		for (int j = 0; j < 3; j++)
			referenceShape.data()[i * 4 + j] = p.data()[j];
		referenceShape.data()[i * 4 + 3] = 1;
	}
	invShapeFunction = referenceShape.inverse();

	phiDerivate.resize(4, 3);

	for (int i = 0; i < 4; i++)
	{
		phiGradient[i].data()[0] = invShapeFunction.data()[0 * 4 + i];
		phiGradient[i].data()[1] = invShapeFunction.data()[1 * 4 + i];
		phiGradient[i].data()[2] = invShapeFunction.data()[2 * 4 + i];

		phiDerivate.data()[0 * 4 + i] = phiGradient[i].data()[0];
		phiDerivate.data()[1 * 4 + i] = phiGradient[i].data()[1];
		phiDerivate.data()[2 * 4 + i] = phiGradient[i].data()[2];
	}
	//init Dm Dm_inverse

	// FEM simulation of 3D deformable solids: a practitioner's guide to theory, discretization and model reduction

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			Dm.data()[j * 3 + i] = points[3 * j + i] - points[3 * 3 + i];
	}
	invDm = Dm.inverse();
	invDmT = invDm.transpose();

	dFdu.resize(9, 12);
	dFdu.setZero();
	qeal Dm11 = invDm.data()[0];
	qeal Dm21 = invDm.data()[1];
	qeal Dm31 = invDm.data()[2];
	qeal Dm12 = invDm.data()[3];
	qeal Dm22 = invDm.data()[4];
	qeal Dm32 = invDm.data()[5];
	qeal Dm13 = invDm.data()[6];
	qeal Dm23 = invDm.data()[7];
	qeal Dm33 = invDm.data()[8];

	dFdu.data()[0] = Dm11;    dFdu.data()[10] = Dm11;    dFdu.data()[20] = Dm11;
	dFdu.data()[3] = Dm12;    dFdu.data()[13] = Dm12;    dFdu.data()[23] = Dm12;
	dFdu.data()[6] = Dm13;    dFdu.data()[16] = Dm13;    dFdu.data()[26] = Dm13;

	dFdu.data()[27] = Dm21;    dFdu.data()[37] = Dm21;    dFdu.data()[47] = Dm21;
	dFdu.data()[30] = Dm22;    dFdu.data()[40] = Dm22;    dFdu.data()[50] = Dm22;
	dFdu.data()[33] = Dm23;    dFdu.data()[43] = Dm23;    dFdu.data()[53] = Dm23;

	dFdu.data()[54] = Dm31;    dFdu.data()[64] = Dm31;    dFdu.data()[74] = Dm31;
	dFdu.data()[57] = Dm32;    dFdu.data()[67] = Dm32;    dFdu.data()[77] = Dm32;
	dFdu.data()[60] = Dm33;    dFdu.data()[70] = Dm33;    dFdu.data()[80] = Dm33;

	dFdu.data()[81] = -1 * (Dm11 + Dm21 + Dm31);    dFdu.data()[91] = -1 * (Dm11 + Dm21 + Dm31);    dFdu.data()[101] = -1 * (Dm11 + Dm21 + Dm31);
	dFdu.data()[84] = -1 * (Dm12 + Dm22 + Dm32);   dFdu.data()[94] = -1 * (Dm12 + Dm22 + Dm32);    dFdu.data()[104] = -1 * (Dm12 + Dm22 + Dm32);
	dFdu.data()[87] = -1 * (Dm13 + Dm23 + Dm33);    dFdu.data()[97] = -1 * (Dm13 + Dm23 + Dm33);    dFdu.data()[107] = -1 * (Dm13 + Dm23 + Dm33);
}

void BaseTetElementParam::computeVolume(qeal * points)
{
	Matrix3 D;
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 3; j++)
			D.data()[3 * i + j] = points[3 * (i + 1) + j] - points[3 * 0 + j];
	volume = abs(D.determinant() / 6.0);
}

void BaseTetElementParam::computeCenter(qeal * points)
{
	center.setZero();
	for (size_t i = 0; i < 3; i++)
		for (size_t j = 0; j < 4; j++)
			center.data()[i] += points[3 * j + i];
	center /= 4.0;
}

void BaseTetElementParam::computeFaceArea(qeal * points)
{
	for (size_t j = 0; j < 4; j++)
	{
		int n0, n1, n2;
		n0 = faceIndices[3 * j + 0];
		n1 = faceIndices[3 * j + 1];
		n2 = faceIndices[3 * j + 2];
		Vector3 p0 = Vector3(points[3 * n0], points[3 * n0 + 1], points[3 * n0 + 2]);
		Vector3 p1 = Vector3(points[3 * n1], points[3 * n1 + 1], points[3 * n1 + 2]);
		Vector3 p2 = Vector3(points[3 * n2], points[3 * n2 + 1], points[3 * n2 + 2]);
		area[j] = TriangleArea(p0, p1, p2);
	}
}

void BaseTetElementParam::computeElementNormal(qeal* points)
{
	for (int j = 0; j < 4; j++)
	{
		int n0, n1, n2;
		n0 = faceIndices[3 * j + 0];
		n1 = faceIndices[3 * j + 1];
		n2 = faceIndices[3 * j + 2];
		Vector3 p0 = Vector3(points[3 * n0], points[3 * n0 + 1], points[3 * n0 + 2]);
		Vector3 p1 = Vector3(points[3 * n1], points[3 * n1 + 1], points[3 * n1 + 2]);
		Vector3 p2 = Vector3(points[3 * n2], points[3 * n2 + 1], points[3 * n2 + 2]);

		TriangleNormal(p0, p1, p2, facesNormal[j]);
	}

	for (int n = 0; n < 4; n++)
	{
		nodesNormal[n].setZero();
		for (int f = 0; f < 4; f++)
		{
			for (int fn = 0; fn < 3; fn++)
			{
				if (faceIndices[f * 3 + fn] == n)
				{
					nodesNormal[n] += area[f] * facesNormal[f];
				}
			}
		}
		nodesNormal[n] /= 3.0;
	}
}

void BaseTetElementParam::computeABCD()
{
	//init ABCD
	//init A
	for (int a = 0; a < 4; a++)
	{
		for (int b = 0; b < 4; b++)
		{
			A[a][b] = phiGradient[a] * phiGradient[b].transpose()*volume;
		}
	}
	//init B
	for (int a = 0; a < 4; a++)
	{
		for (int b = 0; b < 4; b++)
		{
			B[a][b] = phiGradient[a].dot(phiGradient[b])*volume;
		}
	}
	//init C
	for (int a = 0; a < 4; a++)
	{
		for (int b = 0; b < 4; b++)
		{
			for (int c = 0; c < 4; c++)
			{
				C[a][b][c] = phiGradient[a] * (phiGradient[b].dot(phiGradient[c]))*volume;
			}
		}
	}

	//init D
	for (int a = 0; a < 4; a++)
	{
		for (int b = 0; b < 4; b++)
		{
			for (int c = 0; c < 4; c++)
			{
				for (int d = 0; d < 4; d++)
				{
					D[a][b][c][d] = phiGradient[a].dot(phiGradient[b]) * phiGradient[c].dot(phiGradient[d])*volume;
				}
			}
		}
	}
}

BaseTetElementParam * BaseTetMeshHandle::getTetElementParam(const int eid)
{
	return _elementsParam[eid];
}

BaseTetNodeParam * BaseTetMeshHandle::getTetNodeParam(const int nid)
{
	return _nodesParam[nid];
}

BaseTetMeshMaterial * BaseTetMeshHandle::getElementMaterialById(int id) const
{
	return _materials[id];
}

BaseTetMeshMaterial * BaseTetMeshHandle::getElementMaterial(int eid) const
{
	return _materials[_elementMaterials[eid]];
}

BaseTetMeshMaterial * BaseTetMeshHandle::getNodeMaterial(int nid) const
{
	return _materials[_nodeMaterials[nid]];
}

bool BaseTetMeshHandle::readTetMeshSetInfo(std::string& filename)
{
	std::ifstream fin(filename.c_str());
	if (!fin.is_open())
		return false;

	std::vector<int> elementFlag(_tetMesh->tetElementNum, -1);
	int setNum = 0;
	fin >> setNum;
	if (setNum == 0) return true;
	_elementSet.clear();
	_nodeSet.clear();
	for (int i = 0; i < setNum; i++)
	{
		std::string name;
		int eleNum;	
		fin >> name >> eleNum;
		std::cout << name << " " << eleNum << std::endl;
		BaseTetElementSet elementSet(name);
		BaseTetNodeSet pointSet(name);
		for (int j = 0; j < eleNum; j++)
		{
			int eid;
			fin >> eid;
			elementSet.insert(eid);
			elementFlag[eid] = i;

			Vector4i ele = _tetMesh->getTetElement(eid);
			pointSet.insert(ele[0]);
			pointSet.insert(ele[1]);
			pointSet.insert(ele[2]);
			pointSet.insert(ele[3]);
		}
		_elementSet.push_back(elementSet);
		_nodeSet.push_back(pointSet);
	}
	for (int i = 0; i < _tetMesh->tetElementNum; i++)
	{
		if (elementFlag[i] == -1)
		{
			_elementSet[0].insert(i);
			Vector4i ele = _tetMesh->getTetElement(i);
			_nodeSet[0].insert(ele[0]);
			_nodeSet[0].insert(ele[1]);
			_nodeSet[0].insert(ele[2]);
			_nodeSet[0].insert(ele[3]);
		}			
	}
	return true;
}

void BaseTetMeshHandle::writeTetMeshSetInfo(std::string& filename)
{
	std::ofstream fout(filename.c_str());
	fout << _elementSet.size() << std::endl;
	for (int i = 0; i < _elementSet.size(); i++)
	{
		fout << _elementSet[i].getName() << " " << _elementSet[i].getNumElements() << std::endl;
		std::set<int> elementSet;
		_elementSet[i].getElements(elementSet);
		std::set<int>::iterator it = elementSet.begin();
		for (; it != elementSet.end(); ++it)
			fout << *it << std::endl;
	}
}

void BaseTetMeshHandle::writeTetMeshSetInfoDebug(std::string & filename)
{
	std::vector<bool> eleFlag(_tetMesh->tetElementNum, false);
	for (int i = 1; i < _elementSet.size(); i++)
	{
		std::set<int> list;
		_elementSet[i].getElements(list);
		std::set<int>::iterator it = list.begin();
		for (; it != list.end(); it++)
			eleFlag[*it] = true;
	}
	_elementSet[0].clear();
	for (int i = 0; i < _tetMesh->tetElementNum; i++)
	{
		if (eleFlag[i]) continue;
		_elementSet[0].insert(i);
	}
	std::ofstream fout(filename.c_str());
	fout << _elementSet.size() << std::endl;
	for (int i = 0; i < _elementSet.size(); i++)
	{
		fout << _elementSet[i].getName() << " " << _elementSet[i].getNumElements() << std::endl;
		std::set<int> elementSet;
		_elementSet[i].getElements(elementSet);
		std::set<int>::iterator it = elementSet.begin();
		for (; it != elementSet.end(); ++it)
			fout << *it << std::endl;
	}
}

void BaseTetMeshHandle::renderTetMesh(QOpenGLFunctions * f)
{
	if (_tetMesh->_transparent)
	{
		f->glDisable(GL_COLOR_MATERIAL);
		GLfloat matAmbient[] = { 0.2f, 0.2f, 0.2f, 0.5f };
		GLfloat matDiffuse[] = { 1.0f, 1.0f, 1.0f, 0.5f };
		GLfloat matSpecular[] = { 0.0f, 0.0f, 0.0f, 0.5f };
		f->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		f->glEnable(GL_BLEND);

		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, matAmbient);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, matDiffuse);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, matSpecular);

		f->glEnable(GL_CULL_FACE);
		f->glCullFace(GL_BACK);

		f->glEnable(GL_POLYGON_OFFSET_FILL);
		f->glPolygonOffset(0.0, 1.0);
	}

	int selectId = renderElementSetId - 1;
	for (int i = 0; i < _elementSet.size(); i++)
	{
		int flag = 20 + i;
		QColor color;
		//selectId
		if (selectId == -1 || selectId == i)
		{
			color = QColor::colorNames()[flag];
			if (color == Qt::gray) color = QColor::colorNames()[flag + 13];
		}
		else
		{
			color = Qt::gray;
			continue;
		}

		std::set<int> list;
		_elementSet[i].getElements(list);
		std::set<int>::iterator it = list.begin();
		for (; it != list.end(); it++)
		{
			if (!_boundaryElementFlag[*it])
				continue;
			Vector4i indices = _tetMesh->getTetElement(*it);
			BaseTetElementParam* param = getTetElementParam(*it);
			for (int fi = 0; fi < 4; fi++)
			{
				if (_elementNeighbor[*it][fi] != -1) continue;
				int n0 = indices.data()[param->faceIndices[3 * fi + 0]];
				int n1 = indices.data()[param->faceIndices[3 * fi + 1]];
				int n2 = indices.data()[param->faceIndices[3 * fi + 2]];
				Vector3 v0 = _tetMesh->getTetPoint(n0);
				Vector3 v1 = _tetMesh->getTetPoint(n1);
				Vector3 v2 = _tetMesh->getTetPoint(n2);

				Vector3 normal = (v1 - v0).cross(v2 - v0);
				normal.normalize();

				if (selectId != -1 && selectId != i)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				}

				glBegin(GL_TRIANGLES);
				glColor4f(color.red() / 255.0, color.green() / 255.0, color.blue() / 255.0, color.alpha() / 255.0);
				glNormal3f(normal[0], normal[1], normal[2]);
				glVertex3f(v0[0], v0[1], v0[2]);
				glVertex3f(v1[0], v1[1], v1[2]);
				glVertex3f(v2[0], v2[1], v2[2]);
				glEnd();

				glBegin(GL_TRIANGLES);
				glColor4f(color.red() / 255.0, color.green() / 255.0, color.blue() / 255.0, color.alpha() / 255.0);
				glNormal3f(-normal[0], -normal[1], -normal[2]);
				glVertex3f(v0[0], v0[1], v0[2]);
				glVertex3f(v1[0], v1[1], v1[2]);
				glVertex3f(v2[0], v2[1], v2[2]);
				glEnd();


				if (selectId != -1 && selectId != i)
				{
					glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				}

				//if (!_tetMesh->_transparent)
				//{
				//	QColor lineColor = Qt::gray;
				//	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				//	glBegin(GL_POLYGON);
				//	glColor4f(lineColor.red() / 255.0, lineColor.green() / 255.0, lineColor.blue() / 255.0, lineColor.alpha() / 255.0);
				//	glNormal3f(normal[0], normal[1], normal[2]);
				//	glVertex3f(v0[0], v0[1], v0[2]);
				//	glVertex3f(v1[0], v1[1], v1[2]);
				//	glVertex3f(v2[0], v2[1], v2[2]);
				//	glEnd();
				//	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				//}
			}
		}
	}

	if (_tetMesh->_transparent)
	{
		f->glDisable(GL_POLYGON_OFFSET_FILL);
		f->glDisable(GL_BLEND);
		f->glDisable(GL_CULL_FACE);
		f->glDisable(GL_LINE_SMOOTH);
		f->glEnable(GL_COLOR_MATERIAL);
		f->glFlush();
	}
}
