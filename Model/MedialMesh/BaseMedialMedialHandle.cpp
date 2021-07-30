#include "BaseMedialMedialHandle.h"

BaseMedialMeshHandle::BaseMedialMeshHandle(BaseMedialMesh * mesh):
	_medialMesh(mesh),
	renderMedialMeshSetId(0)
{
}

BaseMedialMesh * BaseMedialMeshHandle::getMesh() const
{
	return _medialMesh;
}

void BaseMedialMeshHandle::init()
{
	BaseMedialPointSet pointSet("defalut_sphere");
	BaseMedialPrimitiveSet primitiveSet("defalut_primitive");
	for (size_t id = 0; id < _medialMesh->medialPointsNum; id++)
		pointSet.insert(id);
	for (size_t id = 0; id < _medialMesh->medialPrimitivesNum; id++)
		primitiveSet.insert(id);
	_pointSet.push_back(pointSet);
	_primitiveSet.push_back(primitiveSet);
}

bool BaseMedialMeshHandle::readMedialMeshSetInfo(std::string & filename)
{
	std::ifstream fin(filename.c_str());
	if (!fin.is_open())
		return false;

	std::vector<int> primitiveFlag(_medialMesh->medialPrimitivesNum, -1);
	int setNum = 0;
	fin >> setNum;
	if (setNum == 0) return true;

	_primitiveSet.clear();
	_pointSet.clear();

	for (int i = 0; i < setNum; i++)
	{
		std::string name;
		int spNum, priNum;
		fin >> name >> spNum >> priNum;
		std::string sp_name = name + "_spheres";
		std::string pri_name = name + "_primitives";
		BaseMedialPointSet pointSet(sp_name);
		BaseMedialPrimitiveSet primitiveSet(pri_name);

		for (int j = 0; j < spNum; j++)
		{
			int spid;
			char ch;
			fin >> ch >> spid;
			pointSet.insert(spid);
		}

		for (int j = 0; j < priNum; j++)
		{
			int priid;
			char ch;
			fin >> ch >> priid;
			primitiveSet.insert(priid);
			primitiveFlag[priid] = i;
			Vector3i priIndex = _medialMesh->getMedialPrimitive(priid);
			pointSet.insert(priIndex[0]);
			pointSet.insert(priIndex[1]);
			if(priIndex[2] != -1)
				pointSet.insert(priIndex[2]);
		}
		_pointSet.push_back(pointSet);
		_primitiveSet.push_back(primitiveSet);
	}	

	for (int i = 0; i < _medialMesh->medialPrimitivesNum; i++)
	{
		if (primitiveFlag[i] != -1) continue;
		Vector3i priIndex = _medialMesh->getMedialPrimitive(i);
		_primitiveSet[0].insert(i);
		_pointSet[0].insert(priIndex[0]);
		_pointSet[0].insert(priIndex[1]);
		if(priIndex[2] != -1)
			_pointSet[0].insert(priIndex[2]);
	}

	return true;
}

void BaseMedialMeshHandle::writeMedialMeshSetInfo(std::string & filename)
{
	std::ofstream fout(filename.c_str());
	fout << _pointSet.size() << std::endl;
	for (int i = 0; i < _pointSet.size(); i++)
	{
		size_t sublen = _pointSet[i].getName().find_first_of("_");
		std::string name = _pointSet[i].getName().substr(0, sublen + 1);
		fout << name << " " << _pointSet[i].getNumElements() << " " << _primitiveSet[i].getNumElements() << std::endl;

		std::set<int> pointSet;
		_pointSet[i].getElements(pointSet);
		std::set<int>::iterator it = pointSet.begin();
		for (; it != pointSet.end(); ++it)
			fout << "v " << *it << std::endl;
		std::set<int> primitiveSet;
		_primitiveSet[i].getElements(primitiveSet);
		it = primitiveSet.begin();
		for (; it != primitiveSet.end(); ++it)
			fout << "s " << *it << std::endl;
	}
}

void BaseMedialMeshHandle::renderMedialMesh(QOpenGLFunctions * f)
{
	if (_medialMesh->_transparent)
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

	int selectSetId = renderMedialMeshSetId - 2;

	if (selectSetId == -1)
	{
		for (size_t nid = 0; nid < _medialMesh->medialPointsNum; nid++)
		{			
			SphereElement sphere(_medialMesh->getMedialPoint(nid), _medialMesh->getMedialRadius(nid));
			sphere.show();
		}
	}
	else if(selectSetId == -2)
	{
		for (size_t nid = 0; nid < _medialMesh->medialPointsNum; nid++)
		{
			SphereElement sphere(_medialMesh->getMedialPoint(nid), _medialMesh->getMedialRadius(nid));
			sphere.show();
		}

		for (size_t cid = 0; cid < _medialMesh->medialConesNum; cid++)
		{
			Vector2i indices = _medialMesh->getMedialCone(cid);
			Vector3 cl = _medialMesh->getMedialPoint(indices.data()[0]);
			qeal rl = _medialMesh->getMedialRadius(indices.data()[0]);
			Vector3 cr = _medialMesh->getMedialPoint(indices.data()[1]);
			qeal rr = _medialMesh->getMedialRadius(indices.data()[1]);

			ConeElement cone(cl, rl, cr, rr);
			cone.show();
		}

		for (size_t sid = 0; sid < _medialMesh->medialSlabsNum; sid++)
		{
			Vector3i indices = _medialMesh->getMedialSlab(sid);
			Vector3 c0 = _medialMesh->getMedialPoint(indices.data()[0]);
			qeal r0 = _medialMesh->getMedialRadius(indices.data()[0]);
			Vector3 c1 = _medialMesh->getMedialPoint(indices.data()[1]);
			qeal r1 = _medialMesh->getMedialRadius(indices.data()[1]);
			Vector3 c2 = _medialMesh->getMedialPoint(indices.data()[2]);
			qeal r2 = _medialMesh->getMedialRadius(indices.data()[2]);
			ConeElement cone01(c0, r0, c1, r1);
			cone01.show();
			ConeElement cone02(c0, r0, c2, r2);
			cone02.show();
			ConeElement cone12(c1, r1, c2, r2);
			cone12.show();

			TwoSplintElements st(c0, r0, c1, r1, c2, r2);
			st.show();
		}
	}
	else
	{
		int setId = selectSetId / 2;
		if (selectSetId % 2 == 1)
		{
			std::set<int> set;
			_pointSet[setId].getElements(set);
			std::set<int>::iterator it = set.begin();
			for (; it != set.end(); ++it)
			{
				SphereElement sphere(_medialMesh->getMedialPoint(*it), _medialMesh->getMedialRadius(*it));
				sphere.show();
			}
		}
		else if(selectSetId % 2 == 0) 
		{
			std::set<int> set;
			_primitiveSet[setId].getElements(set);
			std::set<int>::iterator it = set.begin();
			std::set<int> primitiveSphere;
			for (; it != set.end(); ++it)
			{
				if (*it < _medialMesh->medialConesNum)
				{
					// cone
					Vector2i indices = _medialMesh->getMedialCone(*it);
					Vector3 cl = _medialMesh->getMedialPoint(indices.data()[0]);
					qeal rl = _medialMesh->getMedialRadius(indices.data()[0]);
					Vector3 cr = _medialMesh->getMedialPoint(indices.data()[1]);
					qeal rr = _medialMesh->getMedialRadius(indices.data()[1]);

					ConeElement cone(cl, rl, cr, rr);
					cone.show();

					primitiveSphere.insert(indices[0]);
					primitiveSphere.insert(indices[1]);
				}
				else
				{
					// slab
					int slabId = *it - _medialMesh->medialConesNum;
					Vector3i indices = _medialMesh->getMedialSlab(slabId);
					Vector3 c0 = _medialMesh->getMedialPoint(indices.data()[0]);
					qeal r0 = _medialMesh->getMedialRadius(indices.data()[0]);
					Vector3 c1 = _medialMesh->getMedialPoint(indices.data()[1]);
					qeal r1 = _medialMesh->getMedialRadius(indices.data()[1]);
					Vector3 c2 = _medialMesh->getMedialPoint(indices.data()[2]);
					qeal r2 = _medialMesh->getMedialRadius(indices.data()[2]);
					ConeElement cone01(c0, r0, c1, r1);
					cone01.show();
					ConeElement cone02(c0, r0, c2, r2);
					cone02.show();
					ConeElement cone12(c1, r1, c2, r2);
					cone12.show();

					TwoSplintElements st(c0, r0, c1, r1, c2, r2);
					st.show();

					primitiveSphere.insert(indices[0]);
					primitiveSphere.insert(indices[1]);
					primitiveSphere.insert(indices[2]);
				}
			}

			it = primitiveSphere.begin();
			for (; it != primitiveSphere.end(); ++it)
			{
				SphereElement sphere(_medialMesh->getMedialPoint(*it), _medialMesh->getMedialRadius(*it));
				sphere.show();
			}
		}
	}


	if (_medialMesh->_transparent)
	{
		f->glDisable(GL_POLYGON_OFFSET_FILL);
		f->glDisable(GL_BLEND);
		f->glDisable(GL_CULL_FACE);
		f->glDisable(GL_LINE_SMOOTH);
		f->glEnable(GL_COLOR_MATERIAL);
		f->glFlush();
	}
}

