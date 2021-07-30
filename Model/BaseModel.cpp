#include "BaseModel.h"

BaseModel::BaseModel():_valid(false),
	_tx(0),_ty(0),_tz(0),_scale(1.0),_rx(0),_ry(0), _rz(0), _rsita(0), enableGravity(true)
{

}

BaseModel::~BaseModel()
{

}

void BaseModel::refreshBuffer(BaseSurfaceMeshBufferPool* sPool, BaseTetMeshBufferPool* tPool, BaseMedialMeshBufferPool* mPool)
{
	((BaseSurfaceMesh*)this)->refresh(sPool);
	if(isTetMeshValid()) ((BaseTetMesh*)this)->refresh(tPool);
	if(isMedialMeshValid()) ((BaseMedialMesh*)this)->refresh(mPool);
}

void BaseModel::init()
{	
	computeBBox();
	uniform();
	initVBO();
}

void BaseModel::initMeshesHandel()
{
	_smHandle = new BaseSurfaceHandle(this);
	_smHandle->init();
	if (isTetMeshValid())
	{
		_tmHandle = new BaseTetMeshHandle(this);
		_tmHandle->init();
		computeTetMeshInterpolation();
	}
	if (isMedialMeshValid())
	{
		_mmHandle = new BaseMedialMeshHandle(this);
		_mmHandle->init();
	}
}

qeal BaseModel::uniform()
{
	qeal div = BaseSurfaceMesh::uniform();
	if(isTetMeshValid()) BaseTetMesh::uniform(div);
	if (isMedialMeshValid()) BaseMedialMesh::uniform(div);
	return div;
}

void BaseModel::conbineTetAndMat(BaseTetMeshBufferPool * tetPool, BaseMedialMeshBufferPool* matPool, const std::string newTetNodeFilename, const std::string newTetElementFilename, const std::string newMedialMeshFilename)
{
	if (!(isTetMeshValid() && isMedialMeshValid()))
		return;
	if (isBindedTetMesh())
	{
		assert(bindedTM.size() == medialPointsNum);
		for (size_t i = 0; i < medialPointsNum; i++)
		{
			if (bindedTM[i] < tetPointsNum)
			{
				Vector3 c = getMedialPoint(i);
				Vector3 p = getTetPoint(bindedTM[i]);
				if ((c - p).norm() > 1e-12)
				{
					setBindedTetMesh(false);
					break;
				}
			}
			else
			{
				setBindedTetMesh(false);
				break;
			}
		}
	}

	if (isBindedTetMesh())
	{
		bindedInverseTM.resize(tetPointsNum, -1);
		for (size_t i = 0; i < bindedTM.size(); i++)
			bindedInverseTM[bindedTM[i]] = i;
		return;
	}

	for (size_t i = 0; i < medialPointsNum; i++)
	{
		Vector3 c = getMedialPoint(i);
		qeal dist = QEAL_MAX;
		Vector3 np;
		int index = -1;
		for (int j = 0; j < tetPointsNum; j++)
		{
			Vector3 v = getTetPoint(j);
			qeal temp = (c - v).norm();
			if (temp < dist)
			{
				dist = temp;
				np = v;
				index = j;
				
			}
		}

		setMedialPoint(i, np.data());
		bindedTM[i] = index;
	}
	writeMeshToMatFormat(newMedialMeshFilename);

	bindedInverseTM.resize(tetPointsNum, -1);
	for (size_t i = 0; i < bindedTM.size(); i++)
		bindedInverseTM[bindedTM[i]] = i;
}

void BaseModel::alignAllMesh(qeal* tetPointsPtr)
{
	updateSurfaceFromTetMesh(tetPointsPtr + tetPoints.offset);
	updateMedialMeshFromTetMesh(tetPointsPtr + tetPoints.offset);
}

void BaseModel::computeTetMeshInterpolation()
{
	if (!isTetMeshValid())
		return;
	std::string filename = dir + "tetMeshInterpolation.txt";
	tetMeshInterpolationIndices.resize(pointsNum);
	tetMeshInterpolationWeight.resize(4 * pointsNum);
	// read from file
	std::ifstream fin(filename.c_str());
	if (fin.is_open())
	{
		int num;
		fin >> num;
		assert(num == pointsNum);
		for (uint32_t i = 0; i < pointsNum; i++)
		{
			uint32_t index;
			uint32_t eid;
			qeal w0, w1, w2, w3;
			uint32_t temp;
			fin >> index >> eid >> w0 >> w1 >> w2 >> w3 >> temp;
			assert(index == i);
			tetMeshInterpolationIndices[i] = eid;
			tetMeshInterpolationWeight[4 * i] = w0;
			tetMeshInterpolationWeight[4 * i + 1] = w1;
			tetMeshInterpolationWeight[4 * i + 2] = w2;
			tetMeshInterpolationWeight[4 * i + 3] = w3;
		}
		fin.close();
		return;
	}

	for (int i = 0; i < pointsNum; i++)
	{
		Vector3 p = Vector3(points.buffer[3 * i], points.buffer[3 * i + 1], points.buffer[3 * i + 2]);
		int index = -1;
		Vector4 w;
		for (int j = 0; j < tetElementNum; j++)
		{
			w.setZero();
			if (!isInsideTetElement(j, p.data(), w.data()))
				continue;
			index = j;
			break;
		}

		if (index == -1)
		{
			int tetPointId = searchCloseTetNode(p.data());
			points.buffer[3 * i] = tetPoints.buffer[3 * tetPointId];
			points.buffer[3 * i + 1] = tetPoints.buffer[3 * tetPointId + 1];
			points.buffer[3 * i + 2] = tetPoints.buffer[3 * tetPointId + 2];
			p = Vector3(points.buffer[3 * i], points.buffer[3 * i + 1], points.buffer[3 * i + 2]);
			int eid = tetPointElementList[tetPointId].buffer[0];
			w.setZero();
			computeTetElementBarycentricWeights(eid, p.data(), w.data());
			index = eid;
		}
		tetMeshInterpolationIndices[i] = index;
		tetMeshInterpolationWeight[4 * i] = w.data()[0];
		tetMeshInterpolationWeight[4 * i + 1] = w.data()[1];
		tetMeshInterpolationWeight[4 * i + 2] = w.data()[2];
		tetMeshInterpolationWeight[4 * i + 3] = w.data()[3];
	}
	// write file
	std::ofstream fout(filename.c_str());
	fout << pointsNum << std::endl;
	for (uint32_t i = 0; i < pointsNum; i++)
	{
		fout << i << " " << tetMeshInterpolationIndices[i] << " " << tetMeshInterpolationWeight[4 * i] << " " << tetMeshInterpolationWeight[4 * i + 1] << " " << tetMeshInterpolationWeight[4 * i + 2] << " " << tetMeshInterpolationWeight[4 * i + 3] << " " << -1 << std::endl;
	}
	fout.close();
}

void BaseModel::updateSurfaceFromTetMesh()
{
	if (!isTetMeshValid())
		return;
	for (int i = 0; i < pointsNum; i++)
	{
		int eid = tetMeshInterpolationIndices[i];
		for (int j = 0; j < 3; j++)
		{
			points.buffer[3 * i + j] = tetMeshInterpolationWeight[4 * i] * tetPoints.buffer[3 * tetElementIndices.buffer[4 * eid] + j];
		}
		for (int j = 0; j < 3; j++)
		{
			for(int k = 1; k < 4; k++)
				points.buffer[3 * i + j] += tetMeshInterpolationWeight[4 * i + k] * tetPoints.buffer[3 * tetElementIndices.buffer[4 * eid + k] + j];
		}
	}
}

void BaseModel::updateSurfaceFromTetMesh(qeal * currentTetPoints)
{
	if (!isTetMeshValid())
		return;
	for (int i = 0; i < pointsNum; i++)
	{
		int eid = tetMeshInterpolationIndices[i];
		for (int j = 0; j < 3; j++)
		{
			points.buffer[3 * i + j] = tetMeshInterpolationWeight[4 * i] * currentTetPoints[3 * tetElementIndices.buffer[4 * eid] + j];
		}
		for (int j = 0; j < 3; j++)
		{
			for (int k = 1; k < 4; k++)
				points.buffer[3 * i + j] += tetMeshInterpolationWeight[4 * i + k] * currentTetPoints[3 * tetElementIndices.buffer[4 * eid + k] + j];
		}
	}
}

void BaseModel::updateMedialMeshFromTetMesh(qeal * currentTetPoints)
{
	if (!isTetMeshValid() || !isMedialMeshValid() || !isBindedTetMesh())
		return;
	for (int i = 0; i < bindedTM.size(); i++)
	{
		int nid = bindedTM[i];
		setMedialPoint(i, currentTetPoints + 3 * nid);
	}

}

void BaseModel::translateModel(const qeal x, const qeal y, const qeal z)
{	
	BaseSurfaceMesh::translateMesh(x - _tx, y - _ty, z - _tz);
	if (isTetMeshValid()) BaseTetMesh::translateMesh(x - _tx, y - _ty, z - _tz);
	if (isMedialMeshValid()) BaseMedialMesh::translateMesh(x - _tx, y - _ty, z - _tz);

	_tx = x;
	_ty = y;
	_tz = z;
}

void BaseModel::scaleModel(qeal s)
{
	qeal cx, cy, cz;
	getCenter(cx, cy, cz);
	BaseSurfaceMesh::scaleMesh(s / _scale, cx, cy, cz);
	if (isTetMeshValid()) BaseTetMesh::scaleMesh(s / _scale, cx, cy, cz);
	if (isMedialMeshValid()) BaseMedialMesh::scaleMesh(s / _scale, cx, cy, cz);

	_scale = s;
}

void BaseModel::rotateModel(const qeal rx, const qeal ry, const qeal rz, const qeal sita, bool oldClear)
{
	qglviewer::Vec dir(rx, ry, rz);
	if (dir.norm() < 1e-6)
		return;
	dir.normalize();
	qglviewer::Quaternion rm(dir, sita * M_PI / 180.0);

	qglviewer::Quaternion old_rm;
	qglviewer::Vec old_dir(_rx, _ry, _rz);
	if (old_dir.norm() > 1e-6)
	{
		old_dir.normalize();
		old_rm = qglviewer::Quaternion(old_dir, _rsita* M_PI / 180.0);
	}

	old_rm = qglviewer::Quaternion();

	qeal cx, cy, cz;
	getCenter(cx, cy, cz);
	BaseSurfaceMesh::rotateMesh(old_rm, rm, cx, cy, cz);
	if (isTetMeshValid()) BaseTetMesh::rotateMesh(old_rm, rm, cx, cy, cz);
	if (isMedialMeshValid()) BaseMedialMesh::rotateMesh(old_rm, rm, cx, cy, cz);
	if (!oldClear)
	{
		_rx = rx; _ry = ry; _rz = rz; _rsita = sita;
	}	
}

