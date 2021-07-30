#include "BaseTetMesh.h"
#include "Commom\FileIO.h"
#include <assert.h>
#include <algorithm>

bool BaseTetMesh::readMeshFromTetFormat(const std::string filename, BaseTetMeshBufferPool * pool)
{
	std::string dir, name, format;
	getFilenameInfo(filename, dir, name, format);
	if (format == std::string("node"))
	{
		_hasNodes = readNodeFromTetFormat(filename, pool);
		_valid = _hasNodes & _hasElements;
		return _hasNodes;
	}
	else if (format == std::string("ele"))
	{
		_hasElements = readElementFromTetFormat(filename, pool);
		_valid = _hasNodes & _hasElements;
		return _hasElements;
	}		
	return false;
}

bool BaseTetMesh::readNodeFromTetFormat(const std::string filename, BaseTetMeshBufferPool * pool)
{
	std::ifstream fin(filename.c_str());
	std::string line;
	std::stringstream sstream;
	// read nodes
	if (!fin.is_open())
		return false;
	int pointNum, dim, temp1, temp2;
	std::vector<qeal> points;
	std::getline(fin, line);
	sstream << line;
	sstream >> pointNum >> dim >> temp1 >> temp2;
	assert(dim == 3);
	points.resize(pointNum * dim);
	tetPointsNum = pointNum;
	for (int i = 0; i < pointNum; i++)
	{
		sstream.clear();
		std::getline(fin, line);
		int index;
		qeal x, y, z;
		sstream << line;
		sstream >> index >> x >> y >> z;
		points[3 * index + 0] = x;
		points[3 * index + 1] = y;
		points[3 * index + 2] = z;
	}
	sstream.clear();
	fin.close();
	registerToOverallBuffer<qeal>(points, pool->tetPointsBuffer, this->tetPoints);
	pool->totalTetPointsNum += pointNum;
	tetPointsIdOffset = tetPoints.offset / 3;
	return true;
}

bool BaseTetMesh::readElementFromTetFormat(const std::string filename, BaseTetMeshBufferPool * pool)
{
	std::ifstream fin(filename.c_str());
	std::string line;
	std::stringstream sstream;
	// read element
	if (!fin.is_open())
		return false;

	std::getline(fin, line);

	int tetNum, dim, temp1;
	sstream << line;
	sstream >> tetNum >> dim >> temp1;
	assert(dim == 4);
	tetElementNum = tetNum;
	std::vector<int> tetIndices(tetNum * dim);
	std::vector<std::vector<int>> pointTetList;
	std::vector<std::vector<int>> pointNeighborList;
	for (uint32_t i = 0; i < tetNum; i++)
	{
		sstream.clear();
		std::getline(fin, line);
		uint32_t index, v0, v1, v2, v3;
		sstream << line;
		sstream >> index >> v0 >> v1 >> v2 >> v3;

		tetIndices[4 * index + 0] = v0;
		tetIndices[4 * index + 1] = v1;
		tetIndices[4 * index + 2] = v2;
		tetIndices[4 * index + 3] = v3;

		if (v0 >= pointTetList.size())
		{
			pointTetList.resize(v0 + 1);
			pointNeighborList.resize(v0 + 1);
		}
			
		if (v1 >= pointTetList.size())
		{
			pointTetList.resize(v1 + 1);
			pointNeighborList.resize(v1 + 1);
		}
			
		if (v2 >= pointTetList.size())
		{
			pointTetList.resize(v2 + 1);
			pointNeighborList.resize(v2 + 1);
		}
			
		if (v3 >= pointTetList.size())
		{
			pointTetList.resize(v3 + 1);
			pointNeighborList.resize(v3 + 1);
		}
			
		pointTetList[v0].push_back(i);
		pointTetList[v1].push_back(i);
		pointTetList[v2].push_back(i);
		pointTetList[v3].push_back(i);

		pointNeighborList[v0].push_back(v1);
		pointNeighborList[v0].push_back(v2);
		pointNeighborList[v0].push_back(v3);

		pointNeighborList[v1].push_back(v0);
		pointNeighborList[v1].push_back(v2);
		pointNeighborList[v1].push_back(v3);

		pointNeighborList[v2].push_back(v0);
		pointNeighborList[v2].push_back(v1);
		pointNeighborList[v2].push_back(v3);

		pointNeighborList[v3].push_back(v0);
		pointNeighborList[v3].push_back(v1);
		pointNeighborList[v3].push_back(v2);
	}

	for (size_t i = 0; i < pointNeighborList.size(); i++)
	{
		std::sort(pointNeighborList[i].begin(), pointNeighborList[i].end());
		pointNeighborList[i].erase(std::unique(pointNeighborList[i].begin(), pointNeighborList[i].end()), pointNeighborList[i].end());
	}

	sstream.clear();
	fin.close();
	sstream.clear();
	registerToOverallBuffer<int>(tetIndices, pool->tetElementIndicesBuffer, this->tetElementIndices);
	registerToOverallBuffer<int>(pointTetList, pool->tetPointElementListBuffer, this->tetPointElementList);
	registerToOverallBuffer<int>(pointNeighborList, pool->tetPointNeighborListBuffer, this->tetPointNeighborList);
	pool->totalTetElementNum += tetNum;
	tetElementsIdOffset = tetElementIndices.offset / 4;
	return true;
}

bool BaseTetMesh::writeMeshToTetFormat(const std::string filename)
{
	std::string dir, name, format;
	getFilenameInfo(filename, dir, name, format);
	if (format == std::string("node"))
		return writeNodeToTetFormat(filename);
	else if (format == std::string("ele"))
		return writeElementToTetFormat(filename);
	return false;
}

bool BaseTetMesh::writeNodeToTetFormat(const std::string filename)
{
	std::ofstream fout(filename.c_str());
	if (!fout.is_open())
		return false;
	fout << tetPointsNum << " 3 0 0" << std::endl;
	for (int i = 0; i < tetPointsNum; i++)
		fout << i << " " << tetPoints.buffer[3 * i] << " " << tetPoints.buffer[3 * i + 1] << " " << tetPoints.buffer[3 * i + 2] << std::endl;
	fout.close();
	return true;
}

bool BaseTetMesh::writeElementToTetFormat(const std::string filename)
{
	std::ofstream fout(filename.c_str());
	if (!fout.is_open())
		return false;
	fout << tetElementNum << " 4 0" << std::endl;
	for (int i = 0; i < tetElementNum; i++)
		fout << i << " " << tetElementIndices.buffer[4 * i] << " " << tetElementIndices.buffer[4 * i + 1] << " " << tetElementIndices.buffer[4 * i + 2] << " " << tetElementIndices.buffer[4 * i + 3] <<std::endl;
	fout.close();
	return true;
}

bool BaseTetMesh::isInsideTetElement(const int eid, const qeal * point, qeal * weight)
{
	computeTetElementBarycentricWeights(eid, point, weight);

	if (IS_QEAL_ZERO(weight[0]))
		weight[0] = 0.0;
	if (IS_QEAL_ZERO(weight[1]))
		weight[1] = 0.0;
	if (IS_QEAL_ZERO(weight[2]))
		weight[2] = 0.0;
	if (IS_QEAL_ZERO(weight[3]))
		weight[3] = 0.0;

	bool result = false;
	if ((weight[0] >= 0.0) && (weight[1] >= 0.0) && (weight[2] >= 0.0) && (weight[3] >= 0.0))
		result = true;

	return result;
}

void BaseTetMesh::computeTetElementBarycentricWeights(const int eid, const qeal * point, qeal * weight)
{
	Vector4 temp;
	temp.data()[0] = point[0];
	temp.data()[1] = point[1];
	temp.data()[2] = point[2];
	temp.data()[3] = 1;
	Matrix4 referenceShape;
	Matrix4 invShapeFunction;
	Vector4i indices = getTetElement(eid);
	for (int i = 0; i < 4; i++)
	{
		Vector3 p = getTetPoint(indices[i]);
		for (int j = 0; j < 3; j++)
			referenceShape.data()[i * 4 + j] = p.data()[j];
		referenceShape.data()[i * 4 + 3] = 1;
	}
	invShapeFunction = referenceShape.inverse();
	Vector4 w = invShapeFunction * temp;
	for (int i = 0; i < 4; i++)
		weight[i] = w.data()[i];
}

int BaseTetMesh::searchCloseTetNode(const qeal * point)
{
	qeal close = QEAL_MAX;
	int id = -1;
	Vector3 p = Vector3(point[0], point[1], point[2]);
	for (int i = 0; i < tetPointsNum; i++)
	{
		Vector3 tp = Vector3(tetPoints.buffer[3 * i], tetPoints.buffer[3 * i + 1], tetPoints.buffer[3 * i + 2]);
		qeal len = (p - tp).norm();
		if (len < close)
		{
			close = len;
			id = i;
		}
	}
	return id;
}

void BaseTetMesh::uniform(const qeal div)
{
	for (int i = 0; i < tetPointsNum; i++)
	{
		tetPoints.buffer[3 * i] /= div;
		tetPoints.buffer[3 * i + 1] /= div;
		tetPoints.buffer[3 * i + 2] /= div;
	}
}

void BaseTetMesh::translateMesh(const qeal x, const qeal y, const qeal z)
{
	for (int i = 0; i < tetPointsNum; i++)
	{
		tetPoints.buffer[3 * i] += x;
		tetPoints.buffer[3 * i + 1] += y;
		tetPoints.buffer[3 * i + 2] += z;
	}
}

void BaseTetMesh::scaleMesh(const qeal s, const qeal cx, const qeal cy, const qeal cz)
{
	for (int i = 0; i < tetPointsNum; i++)
	{
		tetPoints.buffer[3 * i] = (tetPoints.buffer[3 * i] - cx) * s + cx;
		tetPoints.buffer[3 * i + 1] = (tetPoints.buffer[3 * i + 1] - cy) * s + cy;
		tetPoints.buffer[3 * i + 2] = (tetPoints.buffer[3 * i + 2] - cz) * s + cz;
	}
}

void BaseTetMesh::rotateMesh(const qglviewer::Quaternion oldR, const qglviewer::Quaternion R, const qeal cx, const qeal cy, const qeal cz)
{
	for (int i = 0; i < tetPointsNum; i++)
	{
		qglviewer::Vec pos(tetPoints.buffer[3 * i] - cx, tetPoints.buffer[3 * i + 1] - cy, tetPoints.buffer[3 * i + 2] - cz);
		pos = oldR.inverseRotate(pos);
		pos = R.rotate(pos);

		tetPoints.buffer[3 * i] = pos.x + cx;
		tetPoints.buffer[3 * i + 1] = pos.y + cy;
		tetPoints.buffer[3 * i + 2] = pos.z + cz;
	}
}

void BaseTetMeshBuffer::copyTetPointsToBuffer(qeal * buffer, int size)
{
	if (size == 0)
		std::copy(tetPoints.buffer, tetPoints.buffer + tetPoints.span, buffer);
	else std::copy(tetPoints.buffer, tetPoints.buffer + size, buffer);
}

void BaseTetMeshBuffer::copyBufferToTetPoints(qeal * buffer, int size)
{
	if (size == 0)
		std::copy(buffer, buffer + tetPoints.span, tetPoints.buffer);
	else std::copy(buffer, buffer + size, tetPoints.buffer);
}

Vector3 BaseTetMeshBuffer::getTetPoint(int nid)
{
	return Vector3(tetPoints.buffer[3 * nid], tetPoints.buffer[3 * nid + 1], tetPoints.buffer[3 * nid + 2]);
}

void BaseTetMeshBuffer::setTetPoint(int nid, qeal* p)
{
	for(size_t i = 0; i < 3; i++)
		tetPoints.buffer[3 * nid + i] = p[i];
}

Vector4i BaseTetMeshBuffer::getTetElement(int eleid)
{
	return Vector4i(tetElementIndices.buffer[4*eleid], tetElementIndices.buffer[4 * eleid + 1], tetElementIndices.buffer[4 * eleid + 2], tetElementIndices.buffer[4 * eleid + 3]);
}

void BaseTetMeshBuffer::setTetElement(int eleid, int* indices)
{
	for (size_t i = 0; i < 4; i++)
		tetElementIndices.buffer[4 * eleid + i] = indices[i];
}

VectorX BaseTetMeshBuffer::getTetElementPoint(int eleid)
{
	VectorX points(12);
	Vector4i indices = getTetElement(eleid);
	for (size_t i = 0; i < 4; i++)
		for (size_t j = 0; j < 3; j++)
			points.data()[3 * i + j] = tetPoints.buffer[3 * indices.data()[i] + j];
	return points;
}

Vector4i BaseTetMeshBuffer::getTetElementOverall(int eleid)
{
	return Vector4i(tetElementIndices.buffer[4 * eleid] + tetPointsIdOffset, tetElementIndices.buffer[4 * eleid + 1] + tetPointsIdOffset, tetElementIndices.buffer[4 * eleid + 2] + tetPointsIdOffset, tetElementIndices.buffer[4 * eleid + 3] + tetPointsIdOffset);
}

int BaseTetMeshBuffer::getTetElementOverallId(int eid)
{
	return eid + tetElementsIdOffset;
}

int BaseTetMeshBuffer::getTetPointOverallId(int nid)
{
	return nid + tetPointsIdOffset;
}
