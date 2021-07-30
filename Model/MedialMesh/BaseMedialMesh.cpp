#include "BaseMedialMesh.h"
#include "Commom\FileIO.h"
#include "Commom\GeometryComputation.h"

void BaseMedialMeshBuffer::copyMedialPointsToBuffer(qeal* buffer, int size)
{
	if(size == 0)
		std::copy(medialPoints.buffer, medialPoints.buffer + medialPoints.span, buffer);
	else std::copy(medialPoints.buffer, medialPoints.buffer + size, buffer);
}

void BaseMedialMeshBuffer::copyBufferToMedialPoints(qeal* buffer, int size)
{
	if (size == 0)
		std::copy(buffer, buffer + medialPoints.span, medialPoints.buffer);
	else std::copy(buffer, buffer + size, medialPoints.buffer);
}

void BaseMedialMeshBuffer::copyMedialRadiusToBuffer(qeal * buffer, int size)
{
	if (size == 0)
		std::copy(medialRadius.buffer, medialRadius.buffer + medialRadius.span, buffer);
	else std::copy(medialRadius.buffer, medialRadius.buffer + size, buffer);
}

void BaseMedialMeshBuffer::copyBufferToMedialRadius(qeal * buffer, int size)
{
	if (size == 0)
		std::copy(buffer, buffer + medialRadius.span, medialRadius.buffer);
	else std::copy(buffer, buffer + size, medialRadius.buffer);
}

Vector3 BaseMedialMeshBuffer::getMedialPoint(int nid)
{
	return Vector3(medialPoints.buffer[3 * nid], medialPoints.buffer[3 * nid + 1], medialPoints.buffer[3 * nid + 2]);
}

void BaseMedialMeshBuffer::setMedialPoint(int nid, qeal* p)
{
	medialPoints.buffer[3 * nid] = p[0];
	medialPoints.buffer[3 * nid + 1] = p[1];
	medialPoints.buffer[3 * nid + 2] = p[2];
}

qeal BaseMedialMeshBuffer::getMedialRadius(int nid)
{
	return medialRadius.buffer[nid];
}

void BaseMedialMeshBuffer::setMedialRadius(int nid, qeal& r)
{
	medialRadius.buffer[nid] = r;
}

Vector3i BaseMedialMeshBuffer::getMedialPrimitive(int pid)
{
	return Vector3i(medialPrimitives.buffer[3 * pid], medialPrimitives.buffer[3 * pid + 1], medialPrimitives.buffer[3 * pid + 2]);
}

Vector2i BaseMedialMeshBuffer::getMedialCone(int cid)
{
	int pid = cid;
	return Vector2i(medialPrimitives.buffer[3 * pid], medialPrimitives.buffer[3 * pid + 1]);
}

Vector3i BaseMedialMeshBuffer::getMedialSlab(int sid)
{
	int pid = sid + medialConesNum;
	return Vector3i(medialPrimitives.buffer[3 * pid], medialPrimitives.buffer[3 * pid + 1], medialPrimitives.buffer[3 * pid + 2]);
}

int BaseMedialMeshBuffer::getMedialPointOverallId(int id)
{
	return id + medialPointIdOffset;
}

int BaseMedialMeshBuffer::getMedialPrimitiveOverallId(int pid)
{
	return pid + medialPrimitiveIdOffset;
}

bool BaseMedialMesh::readMeshFromMatFormat(const std::string filename, BaseMedialMeshBufferPool * pool)
{
	std::string dir, name, format;
	getFilenameInfo(filename, dir, name, format);
	if (format == std::string("mat"))
	{
		std::ifstream fin(filename.c_str());
		if (!fin.is_open())
			return false;
		std::string line;
		std::stringstream sstream;
		int pointsNum, conesNum, slabsNum;
		std::getline(fin, line);
		sstream << line;
		sstream >> pointsNum >> conesNum >> slabsNum;
		std::vector<qeal> points(3 * pointsNum);
		std::vector<qeal> radius(pointsNum);
		int primitivesNum = conesNum + slabsNum;

		medialPointsNum = pointsNum;
		medialConesNum = conesNum;
		medialSlabsNum = slabsNum;
		medialPrimitivesNum = primitivesNum;

		std::vector<int> primitives(3 * (conesNum + slabsNum));
		std::vector<std::vector<int>> medialPointsNeighborList(pointsNum);
		std::vector<std::vector<int>> medialPointsCones(pointsNum);
		std::vector<std::vector<int>> medialPointsSlabs(pointsNum);
		bindedTM.resize(pointsNum, -1);
		setBindedTetMesh(true);
		for (size_t i = 0; i < pointsNum; i++)
		{
			sstream.clear();
			std::getline(fin, line);
			char ch;
			qeal x, y, z, r;
			int b;
			sstream << line;
			sstream >> ch >> x >> y >> z >> r >> b;
			assert(ch == 'v');
			points.data()[3 * i + 0] = x;
			points.data()[3 * i + 1] = y;
			points.data()[3 * i + 2] = z;
			radius.data()[i] = r;
			bindedTM[i] = b;
			if(b < 0) setBindedTetMesh(false);
		}

		for (size_t i = 0; i < conesNum; i++)
		{
			sstream.clear();
			std::getline(fin, line);
			char ch;
			int x, y;
			sstream << line;
			sstream >> ch >> x >> y;
			assert(ch == 'c');

			primitives[3 * i + 0] = x;
			primitives[3 * i + 1] = y;
			if (radius[x] < radius[y])
			{
				primitives[3 * i + 0] = y;
				primitives[3 * i + 1] = x;
			}
			primitives[3 * i + 2] = -1;

			medialPointsCones[x].push_back(i);
			medialPointsCones[y].push_back(i);
			medialPointsNeighborList[x].push_back(y);
			medialPointsNeighborList[y].push_back(x);
		}

		for (size_t i = 0; i < slabsNum; i++)
		{
			sstream.clear();
			std::getline(fin, line);
			char ch;
			int x, y, z;
			sstream << line;
			sstream >> ch >> x >> y >> z;
			assert(ch == 's');

			int max = x, mid = y, min = z;
			qeal maxr = radius[x];
			qeal midr = radius[y];
			qeal minr = radius[z];

			if (maxr < midr)
			{
				qeal tempr = maxr;
				maxr = midr;
				midr = tempr;

				int temp = max;
				max = mid;
				mid = temp;
			}

			if (maxr < minr)
			{
				qeal tempr = maxr;
				maxr = minr;
				minr = tempr;

				int temp = max;
				max = min;
				min = temp;
			}

			if (midr < minr)
			{
				qeal tempr = midr;
				midr = minr;
				minr = tempr;

				int temp = mid;
				mid = min;
				min = temp;
			}

			primitives[3 * (i + conesNum) + 0] = max;
			primitives[3 * (i + conesNum) + 1] = mid;
			primitives[3 * (i + conesNum) + 2] = min;
			medialPointsSlabs[x].push_back(i);
			medialPointsSlabs[y].push_back(i);
			medialPointsSlabs[z].push_back(i);
			medialPointsNeighborList[x].push_back(y);
			medialPointsNeighborList[x].push_back(z);
			medialPointsNeighborList[y].push_back(x);
			medialPointsNeighborList[y].push_back(z);
			medialPointsNeighborList[z].push_back(x);
			medialPointsNeighborList[z].push_back(y);
		}
		sstream.clear();
		fin.close();
		for (size_t i = 0; i < pointsNum; i++)
		{
			std::sort(medialPointsNeighborList[i].begin(), medialPointsNeighborList[i].end());
			medialPointsNeighborList[i].erase(std::unique(medialPointsNeighborList[i].begin(), medialPointsNeighborList[i].end()), medialPointsNeighborList[i].end());
		}

		registerToOverallBuffer<qeal>(points, pool->medialPointsBuffer, this->medialPoints);
		registerToOverallBuffer<qeal>(radius, pool->medialRadiusBuffer, this->medialRadius);
		pool->totalMedialPoinsNum += pointsNum;
		registerToOverallBuffer<int>(primitives, pool->medialPrimitiveIndicesBuffer, this->medialPrimitives);
		pool->totalMedialConesNum += conesNum;
		pool->totalMedialSlabsNum += slabsNum;
		pool->totalMedialPrimitivesNum += primitivesNum;

		registerToOverallBuffer<int>(medialPointsNeighborList, pool->medialPointsNeighborListBuffer, this->medialPointsNeighborList);
		registerToOverallBuffer<int>(medialPointsCones, pool->medialPointsConesBuffer, this->medialPointsCones);
		registerToOverallBuffer<int>(medialPointsCones, pool->medialPointsSlabsBuffer, this->medialPointsSlabs);
		_valid = true;
	}
    else if (format == std::string("ma"))
	{
		std::ifstream fin(filename.c_str());
		if (!fin.is_open())
			return false;
		std::string line;
		std::stringstream sstream;
		int pointsNum, edgesNum, facesNum;
		std::getline(fin, line);
		sstream << line;
		sstream >> pointsNum >> edgesNum >> facesNum;
		std::vector<qeal> points(3 * pointsNum);
		std::vector<qeal> radius(pointsNum);

		medialPointsNum = pointsNum;
		for (size_t i = 0; i < pointsNum; i++)
		{
			sstream.clear();
			std::getline(fin, line);
			char ch;
			qeal x, y, z, r;
			sstream << line;
			sstream >> ch >> x >> y >> z >> r;
			assert(ch == 'v');
			points.data()[3 * i + 0] = x;
			points.data()[3 * i + 1] = y;
			points.data()[3 * i + 2] = z;
			radius.data()[i] = r;
		}
		bindedTM.resize(pointsNum, -1);
		setBindedTetMesh(false);
		std::vector<std::pair<bool, Vector2i>> edgeIndices(edgesNum);
		for (size_t i = 0; i < edgesNum; i++)
		{
			sstream.clear();
			std::getline(fin, line);
			char ch;
			int x, y;
			sstream << line;
			sstream >> ch >> x >> y;
			assert(ch == 'e');
			if (radius[x] >= radius[y])
				edgeIndices[i] = std::pair<bool, Vector2i>(true, Vector2i(x, y));
			else edgeIndices[i] = std::pair<bool, Vector2i>(true, Vector2i(y, x));
		}

		std::vector<int> faceIndices(3 * facesNum);
		for (size_t i = 0; i < facesNum; i++)
		{
			sstream.clear();
			std::getline(fin, line);
			char ch;
			int x, y, z;
			sstream << line;
			sstream >> ch >> x >> y >> z;
			assert(ch == 'f');
			faceIndices[3 * i + 0] = x;
			faceIndices[3 * i + 1] = y;
			faceIndices[3 * i + 2] = z;
			std::vector<std::pair<bool, Vector2i>>::iterator it;
			it = std::find(edgeIndices.begin(), edgeIndices.end(), std::pair<bool, Vector2i>(true, Vector2i(x, y)));
			if (it != edgeIndices.end()) it->first = false;
			it = std::find(edgeIndices.begin(), edgeIndices.end(), std::pair<bool, Vector2i>(true, Vector2i(y, x)));
			if (it != edgeIndices.end()) it->first = false;

			it = std::find(edgeIndices.begin(), edgeIndices.end(), std::pair<bool, Vector2i>(true, Vector2i(x, z)));
			if (it != edgeIndices.end()) it->first = false;
			it = std::find(edgeIndices.begin(), edgeIndices.end(), std::pair<bool, Vector2i>(true, Vector2i(z, x)));
			if (it != edgeIndices.end()) it->first = false;

			it = std::find(edgeIndices.begin(), edgeIndices.end(), std::pair<bool, Vector2i>(true, Vector2i(y, z)));
			if (it != edgeIndices.end()) it->first = false;
			it = std::find(edgeIndices.begin(), edgeIndices.end(), std::pair<bool, Vector2i>(true, Vector2i(z, y)));
			if (it != edgeIndices.end()) it->first = false;
		}
		
		medialConesNum = 0;
		std::vector<int> primitives;
		std::vector<std::vector<int>> medialPointsNeighborList(pointsNum);
		std::vector<std::vector<int>> medialPointsCones(pointsNum);
		std::vector<std::vector<int>> medialPointsSlabs(pointsNum);
		for (size_t i = 0; i < edgesNum; i++)
		{
			if (edgeIndices[i].first)
			{
				int x = edgeIndices[i].second.data()[0];
				int y = edgeIndices[i].second.data()[1];
				primitives.push_back(x);
				primitives.push_back(y);
				primitives.push_back(-1);
				medialPointsCones[x].push_back(medialConesNum);
				medialPointsCones[y].push_back(medialConesNum);
				medialConesNum++;
				medialPointsNeighborList[x].push_back(y);
				medialPointsNeighborList[y].push_back(x);
			}
		}
		medialSlabsNum = facesNum;
		for (size_t i = 0; i < facesNum; i++)
		{
			int x = faceIndices[3 * i + 0];
			int y = faceIndices[3 * i + 1];
			int z = faceIndices[3 * i + 2];

			int max = x, mid = y, min = z;
			qeal maxr = radius[x];
			qeal midr = radius[y];
			qeal minr = radius[z];

			if (maxr < midr)
			{
				qeal tempr = maxr;
				maxr = midr;
				midr = tempr;

				int temp = max;
				max = mid;
				mid = temp;
			}

			if (maxr < minr)
			{
				qeal tempr = maxr;
				maxr = minr;
				minr = tempr;

				int temp = max;
				max = min;
				min = temp;
			}

			if (midr < minr)
			{
				qeal tempr = midr;
				midr = minr;
				minr = tempr;

				int temp = mid;
				mid = min;
				min = temp;
			}

			primitives.push_back(max);
			primitives.push_back(mid);
			primitives.push_back(min);
			medialPointsSlabs[x].push_back(i);
			medialPointsSlabs[y].push_back(i);
			medialPointsSlabs[x].push_back(i);
			medialPointsNeighborList[x].push_back(y);
			medialPointsNeighborList[x].push_back(z);
			medialPointsNeighborList[y].push_back(x);
			medialPointsNeighborList[y].push_back(z);
			medialPointsNeighborList[z].push_back(x);
			medialPointsNeighborList[z].push_back(y);
		}
		medialPrimitivesNum = medialConesNum + medialSlabsNum;

		sstream.clear();
		fin.close();

		for (size_t i = 0; i < pointsNum; i++)
		{
			std::sort(medialPointsNeighborList[i].begin(), medialPointsNeighborList[i].end());
			medialPointsNeighborList[i].erase(std::unique(medialPointsNeighborList[i].begin(), medialPointsNeighborList[i].end()), medialPointsNeighborList[i].end());
		}

		registerToOverallBuffer<qeal>(points, pool->medialPointsBuffer, this->medialPoints);
		registerToOverallBuffer<qeal>(radius, pool->medialRadiusBuffer, this->medialRadius);
		pool->totalMedialPoinsNum += pointsNum;
		registerToOverallBuffer<int>(primitives, pool->medialPrimitiveIndicesBuffer, this->medialPrimitives);
		pool->totalMedialConesNum += medialConesNum;
		pool->totalMedialSlabsNum += medialSlabsNum;
		pool->totalMedialPrimitivesNum += medialPrimitivesNum;

		registerToOverallBuffer<int>(medialPointsNeighborList, pool->medialPointsNeighborListBuffer, this->medialPointsNeighborList);
		registerToOverallBuffer<int>(medialPointsCones, pool->medialPointsConesBuffer, this->medialPointsCones);
		registerToOverallBuffer<int>(medialPointsSlabs, pool->medialPointsSlabsBuffer, this->medialPointsSlabs);
		
		_valid = writeMeshToMatFormat(dir + name + ".mat");
	}

	medialPointIdOffset = medialPoints.offset / 3;
	medialPrimitiveIdOffset = medialPrimitives.offset / 3;

	faceList.resize(medialSlabsNum);
	for (int i = 0; i < medialSlabsNum; i++)
	{
		Vector3i face = getMedialSlab(i);
		faceList[i] = face;
	}
	std::set<Vector2i> edgeSet;
	for (int i = 0; i < medialConesNum; i++)
	{
		Vector2i edge = getMedialCone(i);
		bool has = false;
		for (int j = 0; j < edgeList.size(); j++)
		{
			if ((edgeList[j][0] == edge[0] && edgeList[j][1] == edge[1]) || (edgeList[j][1] == edge[0] && edgeList[j][0] == edge[1])) {
				has = true;
				break;
			}
		}
		if(!has)
			edgeList.push_back(edge);
	}
	for (int i = 0; i < medialSlabsNum; i++)
	{
		Vector3i face = getMedialSlab(i);
		Vector2i edge[3];
		edge[0] = Vector2i(face[0], face[1]);
		edge[1] = Vector2i(face[0], face[2]);
		edge[2] = Vector2i(face[1], face[2]);
		for (int k = 0; k < 3; k++)
		{
			bool has = false;
			for (int j = 0; j < edgeList.size(); j++)
			{
				if ((edgeList[j][0] == edge[k][0] && edgeList[j][1] == edge[k][1]) || (edgeList[j][1] == edge[k][0] && edgeList[j][0] == edge[k][1])) {
					has = true;
					break;
				}
			}
			if (!has)edgeList.push_back(edge[k]);
		}
	}


	return _valid;
}

bool BaseMedialMesh::writeMeshToMatFormat(const std::string filename)
{
	std::string dir, name, format;
	getFilenameInfo(filename, dir, name, format);
	if (format == std::string("mat") || format == std::string("ma"))
	{
		std::ofstream fout(filename.c_str());
		fout << medialPointsNum << " " << medialConesNum << " " << medialSlabsNum << std::endl;
		for (size_t i = 0; i < medialPointsNum; i++)
		{
			fout << "v";
			for (size_t j = 0; j < 3; j++) fout << " " << medialPoints.buffer[3 * i + j];
			fout << " " << medialRadius.buffer[i];
			fout << " " << bindedTM[i];
			fout << std::endl;
		}
		for (size_t i = 0; i < medialConesNum; i++)
		{
			Vector2i cone = getMedialCone(i);
			fout << "c " << cone.data()[0] <<" " << cone.data()[1] << std::endl;
		}
		for (size_t i = 0; i < medialSlabsNum; i++)
		{
			Vector3i slab = getMedialSlab(i);
			fout << "s " << slab.data()[0] << " " << slab.data()[1] <<" " << slab.data()[2] << std::endl;
		}
		fout.close();
		return true;
	}
	return false;
}


void BaseMedialMesh::uniform(const qeal div)
{
	for (int i = 0; i < medialPointsNum; i++)
	{
		medialPoints.buffer[3 * i] /= div;
		medialPoints.buffer[3 * i + 1] /= div;
		medialPoints.buffer[3 * i + 2] /= div;
		medialRadius.buffer[i] /= div;
	}
}

void BaseMedialMesh::translateMesh(const qeal x, const qeal y, const qeal z)
{
	for (int i = 0; i < medialPointsNum; i++)
	{
		medialPoints.buffer[3 * i] += x;
		medialPoints.buffer[3 * i + 1] += y;
		medialPoints.buffer[3 * i + 2] += z;
	}
}

void BaseMedialMesh::scaleMesh(qeal s, const qeal cx, const qeal cy, const qeal cz)
{
	for (int i = 0; i < medialPointsNum; i++)
	{
		medialPoints.buffer[3 * i] = (medialPoints.buffer[3 * i] - cx) * s + cx;
		medialPoints.buffer[3 * i + 1] = (medialPoints.buffer[3 * i + 1] - cy) * s + cy;
		medialPoints.buffer[3 * i + 2] = (medialPoints.buffer[3 * i + 2] - cz) * s + cz;
		medialRadius.buffer[i] *= s;
	}
}

void BaseMedialMesh::rotateMesh(const qglviewer::Quaternion oldR, const qglviewer::Quaternion R, const qeal cx, const qeal cy, const qeal cz)
{
	for (int i = 0; i < medialPointsNum; i++)
	{
		qglviewer::Vec pos(medialPoints.buffer[3 * i] - cx, medialPoints.buffer[3 * i + 1] - cy, medialPoints.buffer[3 * i + 2] - cz);
		pos = oldR.inverseRotate(pos);
		pos = R.rotate(pos);

		medialPoints.buffer[3 * i] = pos.x + cx;
		medialPoints.buffer[3 * i + 1] = pos.y + cy;
		medialPoints.buffer[3 * i + 2] = pos.z + cz;
	}
}

void BaseMedialMesh::encloseObject(qeal * objPos, int dim)
{
	int num = dim / 3;
	std::vector<qeal> dr(medialPointsNum, 0);
	for (int i = 0; i < num; i++)
	{
		Vector3 p = Vector3(objPos[3 * i], objPos[3 * i + 1], objPos[3 * i + 2]);
		qeal pr = 0.0;

		qeal minDist = QEAL_MAX;
		int minPrimitiveId = -1;
		for (int c = 0; c < medialConesNum; c++)
		{
			Vector2i cone = getMedialCone(c);
			Vector3 c11 = getMedialPoint(cone[0]);
			Vector3 c12 = getMedialPoint(cone[1]);
			qeal r11 = getMedialRadius(cone[0]);
			qeal r12 = getMedialRadius(cone[1]);
			qeal t, dist;
			getNearestSphereOnMedialCone(p, pr, c11, r11, c12, r12, t, dist);
			if (dist < minDist)
			{
				minDist = dist;
				minPrimitiveId = c;
			}
		}

		for (int s = 0; s < medialSlabsNum; s++)
		{
			Vector3i slab = getMedialSlab(s);
			Vector3 c11 = getMedialPoint(slab[0]);
			Vector3 c12 = getMedialPoint(slab[1]);
			Vector3 c13 = getMedialPoint(slab[2]);
			qeal r11 = getMedialRadius(slab[0]);
			qeal r12 = getMedialRadius(slab[1]);
			qeal r13 = getMedialRadius(slab[2]);
			qeal t0, t1, dist;
			getNearestSphereOnMedialSlab(p, pr, c11, r11, c12, r12, c13, r13, t0, t1, dist);
			if (dist < minDist)
			{
				minDist = dist;
				minPrimitiveId = medialConesNum + s;
			}
		}

		if (minDist <= 0) continue;

		Vector3i primitive = getMedialPrimitive(minPrimitiveId);
		for (int j = 0; j < 3; j++)
		{
			int mi = primitive[j];
			if (mi == -1) continue;
			qeal r = minDist * 1.00005;
			if (dr[mi] < r)
				dr[mi] = r;
		}
	}
	for (int i = 0; i < medialPointsNum; i++)
	{
		qeal r = getMedialRadius(i) + dr[i];
		setMedialRadius(i, r);
	}
}


