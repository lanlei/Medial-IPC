#include "MipcModel.h"

namespace MIPC
{
	void MipcModel::alignAllMesh(qeal * tetPointsPtr)
	{
		updateSurfaceFromTetMesh(tetPointsPtr + tetPoints.offset);
	}

	void MipcModel::computeElementSetFromWeightMatrix(MatrixX& weight)
	{
		std::vector<BaseTetElementSet> elementSet;
		for (int i = 0; i < weight.cols(); i++)
		{
			VectorX col = weight.col(i);
			QString frameIdStr;
			frameIdStr.setNum(i);
			std::string setName = "F" + frameIdStr.toStdString() + "_elements_set";
			BaseTetElementSet set(setName);
			for (int j = 0; j < tetPointsNum; j++)
			{
				if (weight.data()[i * weight.rows() + j] == 0) continue;
				for (int k = 0; k < tetPointElementList[j].span; k++)
				{
					int eleid = tetPointElementList[j].buffer[k];
					set.insert(eleid);
				}
			}
			elementSet.push_back(set);
		}
		_tmHandle->setTetMeshSet(elementSet);
	}

	bool MipcModel::readFrameMatList(std::string filename)
	{
		_reducedTypes.resize(medialPointsNum, ReducedFrameType::STATIC);

		if (!isMedialMeshValid())
			return false;

		std::ifstream fin(filename.c_str());
		if (fin.is_open())
		{
			int linearNum, quadNum, transNum;
			fin >> linearNum >> quadNum >> transNum;

			_staticFramesNum = medialPointsNum - linearNum - quadNum - transNum;
			_linearFramesNum = linearNum;
			_quadFramesNum = quadNum;
			_translaFramesNum = transNum;

			for (size_t i = 0; i < linearNum; i++)
			{
				int id;
				fin >> id;
				_matFrameId.push_back(id);
				_reducedTypes[id] = ReducedFrameType::LINEAR;
			}
			for (size_t i = 0; i < quadNum; i++)
			{
				int id;
				fin >> id;
				_matFrameId.push_back(id);
				_reducedTypes[id] = ReducedFrameType::Quadratic;
			}
			for (size_t i = 0; i < transNum; i++)
			{
				int id;
				fin >> id;
				_matFrameId.push_back(id);
				_reducedTypes[id] = ReducedFrameType::Translation;
			}
		}
		else
		{
			int linearNum = medialPointsNum;
			int quadNum = 0, transNum = 0;
			_staticFramesNum = medialPointsNum - linearNum - quadNum - transNum;
			_linearFramesNum = linearNum;
			_quadFramesNum = quadNum;
			_translaFramesNum = transNum;

			_matFrameId.resize(linearNum);
			for (size_t i = 0; i < _reducedTypes.size(); i++)
				_reducedTypes[i] = ReducedFrameType::LINEAR;
			for (size_t i = 0; i < linearNum; i++)
			{
				_matFrameId[i] = i;
			}
		}

		_reducedDim = 12 * _linearFramesNum + 30 * _quadFramesNum + 3 * _translaFramesNum;

		return true;
	}

	bool MipcModel::writeFrameMatList(std::string filename)
	{
		if (!isMedialMeshValid())
			return false;

		std::ofstream fout(filename.c_str());
		std::vector<int> linearList;
		std::vector<int> quadList;
		for (size_t i = 0; i < _reducedTypes.size(); i++)
		{
			int mid = _matFrameId[i];
			if (_reducedTypes[i] == ReducedFrameType::LINEAR)
			{
				linearList.push_back(mid);
			}
			else if (_reducedTypes[i] == ReducedFrameType::Quadratic)
			{
				quadList.push_back(mid);
			}
		}
		fout << linearList.size() << " " << quadList.size() << std::endl;
		for (size_t i = 0; i < linearList.size(); i++)
			fout << linearList[i] << std::endl;
		for (size_t i = 0; i < quadList.size(); i++)
			fout << quadList[i] << std::endl;
		fout.close();
		return true;
	}

	void MipcModel::createReducedFrame(int& frameOffset, int& bufferOffset, qeal * X, qeal * preX, qeal * Xtilde, qeal * Vel, qeal * preVel, qeal * Acc, qeal * preAcc, std::vector<ReducedFrame*>& frameList)
	{
		for (size_t i = 0; i < medialPointsNum; i++)
		{
			ReducedFrameType frameType = _reducedTypes[i];
			qeal* p = medialPoints.buffer + 3 * i;
			if (frameType == ReducedFrameType::STATIC)
			{
				ReducedFrame* frame = new ReducedFrame(p, -1);
				frameList.push_back(frame);
				_reducedFrames.push_back(frame);
			}
			else if (frameType == ReducedFrameType::LINEAR)
			{
				LinearReducedFrame* frame = new LinearReducedFrame(frameOffset, p, X, preX, Xtilde, Vel, preVel, Acc, preAcc, bufferOffset);
				frameList.push_back(frame);
				_reducedFrames.push_back(frame);
				frameOffset += 1;
				bufferOffset += 12;
			}
			else if (frameType == ReducedFrameType::Quadratic)
			{
				QuadraticReducedFrame* frame = new QuadraticReducedFrame(frameOffset, p, X, preX, Xtilde, Vel, preVel, Acc, preAcc, bufferOffset);
				frameList.push_back(frame);
				_reducedFrames.push_back(frame);
				frameOffset += 1;
				bufferOffset += 30;
			}
			else if (frameType == ReducedFrameType::Translation)
			{
				TranslationReducedFrame* frame = new TranslationReducedFrame(frameOffset, p, X, preX, Xtilde, Vel, preVel, Acc, preAcc, bufferOffset);
				frameList.push_back(frame);
				_reducedFrames.push_back(frame);
				frameOffset += 1;
				bufferOffset += 3;
			}
		}
	}

	void MipcModel::computeLaplacianMatrix(std::vector<TripletX>& matValue)
	{
		for (size_t i = 0; i < tetPointsNum; i++)
		{
			matValue.push_back(Eigen::Triplet<double>(i, i, 1));
			qeal edgeLenTotal = 0.0;
			Vector3 pos = getTetPoint(i);
			int* neighbor = tetPointNeighborList[i].buffer;
			for (int it = 0; it < tetPointNeighborList[i].span; it++)
			{
				Vector3 np = getTetPoint(neighbor[it]);
				qeal edgeLen = (pos - np).norm();
				edgeLenTotal += edgeLen;
			}
			for (int it = 0; it < tetPointNeighborList[i].span; it++)
			{
				Vector3 np = getTetPoint(neighbor[it]);
				qeal edgeLen = (pos - np).norm();
				qeal wij = edgeLen / edgeLenTotal;
				matValue.push_back(TripletX(i, neighbor[it], -wij));
			}
		}
	}

	MatrixX MipcModel::computeHarmonicWeight()
	{
		SparseMatrix coeffMat;
		int tn = tetPointsNum;
		int fn = getNonStaticFramesNum();
		int rowDim = tn + fn;
		int colDim = tn;

		std::vector<TripletX> matValue;
		computeLaplacianMatrix(matValue);

		int rowIndex = tn;
		for (size_t frameId = 0; frameId < _matFrameId.size(); frameId++)
		{
			int mid = _matFrameId[frameId];
			int tetId = bindedTM[mid];
			matValue.push_back(TripletX(rowIndex++, tetId, 1));
		}

		//
		std::string filename = dir + "fixed_tet_points.txt";
		std::ifstream fin(filename.c_str());
		int fixedTetNodesNum = 0;
		std::set<int> cs;
		if (fin.is_open())
		{
			char ch;
			fin >> ch >> fixedTetNodesNum;
			if (ch != 'c')
			{
				fixedTetNodesNum = 0;
			}
			else
			{
				//cs.resize(fixedTetNodesNum);
				for (int j = 0; j < fixedTetNodesNum; j++)
				{
					int tid;
					fin >> tid;
					cs.insert(tid);
				}
			}
		}
		std::vector<int> ccs;
		std::set<int>::iterator it = cs.begin();
		for (; it != cs.end(); ++it)
		{
			ccs.push_back(*it);
		}
		fixedTetNodesNum = ccs.size();
		for (int i = 0; i < fixedTetNodesNum; i++)
		{
			matValue.push_back(TripletX(rowDim + i, ccs[i], 1));
		}

		rowDim += fixedTetNodesNum;
		coeffMat.resize(rowDim, colDim);
		coeffMat.setFromTriplets(matValue.begin(), matValue.end());

		MatrixX weight;
		weight.resize(tn, fn);
		weight.setZero();
		SparseMatrix transpose = coeffMat.transpose();
		SparseMatrix S = transpose * coeffMat;

		Eigen::SimplicialLDLT<SparseMatrix> LDLT(S);
		for (int i = 0; i < fn; i++)
		{
			VectorX d = transpose.col(tn + i);
			VectorX w = LDLT.solve(d);
			weight.col(i) = w;
		}

		if (setWeightLocality)
		{
			for (int i = 0; i < tetPointsNum; i++)
			{
				VectorX ws = weight.row(i);
				qeal sum = 0;
				for (int j = 0; j < ws.size(); j++)
				{
					if (abs(ws.data()[j]) < weightLocalityThreshold)
					{
						ws.data()[j] = 0;
					}

					sum += ws.data()[j];
				}
				for (int j = 0; j < ws.size(); j++)
				{
					if (abs(sum) < 1e-13)
						continue;
					ws.data()[j] = ws.data()[j] / sum;
				}
				weight.row(i) = ws;
			}
			computeElementSetFromWeightMatrix(weight);
		}
		return weight;

	}

	void MipcModel::getGlobalSparseProjectionMatrixTriplet(int rOffset, int cOffset, std::vector<TripletX>& globalTriplet)
	{
		assert(_matFrameId.size() == _reducedFrames.size() && _matFrameId.size() == _reducedTypes.size());
		assert(_matFrameId.size() == (_linearFramesNum + _quadFramesNum + _translaFramesNum));

		size_t tn = tetPointsNum;
		size_t rows = 3 * tn;
		size_t cols = _reducedDim;

		_reducedProjection.resize(rows, cols);

		_tetPointShareFramesList.resize(tetPointsNum);
		_tetElementShareFramesList.resize(tetElementNum);
		_frameShareTetPointsList.resize(_reducedFrames.size());
		_frameShareTetElementList.resize(_reducedFrames.size());

		_harmonicWeight = computeHarmonicWeight();

		for (size_t i = 0; i < tn; i++)
		{
			Vector3 p = getTetPoint(i);
			int offset = 0;
			for (size_t it = 0; it < _reducedFrames.size(); it++)
			{
				qeal w = _harmonicWeight.data()[it * _harmonicWeight.rows() + i];
				if (IS_QEAL_ZERO(w))
				{
					offset += 12;
					continue;
				}
				_tetPointShareFramesList[i].insert(it);
				_frameShareTetPointsList[it].insert(i);

				for (int j = 0; j < tetPointElementList[i].span; j++)
				{
					int eleId = tetPointElementList[i].buffer[j];
					_tetElementShareFramesList[eleId].insert(it);
					_frameShareTetElementList[it].insert(eleId);
				}

				MatrixX u = _reducedFrames[it]->getUMatrix(p.data(), w);
				if (_reducedTypes[it] == ReducedFrameType::LINEAR)
				{

					_reducedProjection.block(3 * i, offset, 3, 12) = u;
					globalTriplet.push_back(TripletX(rOffset + 3 * i + 0, cOffset + offset + 0, u.data()[0]));
					globalTriplet.push_back(TripletX(rOffset + 3 * i + 1, cOffset + offset + 1, u.data()[4]));
					globalTriplet.push_back(TripletX(rOffset + 3 * i + 2, cOffset + offset + 2, u.data()[8]));

					globalTriplet.push_back(TripletX(rOffset + 3 * i + 0, cOffset + offset + 3, u.data()[9]));
					globalTriplet.push_back(TripletX(rOffset + 3 * i + 1, cOffset + offset + 4, u.data()[13]));
					globalTriplet.push_back(TripletX(rOffset + 3 * i + 2, cOffset + offset + 5, u.data()[17]));

					globalTriplet.push_back(TripletX(rOffset + 3 * i + 0, cOffset + offset + 6, u.data()[18]));
					globalTriplet.push_back(TripletX(rOffset + 3 * i + 1, cOffset + offset + 7, u.data()[22]));
					globalTriplet.push_back(TripletX(rOffset + 3 * i + 2, cOffset + offset + 8, u.data()[26]));

					globalTriplet.push_back(TripletX(rOffset + 3 * i + 0, cOffset + offset + 9, u.data()[27]));
					globalTriplet.push_back(TripletX(rOffset + 3 * i + 1, cOffset + offset + 10, u.data()[31]));
					globalTriplet.push_back(TripletX(rOffset + 3 * i + 2, cOffset + offset + 11, u.data()[35]));

					offset += 12;
				}
			}
		}


	}

	MatrixX MipcModel::getReducedProjection()
	{
		MatrixX pro;

		size_t tn = tetPointsNum;
		size_t rows = 3 * tn;
		size_t cols = _reducedDim;
		pro.resize(rows, cols);
		for (size_t i = 0; i < tn; i++)
		{
			Vector3 p = getTetPoint(i);
			int offset = 0;
			for (size_t it = 0; it < _reducedFrames.size(); it++)
			{
				qeal w = _harmonicWeight.data()[it * _harmonicWeight.rows() + i];
				if (IS_QEAL_ZERO(w))
				{
					offset += 12;
					continue;
				}
				MatrixX u = _reducedFrames[it]->getUMatrix(p.data(), w);
				if (_reducedTypes[it] == ReducedFrameType::LINEAR)
				{
					pro.block(3 * i, offset, 3, 12) = u;
					offset += 12;
				}
			}
		}

		return pro;
	}

}