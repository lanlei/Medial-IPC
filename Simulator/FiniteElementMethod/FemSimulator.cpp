#include "FemSimulator.h"
#include "Commom\SPDProjectFunction.h"
#include <QElapsedTimer>

namespace FiniteElementMethod
{
	bool FemSimulator::addModelFromConfigFile(const std::string filename, TiXmlElement * item)
	{
		FemModel* m = new FemModel();
		bool isReadMesh = BaseSimulator::addModelFromConfigFile(filename, item, m);
		if (!isReadMesh)
			return false;
		m = getModel(models.size() - 1);

		qeal scale = 1.0;
		qeal tx = 0.0, ty = 0.0, tz = 0.0;
		qeal rx = 0.0, ry = 1.0, rz = 0.0, rsita = 0.0;
		qeal density = 1.0;
		qeal youngModul = 300.0, poisson = 0.4;
		int enableGravity = 1;
		TiXmlElement* childItem = item->FirstChildElement();
		std::strstream ss;
		while (childItem)
		{
			ss.clear();
			std::string itemName = childItem->Value();

			if (itemName == std::string("scale"))
			{
				std::string str = childItem->GetText();
				ss << str;
				ss >> scale;
			}
			else if (itemName == std::string("translation"))
			{
				std::string str = childItem->GetText();
				ss << str;
				ss >> tx >> ty >> tz;
			}
			else if (itemName == std::string("rotation"))
			{
				std::string str = childItem->GetText();
				ss << str;
				ss >> rx >> ry >> rz >> rsita;
			}
			else if (itemName == std::string("density"))
			{
				std::string str = childItem->GetText();
				ss << str;
				ss >> density;
			}
			else if (itemName == std::string("YoungModul"))
			{
				std::string str = childItem->GetText();
				ss << str;
				ss >> youngModul;
			}
			else if (itemName == std::string("Poisson"))
			{
				std::string str = childItem->GetText();
				ss << str;
				ss >> poisson;
			}
			else if (itemName == std::string("enableGravity"))
			{
				std::string str = childItem->GetText();
				ss << str;
				ss >> enableGravity;
			}

			childItem = childItem->NextSiblingElement();
		}

		m->scaleModel(scale);
		m->translateModel(tx, ty, tz);
		m->rotateModel(rx, ry, rz, rsita);
		m->enableGravityForce(enableGravity);
		m->initMeshesHandel();
		ENuMaterial * material = downcastENuMaterial(m->getTetMeshHandle()->getElementMaterialById(0));
		material->setDensity(density);
		material->setE(youngModul);
		material->setNu(poisson);

		return true;
	}

	bool FemSimulator::addStaticModelFromConfigFile(const std::string filename, TiXmlElement * item)
	{
		BaseModel* m = new BaseModel();
		bool isReadMesh = BaseSimulator::addStaticModelFromConfigFile(filename, item, m);
		if (!isReadMesh)
			return false;
		m = getStaticModel(staticModels.size() - 1);

		qeal scale = 1.0;
		qeal tx = 0.0, ty = 0.0, tz = 0.0;
		qeal rx = 0.0, ry = 1.0, rz = 0.0, rsita = 0.0;

		TiXmlElement* childItem = item->FirstChildElement();
		std::strstream ss;
		while (childItem)
		{
			ss.clear();
			std::string itemName = childItem->Value();
			if (itemName == std::string("scale"))
			{
				std::string str = childItem->GetText();
				ss << str;
				ss >> scale;
			}
			else if (itemName == std::string("translation"))
			{
				std::string str = childItem->GetText();
				ss << str;
				ss >> tx >> ty >> tz;
			}
			else if (itemName == std::string("rotation"))
			{
				std::string str = childItem->GetText();
				ss << str;
				ss >> rx >> ry >> rz >> rsita;
			}
		
			childItem = childItem->NextSiblingElement();
		}

		m->scaleModel(scale);
		m->translateModel(tx, ty, tz);
		m->rotateModel(rx, ry, rz, rsita);
		m->initMeshesHandel();
		return true;
	}

	void FemSimulator::nullAllPointer()
	{

	}

	void FemSimulator::deleteAllPointer()
	{
		for (int i = 0; i < totalTetElementNum; i++)
			free(_row[i]);
		free(_row);

		for (int i = 0; i < totalTetElementNum; i++)
			free(_column[i]);
		free(_column);

		free(solverConfig);
	}

	void FemSimulator::saveFile()
	{
	}

	void FemSimulator::saveSimulator()
	{
	}

	void FemSimulator::initialization()
	{
		Eigen::initParallel();
		omp_set_num_threads(8);
		
		_sysDim = totalTetPointsNum * 3;
		_sysXtilde.resize(_sysDim);
		_sysXtilde.setZero();
		_sysDir.resize(_sysDim);
		_sysDir.setZero();
		_sysMatrix.resize(_sysDim, _sysDim);
		_sysMassMatrix.resize(_sysDim, _sysDim);
		_sysStiffnessMatrix.resize(_sysDim, _sysDim);
		_sysDampingMatrix.resize(_sysDim, _sysDim);

		_sysRhs.resize(_sysDim);
		_sysRhs.setZero();
		_sysX.resize(_sysDim);
		_sysX.setZero();
		_sysXtilde.resize(_sysDim);
		_sysXtilde.setZero();

		_sysXn.resize(_sysDim);
		_sysXn.setZero();
		_sysVn.resize(_sysDim);
		_sysVn.setZero();
		_sysAccn.resize(_sysDim);
		_sysAccn.setZero();

		_sysSurfaceOriginalPosition.resize(3 * totalPointsNum);
		std::copy(pointsBuffer.buffer.begin(), pointsBuffer.buffer.end(), _sysSurfaceOriginalPosition.data());

		_sysOriginalPosition.resize(_sysDim);
		std::copy(tetPointsBuffer.buffer.begin(), tetPointsBuffer.buffer.end(), _sysOriginalPosition.data());
		_sysCurrentPosition.resize(_sysDim);
		_sysCurrentPosition = _sysOriginalPosition;

		_sysVelocity.resize(_sysDim);
		_sysVelocity.setZero();
		_sysAccrelation.resize(_sysDim);
		_sysAccrelation.setZero();
		_sysExternalForce.resize(_sysDim);
		_sysExternalForce.setZero();
		_sysInternalForce.resize(_sysDim);
		_sysInternalForce.setZero();
		_sysMouseForce.resize(_sysDim);
		_sysMouseForce.setZero();
		_sysCollideForce.resize(_sysDim);
		_sysCollideForce.setZero();
		_meshCollisionForce.resize(3 * totalPointsNum);
		dPdF.resize(totalTetElementNum);
		for (size_t i = 0; i < dPdF.size(); i++)
			dPdF[i].resize(9, 9);
		elementsInternalForce.resize(totalTetElementNum);

		_sysMassMatrix = computeMassMatrix();
		_sysInverseMassMatrix = computeInverseMassMatrix();
		computeGravityForce();
		preFillSysStiffnessMatrix();
		_sysMatrix = _sysMassMatrix + _sysStiffnessMatrix;
		handleTetPointsConstraint();

		solverConfig->updateAlphas(_timeStep);

		//
		_tetElementInternalForce.resize(12 * totalTetElementNum);
		_tetElementInternalForce.setZero();

		computePartialPF(_sysX);

		_maxBbox = Vector3(-DBL_MAX, -DBL_MAX, -DBL_MAX);
		_minBbox = Vector3(DBL_MAX, DBL_MAX, DBL_MAX);

		for (int i = 0; i < totalTetPointsNum; i++)
		{
			qeal x = tetPointsBuffer.buffer[3 * i];
			qeal y = tetPointsBuffer.buffer[3 * i + 1];
			qeal z = tetPointsBuffer.buffer[3 * i + 2];
			if (_maxBbox[0] < x) _maxBbox[0] = x;
			if (_maxBbox[1] < y) _maxBbox[1] = y;
			if (_maxBbox[2] < z) _maxBbox[2] = z;
			if (_minBbox[0] > x) _minBbox[0] = x;
			if (_minBbox[1] > y) _minBbox[1] = y;
			if (_minBbox[2] > z) _minBbox[2] = z;
		}

		_diagLen = (_maxBbox - _minBbox).norm();

	}

	void FemSimulator::run(int frame)
	{
		doTime2(frame);
	}

	void FemSimulator::postRun()
	{
		_sysCurrentPosition = _sysOriginalPosition + _sysX;
		std::copy(_sysCurrentPosition.data(), _sysCurrentPosition.data() + _sysCurrentPosition.size(), tetPointsBuffer.buffer.data());

		alignAllMesh(tetPointsBuffer.buffer.data());
	}

	void FemSimulator::reset()
	{
		_sysCurrentPosition = _sysOriginalPosition;	
		std::copy(_sysSurfaceOriginalPosition.data(), _sysSurfaceOriginalPosition.data() + _sysSurfaceOriginalPosition.size(), pointsBuffer.buffer.data());

		std::copy(_sysOriginalPosition.data(), _sysOriginalPosition.data() + _sysOriginalPosition.size(), tetPointsBuffer.buffer.data());

		_sysRhs.setZero();
		_sysX.setZero();
		_sysVelocity.setZero();
		_sysAccrelation.setZero();
		_sysXn.setZero();
		_sysVn.setZero();
		_sysAccn.setZero();
		_sysExternalForce.setZero();
		_sysInternalForce.setZero();
		_sysMouseForce.setZero();
		_sysCollideForce.setZero();
		computeGravityForce();
	}

	void FemSimulator::handleSurfacePointsSelectedEvent()
	{
	}

	void FemSimulator::handleTetPointsSelectedEvent()
	{
	}

	void FemSimulator::addTetPointsConstraint()
	{
		_constrainedTetPointsNum = 0;
		_constrainedTetPoints.clear();
		for (int i = 0; i < models.size(); i++)
		{
			std::string filename = getModel(i)->dir + "fixed_tet_points.txt";
			std::ifstream fin(filename.c_str());
			if (!fin.is_open())
				continue;
			int num = 0;
			char ch;
			fin >> ch >> num;
			if (ch != 'c')
				continue;

			if (num > 0)
			{
				for (int j = 0; j < num; j++)
				{
					int tid;
					fin >> tid;
					_constrainedTetPoints.insert(getModel(i)->getTetPointOverallId(tid));
				}
			}
			else if (num < 0)
			{
				for (int j = 0; j < getModel(i)->tetPointsNum; j++)
				{
					_constrainedTetPoints.insert(getModel(i)->getTetPointOverallId(j));
				}
			}
		}
		_constrainedTetPointsNum = _constrainedTetPoints.size();
		int sysSubDim = 3 * (totalTetPointsNum - _constrainedTetPointsNum);
		_sysSubMatrix.resize(sysSubDim, sysSubDim);
		_sysSubRhs.resize(sysSubDim);
		_sysSubXtilde.resize(sysSubDim);
	}

	void FemSimulator::createMapByConstraints()
	{
		_mapOldNew.clear();
		_mapOldNew.resize(_sysDim);
		std::fill(_mapOldNew.begin(), _mapOldNew.end(), 0);

		int newIndex = 0;

		for (int i = 0; i < totalTetPointsNum; i++)
		{
			if (_constrainedTetPoints.find(i) != _constrainedTetPoints.end())
			{
				_mapOldNew[3 * i] = -1;
				_mapOldNew[3 * i + 1] = -1;
				_mapOldNew[3 * i + 2] = -1;
			}
			else
			{
				_mapOldNew[3 * i] = newIndex++;
				_mapOldNew[3 * i + 1] = newIndex++;
				_mapOldNew[3 * i + 2] = newIndex++;
			}
		}
	}

	void FemSimulator::preFillSysStiffnessMatrix()
	{
		_row = (int **)malloc(sizeof(int*)* totalTetElementNum);
		_column = (int **)malloc(sizeof(int*)*totalTetElementNum);

		_stiffnessMatrixTopology = new SparseMatrix;
		getSparseMatrixTopologyTYPE(*_stiffnessMatrixTopology);
		SparseMatrixTopologyTYPE<qeal>* sparseMatrixTopology = new SparseMatrixTopologyTYPE<qeal>(_stiffnessMatrixTopology);

		for (size_t mid = 0; mid < models.size(); mid++)
		{
			BaseModel* m = getModel(mid);
			for (size_t eid = 0; eid < m->tetElementNum; eid++)
			{
				int geid = m->getTetElementOverallId(eid);
				_row[geid] = (int*)malloc(sizeof(int) * 4);
				_column[geid] = (int*)malloc(sizeof(int) * 4 * 4 * 3 * 3);

				Vector4i indices = m->getTetElement(eid);
				for (size_t ver = 0; ver < 4; ver++)
					_row[geid][ver] = m->getTetPointOverallId(indices.data()[ver]);

				for (size_t i = 0; i < 4; i++)
					for (size_t j = 0; j < 4; j++)
					{
						for (size_t k = 0; k < 3; k++)
						{
							for (size_t l = 0; l < 3; l++)
							{
								int block_r = i * 3 + k;
								int block_c = j * 3 + l;
								_column[geid][3 * 4 * block_c + block_r] = sparseMatrixTopology->getValueIndex(3 * _row[geid][i] + k, 3 * _row[geid][j] + l);
							}
						}
					}
			}
		}

		getSparseMatrixTopologyTYPE(_sysStiffnessMatrix);
	}

	void FemSimulator::handleTetPointsConstraint()
	{
		addTetPointsConstraint();
		createMapByConstraints();
		createSparseMapbyTopology(_sysMatrix, _sysSubMatrix);
		sparseMatrixRemoveRows(_sysMatrix, _sysSubMatrix);
		createSparseMapbyTopology(_sysMassMatrix, _sysSubMassMatrix);
		sparseMatrixRemoveRows(_sysMassMatrix, _sysSubMassMatrix);
		createSparseMapbyTopology(_sysStiffnessMatrix, _sysSubStiffnessMatrix);
	}

	void FemSimulator::computeGravityForce()
	{
		_sysGravity.resize(_sysDim);
		_sysGravityForce.resize(_sysDim);
		_sysGravityForce.setZero();
		_sysGravity.setZero();
		for (int i = 0; i < models.size(); i++)
		{
			qeal gravity = _gravity;
			if (!models[i]->enableGravity)
				gravity = 0;
			for (int j = 0; j < models[i]->tetPointsNum; j++)
			{
				int gid = models[i]->getTetPointOverallId(j);
				_sysGravity.data()[3 * gid + 1] = gravity;
			}
		}
		_sysGravityForce = _sysMassMatrix * _sysGravity;
	}

	void FemSimulator::computeExternalForce()
	{
		_sysExternalForce.setZero();
		if (isUseGravity()) _sysExternalForce += _sysGravityForce;

		_sysExternalForce += _sysCollideForce;
		_sysCollideForce.setZero();
		_sysExternalForce += _sysMouseForce;
		_sysMouseForce.setZero();
	}

	void FemSimulator::computeInternalForce()
	{
		_sysInternalForce.setZero();
		for (size_t i = 0; i < models.size(); i++)
		{
			BaseModel* m = getModel(i);
			for (size_t eid = 0; eid < m->tetElementNum; eid++) 
			{
				int geid = m->getTetElementOverallId(eid);
				Vector4i indices = m->getTetElementOverall(eid);
				Vector3 f[4];
				for (int j = 0; j < 4; j++)
				{
					int nid = indices.data()[j];
					if (j != 3)
					{
						f[j].data()[0] = elementsInternalForce[geid].data()[j * 3];
						f[j].data()[1] = elementsInternalForce[geid].data()[j * 3 + 1];
						f[j].data()[2] = elementsInternalForce[geid].data()[j * 3 + 2];
					}
					else
					{
						f[j].data()[0] = -1.0*(f[0].data()[0] + f[1].data()[0] + f[2].data()[0]);
						f[j].data()[1] = -1.0*(f[0].data()[1] + f[1].data()[1] + f[2].data()[1]);
						f[j].data()[2] = -1.0*(f[0].data()[2] + f[1].data()[2] + f[2].data()[2]);
					}					
					_sysInternalForce.data()[3 * nid] -= f[j].data()[0];
					_sysInternalForce.data()[3 * nid + 1] -= f[j].data()[1];
					_sysInternalForce.data()[3 * nid + 2] -= f[j].data()[2];
				}

				//_tetElementInternalForce[12 * geid + 0] = -elementsInternalForce[geid].data()[0];
				//_tetElementInternalForce[12 * geid + 1] = -elementsInternalForce[geid].data()[1];
				//_tetElementInternalForce[12 * geid + 2] = -elementsInternalForce[geid].data()[2];

				//_tetElementInternalForce[12 * geid + 3] = -elementsInternalForce[geid].data()[3];
				//_tetElementInternalForce[12 * geid + 4] = -elementsInternalForce[geid].data()[4];
				//_tetElementInternalForce[12 * geid + 5] = -elementsInternalForce[geid].data()[5];

				//_tetElementInternalForce[12 * geid + 6] = -elementsInternalForce[geid].data()[6];
				//_tetElementInternalForce[12 * geid + 7] = -elementsInternalForce[geid].data()[7];
				//_tetElementInternalForce[12 * geid + 8] = -elementsInternalForce[geid].data()[8];

				//_tetElementInternalForce[12 * geid + 9] = (elementsInternalForce[geid].data()[0] + elementsInternalForce[geid].data()[3] + elementsInternalForce[geid].data()[6]);
				//_tetElementInternalForce[12 * geid + 10] = (elementsInternalForce[geid].data()[1] + elementsInternalForce[geid].data()[4] + elementsInternalForce[geid].data()[7]);
				//_tetElementInternalForce[12 * geid + 11] = (elementsInternalForce[geid].data()[2] + elementsInternalForce[geid].data()[5] + elementsInternalForce[geid].data()[8]);

				//if (geid == 1000)
				//{
				//	for (int c = 0; c < 12; c++)
				//		printf("cpu %d, %f\n", c, _tetElementInternalForce[12 * geid + c]);
				//}

			}
		}
	}

	void FemSimulator::computeStiffnessMatrix()
	{
		resetSparseMatrix(_sysStiffnessMatrix);
		for (size_t i = 0; i < models.size(); i++)
		{
			BaseModel* m = getModel(i);
			for (size_t eid = 0; eid < m->tetElementNum; eid++) 
			{
				int geid = m->getTetElementOverallId(eid);
				Vector4i indices = m->getTetElementOverall(eid);

				MatrixX Kij;
				MatrixX dPdu = dPdF[geid] * m->getTetMeshHandle()->getTetElementParam(eid)->dFdu;
				Kij = dPdu.transpose() * m->getTetMeshHandle()->getTetElementParam(eid)->dFdu;
				Kij *= m->getTetMeshHandle()->getTetElementParam(eid)->volume;

				//makePD(Kij);

				for (int c = 0; c < 4; c++)
				{
					for (int a = 0; a < 4; a++)
					{
						int row = indices.data()[c] * 3;
						int col = indices.data()[a] * 3;

						int vid_i = indices.data()[c];
						int vid_j = indices.data()[a];

						for (int r = 0; r < 3; r++)
						{
							for (int l = 0; l < 3; l++)
							{
								int block_r = c * 3 + r;
								int block_c = a * 3 + l;
								int dd = block_c * 12 + block_r;
								qeal data = Kij.data()[block_c * 12 + block_r];
								int index = _column[geid][3 * 4 * block_c + block_r];
								_sysStiffnessMatrix.valuePtr()[index] += data;
							}
						}
					}
				}
			}
		}
	}

	void FemSimulator::computePartialPF(VectorX& displacement)
	{
		for (size_t i = 0; i < models.size(); i++)
		{
			BaseModel* m = getModel(i);
			for (size_t eid = 0; eid < m->tetElementNum; eid++)
			{
				int geid = m->getTetElementOverallId(eid);
				Vector4i indices = m->getTetElementOverall(eid);
				BaseTetElementParam* para = m->getTetMeshHandle()->getTetElementParam(eid);
				Matrix3 Ds;
				para->computeElementDeformationByShapeMatrix(Ds, displacement.data(), indices.data());
				Matrix3 F = Ds * para->invDm;

				//Matrix3 mU, mV, mSingularF;
				//computeSVD(F, mU, mV, mSingularF, 1e-8, 1);
				//for (int i = 0; i < 3; i++)
				//{
				//	if (mSingularF.data()[i * 3 + i] < 1e-14)
				//	{
				//		mSingularF.data()[i * 3 + i] = 1e-14;
				//	}
				//}
				//F = mU * mSingularF * mV.transpose();

				Matrix3 FPK;
				stableNeoHookean(eid, m->getTetMeshHandle(), F, FPK, dPdF[geid]);
				//commonNeoHookean(eid, m->getTetMeshHandle(), F, FPK, dPdF[geid]);
				elementsInternalForce[geid] = (-1.0 * para->volume) * FPK * para->invDmT;


			}
		}
	}

	void FemSimulator::computeElasticsEnergy(qeal & energy, VectorX& X, VectorX& V)
	{
		energy = 0.0;
		for (size_t i = 0; i < models.size(); i++)
		{
			BaseModel* m = getModel(i);

			for (size_t eid = 0; eid < m->tetElementNum; eid++)
			{
				int geid = m->getTetElementOverallId(eid);
				Vector4i indices = m->getTetElementOverall(eid);
				BaseTetElementParam* para = m->getTetMeshHandle()->getTetElementParam(eid);
				Matrix3 Ds;
				VectorX displacement = X - _sysOriginalPosition;
				para->computeElementDeformationByShapeMatrix(Ds, displacement.data(), indices.data());
				Matrix3 F = Ds * para->invDm;
				qeal e = computeNeoHookeanEnergy(eid, m->getTetMeshHandle(), F);
				energy += _timeStep * _timeStep * e;
			}
		}

		VectorX Sn = _sysCurrentPosition + _timeStep * V + _timeStep * _timeStep * _sysInverseMassMatrix * _sysExternalForce;

		energy += (X - Sn).transpose() * _sysMassMatrix * (X - Sn);

		energy += X.transpose() * (_timeStep * _timeStep *_sysDampingMatrix * _sysVelocity);
	}

	void FemSimulator::computeElasticsEnergy(qeal & energy, SparseMatrix& A, VectorX & rhs, VectorX & X)
	{
		energy = X.transpose() * A * X;
		energy *= 0.5;
		energy -= X.transpose() * rhs;
	}

	void FemSimulator::computeNewtonSolverEnergy(qeal & E, qeal & Eprev, SparseMatrix & A, VectorX & rhs, VectorX& x)
	{
		// Ax = rhs
		// E = Eprev + x^T (-rhs) + 0.5 * x^T A x; 
		E = Eprev + x.transpose() * (0.5 * A * x - rhs);
	}

	void FemSimulator::computeNewtonSolverEnergy(qeal & E, qeal & Eprev, MatrixX & A, VectorX & rhs, VectorX & x, qeal linearSearch)
	{
		// Ax = rhs
		// E = Eprev + x^T (-rhs) + 0.5 * x^T A x; 
		E = Eprev + (x.transpose() * (linearSearch * (linearSearch * 0.5 * A * x - rhs)));
	}

	void FemSimulator::computeNewtonSolverEnergy(qeal & E, MatrixX & A, VectorX & rhs, VectorX & x)
	{
		E = x.transpose() * (0.5 * A * x - rhs);
	}

	qeal FemSimulator::computeNeoHookeanEnergy(int eid, BaseTetMeshHandle* tet, Matrix3& F)
	{
		qeal F11 = F.data()[0];  qeal F12 = F.data()[3];  qeal F13 = F.data()[6];
		qeal F21 = F.data()[1];  qeal F22 = F.data()[4];  qeal F23 = F.data()[7];
		qeal F31 = F.data()[2];  qeal F32 = F.data()[5];  qeal F33 = F.data()[8];
		
		qeal Ic = pow(F11, 2) + pow(F12, 2) + pow(F13, 2) +
			pow(F21, 2) + pow(F22, 2) + pow(F23, 2) +
			pow(F31, 2) + pow(F32, 2) + pow(F33, 2);
		qeal J = F.determinant();
		ENuMaterial * material = downcastENuMaterial(tet->getElementMaterial(eid));
		qeal E = material->getE();
		qeal Nu = material->getNu();
		
		qeal lameMiu = material->getMu();

		qeal lameLamda = material->getLambda();

		qeal miu = (4.0 * lameMiu) / 3.0;
		qeal lamda = lameLamda + (5.0 * lameMiu) / 6.0;

		qeal alpha = 1.0 + (3.0 * miu) / (4.0 * lamda);

		qeal e = 0.5 * miu * (Ic - 3.0) + 0.5 * lamda * (J - alpha) *  (J - alpha) - 0.5 * miu * log(Ic + 1);

		return e;
	}

	qeal FemSimulator::computeNeoHookeanEnergy(VectorX & displacement)
	{
		qeal potential = 0.0;
		for (size_t i = 0; i < models.size(); i++)
		{
			BaseModel* m = getModel(i);
			for (size_t eid = 0; eid < m->tetElementNum; eid++)
			{
				int geid = m->getTetElementOverallId(eid);
				Vector4i indices = m->getTetElementOverall(eid);
				BaseTetElementParam* para = m->getTetMeshHandle()->getTetElementParam(eid);
				Matrix3 Ds;
				para->computeElementDeformationByShapeMatrix(Ds, displacement.data(), indices.data());
				Matrix3 F = Ds * para->invDm;

				qeal F11 = F.data()[0];  qeal F12 = F.data()[3];  qeal F13 = F.data()[6];
				qeal F21 = F.data()[1];  qeal F22 = F.data()[4];  qeal F23 = F.data()[7];
				qeal F31 = F.data()[2];  qeal F32 = F.data()[5];  qeal F33 = F.data()[8];

				qeal Ic = pow(F11, 2) + pow(F12, 2) + pow(F13, 2) +
					pow(F21, 2) + pow(F22, 2) + pow(F23, 2) +
					pow(F31, 2) + pow(F32, 2) + pow(F33, 2);
				qeal J = F.determinant();
				ENuMaterial * material = downcastENuMaterial(m->getTetMeshHandle()->getElementMaterial(eid));

				qeal miu = material->getMu();
				qeal lamda = material->getLambda();
				qeal alpha = 1.0 + (3.0 * miu) / (4.0 * lamda);

				qeal v = m->getTetMeshHandle()->getTetElementParam(eid)->volume;
				qeal e = 0.5 * miu * (Ic - 3.0) + 0.5 * lamda * (J - alpha) *  (J - alpha) - 0.5 * miu * log(Ic + 1);
				potential += v * e;
			}
		}

		return potential;
	}

	void FemSimulator::computeElasticsEnergy(VectorX & X, VectorX & Xtilde, SparseMatrix & Mass, qeal & h, qeal & E)
	{
		E = 0.5 * (X).transpose() * (Mass * (X));
		E += computeNeoHookeanEnergy(X);
		E -= _sysExternalForce.dot(X);
	}

	void FemSimulator::stableNeoHookean(int eid, BaseTetMeshHandle* tet, Matrix3& F, Matrix3& FPK, MatrixX&PF)
	{
		// 3*3*3*3 tensor is represented as 9*9 matrix;
		qeal F11 = F.data()[0];  qeal F12 = F.data()[3];  qeal F13 = F.data()[6];
		qeal F21 = F.data()[1];  qeal F22 = F.data()[4];  qeal F23 = F.data()[7];
		qeal F31 = F.data()[2];  qeal F32 = F.data()[5];  qeal F33 = F.data()[8];
		qeal Ic = pow(F11, 2) + pow(F12, 2) + pow(F13, 2) +
			pow(F21, 2) + pow(F22, 2) + pow(F23, 2) +
			pow(F31, 2) + pow(F32, 2) + pow(F33, 2);

		qeal J = F.determinant();

		Matrix3 dJdF;
		qeal dJdF11 = F22 * F33 - F23 * F32;   qeal dJdF12 = F23 * F31 - F21 * F33;   qeal dJdF13 = F21 * F32 - F22 * F31;
		qeal dJdF21 = F13 * F32 - F12 * F33;   qeal dJdF22 = F11 * F33 - F13 * F31;   qeal dJdF23 = F12 * F31 - F11 * F32;
		qeal dJdF31 = F12 * F23 - F13 * F22;   qeal dJdF32 = F13 * F21 - F11 * F23;   qeal dJdF33 = F11 * F22 - F12 * F21;

		dJdF.data()[0] = dJdF11;   dJdF.data()[3] = dJdF12;   dJdF.data()[6] = dJdF13;
		dJdF.data()[1] = dJdF21;   dJdF.data()[4] = dJdF22;   dJdF.data()[7] = dJdF23;
		dJdF.data()[2] = dJdF31;   dJdF.data()[5] = dJdF32;   dJdF.data()[8] = dJdF33;

		ENuMaterial * material = downcastENuMaterial(tet->getElementMaterial(eid));

		qeal lameMiu = material->getMu();
		qeal lameLamda = material->getLambda();


		//qeal miu = (4.0 * lameMiu) / 3.0;
	//	qeal lamda = lameLamda + (5.0 * lameMiu) / 6.0;

		qeal miu = lameMiu;
		qeal lamda = lameLamda;

		qeal alpha = 1.0 + (3.0 * miu) / (4.0 * lamda);
		qeal a0 = miu * (1.0 - 1.0 / (Ic + 1.0));
		qeal a1 = 2.0 * miu / (pow(Ic + 1.0, 2));
		qeal a2 = lamda;
		qeal a3 = lamda * (J - alpha);

		FPK.setZero();
		FPK = a0 * F + a3 * dJdF;
		Matrix3 mU, mV, mSingularF;
		computeSVD(F, mU, mV, mSingularF, 1e-8, 1);
		qeal singularF0 = mSingularF.data()[0], singularF1 = mSingularF.data()[4], singularF2 = mSingularF.data()[8];

		MatrixX eigenValues(9, 9);
		eigenValues.setZero();
		MatrixX eigenVector(9, 9);
		eigenVector.setZero();
		//
		Matrix3 A1;
		A1.setIdentity();
		Matrix3 A2;
		A2(0, 0) = singularF0 * singularF0;
		A2(1, 1) = singularF1 * singularF1;
		A2(2, 2) = singularF2 * singularF2;
		A2(0, 1) = A2(1, 0) = singularF0 * singularF1;
		A2(0, 2) = A2(2, 0) = singularF0 * singularF2;
		A2(2, 1) = A2(1, 2) = singularF2 * singularF1;

		Matrix3 A3;
		A3(0, 0) = singularF1 * singularF1 *  singularF2 * singularF2;
		A3(1, 1) = singularF0 * singularF0 *  singularF2 * singularF2;
		A3(2, 2) = singularF1 * singularF1 *  singularF0 * singularF0;
		A3(0, 1) = A3(1, 0) = singularF0 * singularF1 * singularF2 * singularF2;
		A3(0, 2) = A3(2, 0) = singularF0 * singularF2 * singularF1 * singularF1;
		A3(2, 1) = A3(1, 2) = singularF2 * singularF1 * singularF0 * singularF0;

		Matrix3 A4;
		A4.setZero();
		A4(0, 1) = A4(1, 0) = singularF2;
		A4(0, 2) = A4(2, 0) = singularF1;
		A4(1, 2) = A4(2, 1) = singularF0;

		Matrix3 Ax = a0 * A1 + a1 * A2 + a2 * A3 + a3 * A4;
		const Eigen::SelfAdjointEigenSolver<Matrix3> Aeigs(Ax);
		Vector3 evs = Aeigs.eigenvalues();
		eigenValues(0, 0) = evs[0];
		eigenValues(1, 1) = evs[1];
		eigenValues(2, 2) = evs[2];
		Matrix3 eUs = Aeigs.eigenvectors();

		Matrix3 eUs0, eUs1, eUs2;
		eUs0.setZero();
		eUs1.setZero();
		eUs2.setZero();
		eUs0(0, 0) = eUs.data()[0];
		eUs0(1, 1) = eUs.data()[1];
		eUs0(2, 2) = eUs.data()[2];

		eUs1(0, 0) = eUs.data()[3];
		eUs1(1, 1) = eUs.data()[4];
		eUs1(2, 2) = eUs.data()[5];

		eUs2(0, 0) = eUs.data()[6];
		eUs2(1, 1) = eUs.data()[7];
		eUs2(2, 2) = eUs.data()[8];

		Matrix3 ev0 = mU * eUs0 * mV.transpose();
		Matrix3 ev1 = mU * eUs1 * mV.transpose();
		Matrix3 ev2 = mU * eUs2 * mV.transpose();

		for (int i = 0; i < 9; i++)
			eigenVector.data()[i] = ev0.data()[i];
		for (int i = 0; i < 9; i++)
			eigenVector.data()[9 + i] = ev1.data()[i];
		for (int i = 0; i < 9; i++)
			eigenVector.data()[18 + i] = ev2.data()[i];

		eigenValues(3, 3) = lamda * (J - alpha) * singularF0 + a0;
		eigenVector.data()[3 * 9 + 0] = -mU.data()[6] * mV.data()[3] + mU.data()[3] * mV.data()[6];
		eigenVector.data()[3 * 9 + 1] = -mU.data()[7] * mV.data()[3] + mU.data()[4] * mV.data()[6];
		eigenVector.data()[3 * 9 + 2] = -mU.data()[8] * mV.data()[3] + mU.data()[5] * mV.data()[6];

		eigenVector.data()[3 * 9 + 3] = -mU.data()[6] * mV.data()[4] + mU.data()[3] * mV.data()[7];
		eigenVector.data()[3 * 9 + 4] = -mU.data()[7] * mV.data()[4] + mU.data()[4] * mV.data()[7];
		eigenVector.data()[3 * 9 + 5] = -mU.data()[8] * mV.data()[4] + mU.data()[5] * mV.data()[7];

		eigenVector.data()[3 * 9 + 6] = -mU.data()[6] * mV.data()[5] + mU.data()[3] * mV.data()[8];
		eigenVector.data()[3 * 9 + 7] = -mU.data()[7] * mV.data()[5] + mU.data()[4] * mV.data()[8];
		eigenVector.data()[3 * 9 + 8] = -mU.data()[8] * mV.data()[5] + mU.data()[5] * mV.data()[8];
		eigenVector.col(3) *= 1.0 / sqrt(2.0);

		eigenValues(4, 4) = lamda * (J - alpha) * singularF1 + a0;
		eigenVector.data()[4 * 9 + 0] = -mU.data()[6] * mV.data()[0] + mU.data()[0] * mV.data()[6];
		eigenVector.data()[4 * 9 + 1] = -mU.data()[7] * mV.data()[0] + mU.data()[1] * mV.data()[6];
		eigenVector.data()[4 * 9 + 2] = -mU.data()[8] * mV.data()[0] + mU.data()[2] * mV.data()[6];

		eigenVector.data()[4 * 9 + 3] = -mU.data()[6] * mV.data()[1] + mU.data()[0] * mV.data()[7];
		eigenVector.data()[4 * 9 + 4] = -mU.data()[7] * mV.data()[1] + mU.data()[1] * mV.data()[7];
		eigenVector.data()[4 * 9 + 5] = -mU.data()[8] * mV.data()[1] + mU.data()[2] * mV.data()[7];

		eigenVector.data()[4 * 9 + 6] = -mU.data()[6] * mV.data()[2] + mU.data()[0] * mV.data()[8];
		eigenVector.data()[4 * 9 + 7] = -mU.data()[7] * mV.data()[2] + mU.data()[1] * mV.data()[8];
		eigenVector.data()[4 * 9 + 8] = -mU.data()[8] * mV.data()[2] + mU.data()[2] * mV.data()[8];
		eigenVector.col(4) *= 1.0 / sqrt(2.0);
		//

		eigenValues(5, 5) = lamda * (J - alpha) * singularF2 + a0;
		eigenVector.data()[5 * 9 + 0] = -mU.data()[3] * mV.data()[0] + mU.data()[0] * mV.data()[3];
		eigenVector.data()[5 * 9 + 1] = -mU.data()[4] * mV.data()[0] + mU.data()[1] * mV.data()[3];
		eigenVector.data()[5 * 9 + 2] = -mU.data()[5] * mV.data()[0] + mU.data()[2] * mV.data()[3];

		eigenVector.data()[5 * 9 + 3] = -mU.data()[3] * mV.data()[1] + mU.data()[0] * mV.data()[4];
		eigenVector.data()[5 * 9 + 4] = -mU.data()[4] * mV.data()[1] + mU.data()[1] * mV.data()[4];
		eigenVector.data()[5 * 9 + 5] = -mU.data()[5] * mV.data()[1] + mU.data()[2] * mV.data()[4];

		eigenVector.data()[5 * 9 + 6] = -mU.data()[3] * mV.data()[2] + mU.data()[0] * mV.data()[5];
		eigenVector.data()[5 * 9 + 7] = -mU.data()[4] * mV.data()[2] + mU.data()[1] * mV.data()[5];
		eigenVector.data()[5 * 9 + 8] = -mU.data()[5] * mV.data()[2] + mU.data()[2] * mV.data()[5];
		eigenVector.col(5) *= 1.0 / sqrt(2.0);
		//

		eigenValues(6, 6) = -lamda * (J - alpha) * singularF0 + a0;
		eigenVector.data()[6 * 9 + 0] = mU.data()[6] * mV.data()[3] + mU.data()[3] * mV.data()[6];
		eigenVector.data()[6 * 9 + 1] = mU.data()[7] * mV.data()[3] + mU.data()[4] * mV.data()[6];
		eigenVector.data()[6 * 9 + 2] = mU.data()[8] * mV.data()[3] + mU.data()[5] * mV.data()[6];

		eigenVector.data()[6 * 9 + 3] = mU.data()[6] * mV.data()[4] + mU.data()[3] * mV.data()[7];
		eigenVector.data()[6 * 9 + 4] = mU.data()[7] * mV.data()[4] + mU.data()[4] * mV.data()[7];
		eigenVector.data()[6 * 9 + 5] = mU.data()[8] * mV.data()[4] + mU.data()[5] * mV.data()[7];

		eigenVector.data()[6 * 9 + 6] = mU.data()[6] * mV.data()[5] + mU.data()[3] * mV.data()[8];
		eigenVector.data()[6 * 9 + 7] = mU.data()[7] * mV.data()[5] + mU.data()[4] * mV.data()[8];
		eigenVector.data()[6 * 9 + 8] = mU.data()[8] * mV.data()[5] + mU.data()[5] * mV.data()[8];
		eigenVector.col(6) *= 1.0 / sqrt(2.0);
		//
		eigenValues(7, 7) = -lamda * (J - alpha) * singularF1 + a0;
		eigenVector.data()[7 * 9 + 0] = mU.data()[6] * mV.data()[0] + mU.data()[0] * mV.data()[6];
		eigenVector.data()[7 * 9 + 1] = mU.data()[7] * mV.data()[0] + mU.data()[1] * mV.data()[6];
		eigenVector.data()[7 * 9 + 2] = mU.data()[8] * mV.data()[0] + mU.data()[2] * mV.data()[6];

		eigenVector.data()[7 * 9 + 3] = mU.data()[6] * mV.data()[1] + mU.data()[0] * mV.data()[7];
		eigenVector.data()[7 * 9 + 4] = mU.data()[7] * mV.data()[1] + mU.data()[1] * mV.data()[7];
		eigenVector.data()[7 * 9 + 5] = mU.data()[8] * mV.data()[1] + mU.data()[2] * mV.data()[7];

		eigenVector.data()[7 * 9 + 6] = mU.data()[6] * mV.data()[2] + mU.data()[0] * mV.data()[8];
		eigenVector.data()[7 * 9 + 7] = mU.data()[7] * mV.data()[2] + mU.data()[1] * mV.data()[8];
		eigenVector.data()[7 * 9 + 8] = mU.data()[8] * mV.data()[2] + mU.data()[2] * mV.data()[8];
		eigenVector.col(7) *= 1.0 / sqrt(2.0);

		//
		eigenValues(8, 8) = -lamda * (J - alpha) * singularF2 + a0;
		eigenVector.data()[8 * 9 + 0] = mU.data()[3] * mV.data()[0] + mU.data()[0] * mV.data()[3];
		eigenVector.data()[8 * 9 + 1] = mU.data()[4] * mV.data()[0] + mU.data()[1] * mV.data()[3];
		eigenVector.data()[8 * 9 + 2] = mU.data()[5] * mV.data()[0] + mU.data()[2] * mV.data()[3];

		eigenVector.data()[8 * 9 + 3] = mU.data()[3] * mV.data()[1] + mU.data()[0] * mV.data()[4];
		eigenVector.data()[8 * 9 + 4] = mU.data()[4] * mV.data()[1] + mU.data()[1] * mV.data()[4];
		eigenVector.data()[8 * 9 + 5] = mU.data()[5] * mV.data()[1] + mU.data()[2] * mV.data()[4];

		eigenVector.data()[8 * 9 + 6] = mU.data()[3] * mV.data()[2] + mU.data()[0] * mV.data()[5];
		eigenVector.data()[8 * 9 + 7] = mU.data()[4] * mV.data()[2] + mU.data()[1] * mV.data()[5];
		eigenVector.data()[8 * 9 + 8] = mU.data()[5] * mV.data()[2] + mU.data()[2] * mV.data()[5];
		eigenVector.col(8) *= 1.0 / sqrt(2.0);

		for (int i = 0; i < 9; i++)
			if (eigenValues(i, i) < 1e-13) eigenValues(i, i) = 0.0;
		PF = eigenVector * eigenValues * eigenVector.transpose();		
	}		

	void FemSimulator::commonNeoHookean(int eid, BaseTetMeshHandle * tet, Matrix3 & F, Matrix3 & FPK, MatrixX & PF)
	{
		AutoFlipSVD<Matrix3> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Vector3 S = svd.singularValues();
		if (S[2] < 1e-14)
		{
			S[2] *= -1;
			Matrix3 Si;
			Si.setZero();
			Si.data()[0] = S[0];
			Si.data()[4] = S[1];
			Si.data()[8] = S[2];
			Matrix3 U = svd.matrixU();
			Matrix3 V = svd.matrixV();
			F = U * Si * V.transpose();
			svd.compute(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		}

		ENuMaterial * material = downcastENuMaterial(tet->getElementMaterial(eid));
		qeal lameMiu = material->getMu();
		qeal lameLamda = material->getLambda();
		compute_dE_div_dF(F, svd, lameMiu, lameLamda, FPK);
		compute_dP_div_dF(svd, lameMiu, lameLamda, PF);	
	}

	void FemSimulator::computeCommonNeoHookeanPartialPF(VectorX & displacement)
	{
		for (size_t i = 0; i < models.size(); i++)
		{
			BaseModel* m = getModel(i);
			for (size_t eid = 0; eid < m->tetElementNum; eid++)
			{
				int geid = m->getTetElementOverallId(eid);
				Vector4i indices = m->getTetElementOverall(eid);
				BaseTetElementParam* para = m->getTetMeshHandle()->getTetElementParam(eid);
				Matrix3 Ds;
				para->computeElementDeformationByShapeMatrix(Ds, displacement.data(), indices.data());

				Matrix3 F = Ds * para->invDm;

				Matrix3 FPK;
				commonNeoHookean(eid, m->getTetMeshHandle(), F, FPK, dPdF[geid]);
				elementsInternalForce[geid] = (-1.0 * para->volume) * FPK * para->invDmT;
			}
		}

	}

	void FemSimulator::doTime(int frame)
	{
		int numIter = 0;
		qeal TOL = solverConfig->epsilon * solverConfig->epsilon;
		qeal errorAtFirstIter = 0;
		qeal errorQuotient;
		qeal preerror = QEAL_MAX;
		_sysXn = _sysX;
		_sysVn = _sysVelocity;
		_sysAccn = _sysAccrelation;

		//predictive vel and acc at t+1 time
		_sysAccrelation = - solverConfig->alpha[1] * _sysVn - solverConfig->alpha[2] * _sysAccn;
		_sysVelocity = solverConfig->alpha[4] * _sysVn + solverConfig->alpha[5] * _sysAccn;

		computeExternalForce();
		VectorX groundCollideForce(_sysDim);
		groundCollideForce.setZero();
		for (int i = 0; i < totalTetPointsNum; i++)
		{
			qeal x, y, z;
			getTetPoint(i, x, y, z);
			if (y < 0)
			{
				qeal force = -y * 100.0;
				groundCollideForce[3 * i + 1] += force;
				groundCollideForce[3 * i + 1] += 2.0 * _sysVelocity[3 * i + 1];
			}
		}

		do 
		{
			computePartialPF(_sysX);
			computeInternalForce();
			computeStiffnessMatrix();

			_sysDampingMatrix = solverConfig->dampingMassCoef * _sysMassMatrix + solverConfig->dampingStiffnessCoef * _sysStiffnessMatrix;

			_sysMatrix = solverConfig->alpha[0] * _sysMassMatrix + solverConfig->alpha[3] * _sysDampingMatrix + _sysStiffnessMatrix;

			_sysRhs = _sysMassMatrix * _sysAccrelation + _sysDampingMatrix * _sysVelocity + _sysInternalForce - _sysExternalForce - groundCollideForce;
			_sysRhs *= -1;

			qeal error = _sysRhs.norm() * _sysRhs.norm();

			if (numIter == 0)
			{
				errorAtFirstIter = error;
				errorQuotient = 1;
				//computeNewtonSolverEnergy(Eprev, E, _sysMatrix, _sysRhs, _sysCurrentPosition);
			}
			else
			{
				errorQuotient = error / errorAtFirstIter;
			}
			//cout << "Iter " << numIter << " :  " << error << endl;
			if (std::abs(error - preerror) / std::abs(error) < solverConfig->epsilon || abs(error) < TOL)
			{
				break;
			}

			preerror = error;

			if (errorQuotient < TOL)
			{
				break;
			}

			sparseMatrixRemoveRows(_sysMatrix, _sysSubMatrix);
			vectorRemoveRows(_sysRhs, _sysSubRhs);
			_sparseLDLT.compute(_sysSubMatrix);

			_sysSubXtilde = _sparseLDLT.solve(_sysSubRhs);
			vectorInsertRows(_sysSubXtilde, _sysXtilde);

			_sysX += _sysXtilde;
			
			_sysAccrelation = solverConfig->alpha[0] * (_sysX - _sysXn) - solverConfig->alpha[1] * _sysVn - solverConfig->alpha[2] * _sysAccn;
			_sysVelocity = solverConfig->alpha[3] * (_sysX - _sysXn) + solverConfig->alpha[4] * _sysVn + solverConfig->alpha[5] * _sysAccn;
			numIter++;

		} while (numIter < solverConfig->maxIterations);
	}



	void FemSimulator::doTime2(int frame)
	{
		int newton_iter = 0;
		qeal tol = 1e-2 * (_diagLen / 2);
		qeal TOLs = solverConfig->epsilon * solverConfig->epsilon;
		qeal errorAtFirstIter = 0;
		qeal errorQuotient;
		qeal preerror = QEAL_MAX;

		qeal dt = _timeStep;
		VectorX xn = _sysX;
		// update_predictive_pos
		VectorX x_tilde(xn);
		x_tilde = xn + _sysVelocity * dt;
		computeExternalForce();
		for (; newton_iter < 200; ++newton_iter)
		{
			//¼ÆËãgradient ºÍ hessain
			computePartialPF(xn);
			computeInternalForce();
			computeStiffnessMatrix();
			_sysMatrix = _sysMassMatrix + dt * dt * _sysStiffnessMatrix;
			_sysRhs = _sysMassMatrix * (xn - x_tilde) + dt * dt * (_sysInternalForce - _sysExternalForce);
			_sysRhs *= -1;

			//solver
			Eigen::SimplicialLLT<SparseMatrix> solver;
			solver.compute(_sysMatrix);
			VectorX direction = solver.solve(_sysRhs);
			//
			qeal res = direction.cwiseAbs().maxCoeff() / dt;
			if (newton_iter > 0 && res <= tol)
			{
				std::cout << "newton_iter " << newton_iter << " " << res << " " << tol << res - tol << std::endl;
				break;
			}
			//
			//
			qeal E0, E;
			E0 = getEnergy(xn, x_tilde);
			qeal alpha = 1.0;
			VectorX new_x = xn + alpha * direction;
			E = getEnergy(new_x, x_tilde);

			while (E > E0)
			{
				std::cout << frame <<" " << newton_iter << " E0 " << E0 << " E " << E <<" " << E - E0 << " " << alpha<< std::endl;
				alpha *= 0.5;
				new_x = xn + alpha * direction;				
				E = getEnergy(new_x, x_tilde);
			}
			E0 = E;
			xn = new_x;
			_sysX = xn;
		}
		
		VectorX new_a = (1.0 / (dt * dt))* (xn - x_tilde);
		_sysVelocity += dt * new_a;
	}

	void FemSimulator::doComomNHTime(int frame)
	{
		std::cout << "----- doComomNHTime ----" << std::endl;
		int numIter = 0;
		qeal TOL = solverConfig->epsilon * solverConfig->epsilon;
		qeal errorAtFirstIter = 0;
		qeal errorQuotient;
		qeal preerror = QEAL_MAX;
		_sysXn = _sysX;
		_sysVn = _sysVelocity;

		computeExternalForce();
		VectorX groundCollideForce(_sysDim);
		groundCollideForce.setZero();
		for (int i = 0; i < totalTetPointsNum; i++)
		{
			qeal x, y, z;
			getTetPoint(i, x, y, z);
			if (y < 0)
			{
				qeal force = -y * 100.0;
				groundCollideForce[3 * i + 1] += force;
				groundCollideForce[3 * i + 1] -= 2.0 * _sysVelocity[3 * i + 1];
			}
		}

		//_sysExternalForce += groundCollideForce;

		VectorX Xtilde = _sysXn + _timeStep * _sysVelocity + _timeStep * _timeStep * (_sysInverseMassMatrix * _sysExternalForce);


		do
		{
			//computeCommonNeoHookeanPartialPF(_sysX);
			computePartialPF(_sysX);
			computeInternalForce();
			computeStiffnessMatrix();

			_sysMatrix = _sysMassMatrix + _timeStep * _timeStep * _sysStiffnessMatrix;

			_sysRhs = _sysMassMatrix * (_sysX - _sysXn - _timeStep * _sysVn) + (_timeStep * _timeStep * (_sysInternalForce - _sysExternalForce));

			_sysRhs *= -1.0;

			qeal error = _sysRhs.norm() * _sysRhs.norm();

			//if (numIter == 0)
			//{
			//	errorAtFirstIter = error;
			//	errorQuotient = 1;
			//	//computeNewtonSolverEnergy(Eprev, E, _sysMatrix, _sysRhs, _sysCurrentPosition);
			//}
			//else
			//{
			//	errorQuotient = error / errorAtFirstIter;
			//}
			std::cout << "Iter " << numIter << " :  " << error << std::endl;
			if (abs(error) < TOL)
			{
				break;
			}

			//preerror = error;

			//if (errorQuotient < TOL)
			//{
			//	break;
			//}

			sparseMatrixRemoveRows(_sysMatrix, _sysSubMatrix);
			vectorRemoveRows(_sysRhs, _sysSubRhs);
			_sparseLDLT.compute(_sysSubMatrix);

			_sysSubXtilde = _sparseLDLT.solve(_sysSubRhs);
			vectorInsertRows(_sysSubXtilde, _sysXtilde);

			_sysX += _sysXtilde;

			//for (int i = 0; i < _sysX.size(); i++)
			//{
			//	std::cout << i << " " << _sysXtilde[i] << " " << _sysRhs[i] << std::endl;
			//}


			numIter++;

			qeal infNorm = 0.0;
			for (int i = 0; i < _sysXtilde.size(); ++i) {
				if (infNorm < std::abs(_sysXtilde[i])) {
					infNorm = std::abs(_sysXtilde[i]);
				}
			}

			std::cout << "inNorm: " << (infNorm / _timeStep) << " " << (1e-2 * _diagLen) << std::endl;
			if ((infNorm / _timeStep) < (1e-2 * _diagLen))
				break;


		} while (true);

		_sysVelocity = (1.0 / (_timeStep)) * (_sysX - _sysXn);
	}

    void FemSimulator::NewtonLineSearch(qeal & E, qeal & Eprev, qeal& alpha, SparseMatrix & A, VectorX & rhs, VectorX & x)
	{
		alpha = 1.0;
		do
		{
			x = alpha * x;
			computeNewtonSolverEnergy(E, Eprev, A, rhs, x);
			alpha *= 0.5;
		} while (E > Eprev);
		Eprev = E;
	}

	void FemSimulator::saveFrameStatus(int frame, std::ofstream& fout)
	{
		EigenMatrixIO::write_binary(fout, _sysSurfaceOriginalPosition);
		EigenMatrixIO::write_binary(fout, _sysOriginalPosition);
		EigenMatrixIO::write_binary(fout, _sysCurrentPosition);
		EigenMatrixIO::write_binary(fout, _sysVelocity);
		EigenMatrixIO::write_binary(fout, _sysAccrelation);
		EigenMatrixIO::write_binary(fout, _sysExternalForce);
		EigenMatrixIO::write_binary(fout, _sysInternalForce);
		EigenMatrixIO::write_binary(fout, _sysDampingForce);
		EigenMatrixIO::write_binary(fout, _sysGravityForce);
		EigenMatrixIO::write_binary(fout, _sysMouseForce);
		EigenMatrixIO::write_binary(fout, _sysCollideForce);
		EigenMatrixIO::write_binary(fout, _meshCollisionForce);
		EigenMatrixIO::write_binary(fout, _sysX);
		EigenMatrixIO::write_binary(fout, _sysXtilde);
	}

	void FemSimulator::recoverFrameStatus(int frame, std::ifstream& fin)
	{
		EigenMatrixIO::read_binary(fin, _sysSurfaceOriginalPosition);
		EigenMatrixIO::read_binary(fin, _sysOriginalPosition);
		EigenMatrixIO::read_binary(fin, _sysCurrentPosition);
		EigenMatrixIO::read_binary(fin, _sysVelocity);
		EigenMatrixIO::read_binary(fin, _sysAccrelation);
		EigenMatrixIO::read_binary(fin, _sysExternalForce);
		EigenMatrixIO::read_binary(fin, _sysInternalForce);
		EigenMatrixIO::read_binary(fin, _sysDampingForce);
		EigenMatrixIO::read_binary(fin, _sysGravityForce);
		EigenMatrixIO::read_binary(fin, _sysMouseForce);
		EigenMatrixIO::read_binary(fin, _sysCollideForce);
		EigenMatrixIO::read_binary(fin, _meshCollisionForce);
		EigenMatrixIO::read_binary(fin, _sysX);
		EigenMatrixIO::read_binary(fin, _sysXtilde);
	}

	void FemSimulator::computeCommonNHElasticsEnergy(VectorX & X, VectorX & Xtilde, SparseMatrix & Mass, qeal & h, qeal & E)
	{
		E = 0.5 * (X - Xtilde).transpose() * (Mass * (X - Xtilde));
		E += h * h * computeCommonNeoHookeanEnergy(X);
	}

	qeal FemSimulator::computeCommonNeoHookeanEnergy(VectorX & displacement)
	{
		qeal potential = 0.0;
		for (size_t i = 0; i < models.size(); i++)
		{
			BaseModel* m = getModel(i);
			for (size_t eid = 0; eid < m->tetElementNum; eid++)
			{
				int geid = m->getTetElementOverallId(eid);
				Vector4i indices = m->getTetElementOverall(eid);
				BaseTetElementParam* para = m->getTetMeshHandle()->getTetElementParam(eid);
				Matrix3 Ds;
				para->computeElementDeformationByShapeMatrix(Ds, displacement.data(), indices.data());
				Matrix3 F = Ds * para->invDm;

				potential += computeCommonNeoHookeanEnergy(eid, m->getTetMeshHandle(), F);
			}
		}
		return potential;
	}

	qeal FemSimulator::computeCommonNeoHookeanEnergy(int eid, BaseTetMeshHandle * tet, Matrix3 & F)
	{
		ENuMaterial * material = downcastENuMaterial(tet->getElementMaterial(eid));
		qeal lameMiu = material->getMu();
		qeal lameLamda = material->getLambda();
		AutoFlipSVD<Matrix3> svd(F);

		qeal E = 0;
		Vector3 S = svd.singularValues();
		computeElementCNHEnergy(S, lameMiu, lameLamda, E);

		return E;
	}

	void FemSimulator::compute_dP_div_dF(const AutoFlipSVD<Matrix3>& svd, qeal miu, qeal lambda, MatrixX& dP_div_dF, bool projectSPD)
	{
		int dim = 3;
		// compute A
		Vector3 sigma = svd.singularValues();
		Vector3 dE_div_dsigma;
		compute_dE_div_dsigma(sigma, miu, lambda, dE_div_dsigma);
		Matrix3 d2E_div_dsigma2;
		compute_d2E_div_dsigma2(sigma, miu, lambda, d2E_div_dsigma2);
	
		makePD3d(d2E_div_dsigma2); //TODO: use implicit QR to accelerate

		// compute B
		const int Cdim2 = dim * (dim - 1) / 2;
		VectorX BLeftCoef(Cdim2);
		compute_BLeftCoef(sigma, miu, lambda, BLeftCoef);

		std::vector<Matrix2> B(Cdim2);
		for (int cI = 0; cI < Cdim2; cI++) {
			int cI_post = (cI + 1) % dim;

			qeal rightCoef = dE_div_dsigma[cI] + dE_div_dsigma[cI_post];
			qeal sum_sigma = sigma[cI] + sigma[cI_post];
			const qeal eps = 1.0e-6;
			if (sum_sigma < eps) {
				rightCoef /= 2.0 * eps;
			}
			else {
				rightCoef /= 2.0 * sum_sigma;
			}

			const qeal& leftCoef = BLeftCoef[cI];
			B[cI](0, 0) = B[cI](1, 1) = leftCoef + rightCoef;
			B[cI](0, 1) = B[cI](1, 0) = leftCoef - rightCoef;
			//if (projectSPD) {
				makePD2d(B[cI]);
			//}
		}

		// compute M using A(d2E_div_dsigma2) and B
		MatrixX M(9, 9);
		M.setZero();
		// A
		M(0, 0) = d2E_div_dsigma2(0, 0);
		M(0, 4) = d2E_div_dsigma2(0, 1);
		M(0, 8) = d2E_div_dsigma2(0, 2);
		M(4, 0) = d2E_div_dsigma2(1, 0);
		M(4, 4) = d2E_div_dsigma2(1, 1);
		M(4, 8) = d2E_div_dsigma2(1, 2);
		M(8, 0) = d2E_div_dsigma2(2, 0);
		M(8, 4) = d2E_div_dsigma2(2, 1);
		M(8, 8) = d2E_div_dsigma2(2, 2);
		// B01
		M(1, 1) = B[0](0, 0);
		M(1, 3) = B[0](0, 1);
		M(3, 1) = B[0](1, 0);
		M(3, 3) = B[0](1, 1);
		// B12
		M(5, 5) = B[1](0, 0);
		M(5, 7) = B[1](0, 1);
		M(7, 5) = B[1](1, 0);
		M(7, 7) = B[1](1, 1);
		// B20
		M(2, 2) = B[2](1, 1);
		M(2, 6) = B[2](1, 0);
		M(6, 2) = B[2](0, 1);
		M(6, 6) = B[2](0, 0);

		// compute dP_div_dF
		dP_div_dF.resize(9, 9);
		Matrix3 U, V;
		U = svd.matrixU();
		V = svd.matrixV();

		for (int i = 0; i < dim; i++) {
			int _dim_i = i * dim;
			for (int j = 0; j < dim; j++) {
				int ij = _dim_i + j;
				for (int r = 0; r < dim; r++) {
					int _dim_r = r * dim;
					for (int s = 0; s < dim; s++) {
						int rs = _dim_r + s;
						if (ij > rs) {
							// bottom left, same as upper right
							continue;
						}

						dP_div_dF(ij, rs) = 
							M(0, 0) * U(i, 0) * V(j, 0) * U(r, 0) * V(s, 0) 
							+ M(0, 4) * U(i, 0) * V(j, 0) * U(r, 1) * V(s, 1) 
							+ M(0, 8) * U(i, 0) * V(j, 0) * U(r, 2) * V(s, 2) 
							+ M(4, 0) * U(i, 1) * V(j, 1) * U(r, 0) * V(s, 0) 
							+ M(4, 4) * U(i, 1) * V(j, 1) * U(r, 1) * V(s, 1) 
							+ M(4, 8) * U(i, 1) * V(j, 1) * U(r, 2) * V(s, 2) 
							+ M(8, 0) * U(i, 2) * V(j, 2) * U(r, 0) * V(s, 0) 
							+ M(8, 4) * U(i, 2) * V(j, 2) * U(r, 1) * V(s, 1) 
							+ M(8, 8) * U(i, 2) * V(j, 2) * U(r, 2) * V(s, 2) 
							+ M(1, 1) * U(i, 0) * V(j, 1) * U(r, 0) * V(s, 1) 
							+ M(1, 3) * U(i, 0) * V(j, 1) * U(r, 1) * V(s, 0) 
							+ M(3, 1) * U(i, 1) * V(j, 0) * U(r, 0) * V(s, 1) 
							+ M(3, 3) * U(i, 1) * V(j, 0) * U(r, 1) * V(s, 0) 
							+ M(5, 5) * U(i, 1) * V(j, 2) * U(r, 1) * V(s, 2) 
							+ M(5, 7) * U(i, 1) * V(j, 2) * U(r, 2) * V(s, 1) 
							+ M(7, 5) * U(i, 2) * V(j, 1) * U(r, 1) * V(s, 2) 
							+ M(7, 7) * U(i, 2) * V(j, 1) * U(r, 2) * V(s, 1) 
							+ M(2, 2) * U(i, 0) * V(j, 2) * U(r, 0) * V(s, 2) 
							+ M(2, 6) * U(i, 0) * V(j, 2) * U(r, 2) * V(s, 0) 
							+ M(6, 2) * U(i, 2) * V(j, 0) * U(r, 0) * V(s, 2) 
							+ M(6, 6) * U(i, 2) * V(j, 0) * U(r, 2) * V(s, 0);

						if (ij < rs) {
							dP_div_dF(rs, ij) = dP_div_dF(ij, rs);
						}
					}
				}
			}
		}

	}

	void FemSimulator::compute_dE_div_dsigma(Vector3 & singularValues, qeal miu, qeal lambda, Vector3 & dE_div_dsigma)
	{
		if (miu == 0.0 && lambda == 0.0) {
			dE_div_dsigma.setZero();
			return;
		}

		const qeal log_sigmaProd = std::log(singularValues.prod());
		const qeal inv0 = 1.0 / singularValues[0];
		dE_div_dsigma[0] = miu * (singularValues[0] - inv0) + lambda * inv0 * log_sigmaProd;
		const qeal inv1 = 1.0 / singularValues[1];
		dE_div_dsigma[1] = miu * (singularValues[1] - inv1) + lambda * inv1 * log_sigmaProd;
		const qeal inv2 = 1.0 / singularValues[2];
		dE_div_dsigma[2] = miu * (singularValues[2] - inv2) + lambda * inv2 * log_sigmaProd;
	}

	void FemSimulator::compute_d2E_div_dsigma2(Vector3 & singularValues, qeal miu, qeal lambda, Matrix3 & d2E_div_dsigma2)
	{
		if (miu == 0.0 && lambda == 0.0) {
			d2E_div_dsigma2.setZero();
			return;
		}

		const qeal log_sigmaProd = std::log(singularValues.prod());

		const qeal inv2_0 = 1.0 / singularValues[0] / singularValues[0];
		d2E_div_dsigma2(0, 0) = miu * (1.0 + inv2_0) - lambda * inv2_0 * (log_sigmaProd - 1.0);
		const qeal inv2_1 = 1.0 / singularValues[1] / singularValues[1];
		d2E_div_dsigma2(1, 1) = miu * (1.0 + inv2_1) - lambda * inv2_1 * (log_sigmaProd - 1.0);
		d2E_div_dsigma2(0, 1) = d2E_div_dsigma2(1, 0) = lambda / singularValues[0] / singularValues[1];
		const qeal inv2_2 = 1.0 / singularValues[2] / singularValues[2];
		d2E_div_dsigma2(2, 2) = miu * (1.0 + inv2_2) - lambda * inv2_2 * (log_sigmaProd - 1.0);
		d2E_div_dsigma2(1, 2) = d2E_div_dsigma2(2, 1) = lambda / singularValues[1] / singularValues[2];
		d2E_div_dsigma2(2, 0) = d2E_div_dsigma2(0, 2) = lambda / singularValues[2] / singularValues[0];
	}

	void FemSimulator::compute_BLeftCoef(Vector3 & singularValues, qeal miu, qeal lambda, VectorX & BLeftCoef)
	{
		if (miu == 0.0 && lambda == 0.0) {
			BLeftCoef.setZero();
			return;
		}

		//TODO: right coef also has analytical form
		const double sigmaProd = singularValues.prod();
		const double middle = miu - lambda * std::log(sigmaProd);
		BLeftCoef[0] = (miu + middle / singularValues[0] / singularValues[1]) / 2.0;
		BLeftCoef[1] = (miu + middle / singularValues[1] / singularValues[2]) / 2.0;
		BLeftCoef[2] = (miu + middle / singularValues[2] / singularValues[0]) / 2.0;
	}

	void FemSimulator::compute_dE_div_dF(Matrix3 & F, const AutoFlipSVD<Matrix3>& svd, qeal miu, qeal lambda, Matrix3 & dE_div_dF)
	{
		if (miu == 0.0 && lambda == 0.0) {
			dE_div_dF.setZero();
			return;
		}

		const qeal J = svd.singularValues().prod();
		Matrix3 FInvT;
		computeCofactorMtr(F, FInvT);
		FInvT /= J;
		dE_div_dF = miu * (F - FInvT) + lambda * std::log(J) * FInvT;
	}

	void FemSimulator::computeElementCNHEnergy(Vector3 & singularValues, qeal miu, qeal lambda, qeal & E)
	{
		if (miu == 0.0 && lambda == 0.0) {
			E = 0.0;
			return;
		}

		const qeal sigma2Sum = singularValues.squaredNorm();
		const qeal sigmaProd = singularValues.prod();
		const qeal log_sigmaProd = std::log(sigmaProd);

		E = miu / 2.0 * (sigma2Sum - 3) - (miu - lambda / 2.0 * log_sigmaProd) * log_sigmaProd;

	}

	void FemSimulator::computeCofactorMtr(const Matrix3 & F, Matrix3 & A)
	{
		A(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
		A(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
		A(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);
		A(1, 0) = F(0, 2) * F(2, 1) - F(0, 1) * F(2, 2);
		A(1, 1) = F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0);
		A(1, 2) = F(0, 1) * F(2, 0) - F(0, 0) * F(2, 1);
		A(2, 0) = F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1);
		A(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
		A(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
	}

	void FemSimulator::makePD3d(Matrix3 & symMtr)
	{
		Eigen::SelfAdjointEigenSolver<Matrix3> eigenSolver(symMtr);
		if (eigenSolver.eigenvalues()[0] >= 0.0) {
			return;
		}
		Eigen::DiagonalMatrix<qeal, 3> D(eigenSolver.eigenvalues());
		int rows = ((3 == Eigen::Dynamic) ? symMtr.rows() : 3);
		for (int i = 0; i < rows; i++) {
			if (D.diagonal()[i] < 0.0) {
				D.diagonal()[i] = 0.0;
			}
			else {
				break;
			}
		}
		symMtr = eigenSolver.eigenvectors() * D * eigenSolver.eigenvectors().transpose();
	}

	void FemSimulator::makePD2d(Matrix2 & symMtr)
	{
		// based on http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/

		const qeal a = symMtr(0, 0);
		const qeal b = (symMtr(0, 1) + symMtr(1, 0)) / 2.0;
		const qeal d = symMtr(1, 1);

		qeal b2 = b * b;
		qeal D = a * d - b2;
		qeal T_div_2 = (a + d) / 2.0;
		qeal sqrtTT4D = std::sqrt(T_div_2 * T_div_2 - D);
		qeal L2 = T_div_2 - sqrtTT4D;
		if (L2 < 0.0) {
			const qeal L1 = T_div_2 + sqrtTT4D;
			if (L1 <= 0.0) {
				symMtr.setZero();
			}
			else {
				if (b2 == 0.0) {
					symMtr << L1, 0.0, 0.0, 0.0;
				}
				else {
					const qeal L1md = L1 - d;
					const qeal L1md_div_L1 = L1md / L1;
					symMtr(0, 0) = L1md_div_L1 * L1md;
					symMtr(0, 1) = symMtr(1, 0) = b * L1md_div_L1;
					symMtr(1, 1) = b2 / L1;
				}
			}
		}
	}

	qeal FemSimulator::getEnergy(VectorX & xn, VectorX & x_tilde)
	{
		qeal E = 0;
		qeal dt = _timeStep;
		// compute static force energy
		qeal e0 = dt * dt * _sysExternalForce.dot(xn);
		e0 *= -1.0;
		//compute inertial energy
		qeal e1 = (xn - x_tilde).dot(_sysMassMatrix *  (xn - x_tilde));
		e1 *= 0.5;
		//compute elastics energy
		qeal e2 = computeNeoHookeanEnergy(xn);
		e2 *= dt * dt;
		E = e0 + e1 + e2;

		std::cout << e0 << " " << e1 << std::endl;
		system("pause");
		return E;
	}

	void FemSimulator::createSparseMapbyTopology(SparseMatrix& target, SparseMatrix& result)
	{
		int sysSubDim = 3 * (totalTetPointsNum - _constrainedTetPointsNum);

		result.resize(sysSubDim, sysSubDim);

		std::vector<TripletX> subCoef;

		for (int j = 0; j < target.outerSize(); ++j)
			for (SparseMatrix::InnerIterator it(target, j); it; ++it)
			{
				int i_p = _mapOldNew[it.row()];
				int j_p = _mapOldNew[it.col()];
				if (!(i_p == -1 || j_p == -1))
				{
					subCoef.push_back(TripletX(i_p, j_p, it.value()));
				}
			}

		result.setFromTriplets(subCoef.begin(), subCoef.end());

		int supersize = target.nonZeros();
		int subsize = result.nonZeros();

		_mapSparseMatrixEntryOldNew.resize(supersize);
		std::fill(_mapSparseMatrixEntryOldNew.begin(), _mapSparseMatrixEntryOldNew.end(), -1);

		SparseMatrixTopologyTYPE<qeal> supermatrix(&target);
		SparseMatrixTopologyTYPE<qeal> submatrix(&result);

		for (int j = 0; j < target.outerSize(); ++j)
			for (SparseMatrix::InnerIterator it(target, j); it; ++it)
			{
				int i_p = _mapOldNew[it.row()];
				int j_p = _mapOldNew[it.col()];
				if (!(i_p == -1 || j_p == -1))
				{
					_mapSparseMatrixEntryOldNew[supermatrix.getValueIndex(it.row(), it.col())] =
						submatrix.getValueIndex(i_p, j_p);
				}
			}
	}

	void FemSimulator::sparseMatrixRemoveRows(SparseMatrix& target, SparseMatrix& result)
	{
		int i = 0;
		for (int j = 0; j < target.outerSize(); ++j)
			for (SparseMatrix::InnerIterator it(target, j); it; ++it, ++i)
			{
				if (_mapSparseMatrixEntryOldNew[i] != -1)
				{
					result.valuePtr()[_mapSparseMatrixEntryOldNew[i]] = it.value();
				}
			}
	}

	void FemSimulator::vectorRemoveRows(VectorX & target, VectorX & result)
	{
		for (int i = 0; i < target.size(); i++)
		{
			if (_mapOldNew[i] != -1)
			{
				result.data()[_mapOldNew[i]] = target.data()[i];
			}
		}
	}

	void FemSimulator::vectorInsertRows(VectorX & target, VectorX & result)
	{
		result.setZero();
		for (int i = 0; i < result.size(); i++)
		{
			if (_mapOldNew[i] != -1)
			{
				result.data()[i] = target.data()[_mapOldNew[i]];
			}
		}
	}

	void FemSimulator::resetSparseMatrix(SparseMatrix& sparseMatrix)
	{
		for (int i = 0; i < sparseMatrix.outerSize(); ++i)
			for (SparseMatrix::InnerIterator it(sparseMatrix, i); it; ++it)
			{
				it.valueRef() = 0;
			}
	}

	void FemSimulator::getSparseMatrixTopologyTYPE(SparseMatrix& sparesMatrixTopology)
	{
		std::vector<TripletX> entrys;

		for (size_t mid = 0; mid < models.size(); mid++)
		{
			BaseModel* m = getModel(mid);
			for (size_t eid = 0; eid < m->tetElementNum; eid++)
			{
				Vector4i indices = m->getTetElementOverall(eid);
				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						for (int k = 0; k < 3; k++)
							for (int l = 0; l < 3; l++)
								entrys.push_back(TripletX(3 * indices.data()[i] + k, 3 * indices.data()[j] + l, 0.0));
					}
				}
			}
		}

		sparesMatrixTopology.resize(_sysDim, _sysDim);
		sparesMatrixTopology.setFromTriplets(entrys.begin(), entrys.end());
	}

}