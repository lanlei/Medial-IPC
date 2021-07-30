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
		doTime(frame);
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

	void FemSimulator::doTime(int frame){}

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