#include "MipcSimulator.h"

bool  MIPC::MipcSimulator::addModelFromConfigFile(const std::string filename, TiXmlElement * item)
{
	MipcModel* m = new MipcModel();
	bool isReadMesh = BaseSimulator::addModelFromConfigFile(filename, item, m);
	if (!isReadMesh)
		return false;
	m = getModel(models.size() - 1);
	qeal scale = 1.0;
	qeal tx = 0.0, ty = 0.0, tz = 0.0;
	qeal rx = 0.0, ry = 1.0, rz = 0.0, rsita = 0.0;
	qeal density = 1.0;
	qeal youngModulus = 300.0, poisson = 0.35;
	int enableGravity = 1;
	bool locality = false;
	qeal localityThreshold = 0.0;
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
		else if (itemName == std::string("YoungModulus"))
		{
			std::string str = childItem->GetText();
			ss << str;
			ss >> youngModulus;
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
		else if (itemName == std::string("localityWeight"))
		{
			std::string str = childItem->GetText();
			ss << str;
			ss >> localityThreshold;
			if (localityThreshold != 0.0) locality = true;
		}

		childItem = childItem->NextSiblingElement();
	}
	m->scaleModel(scale);
	m->translateModel(tx, ty, tz);
	m->rotateModel(rx, ry, rz, rsita);
	m->enableGravityForce(enableGravity);
	m->setWeightLocality = locality;
	m->weightLocalityThreshold = localityThreshold;
	m->initMeshesHandel();
	ENuMaterial * material = downcastENuMaterial(m->getTetMeshHandle()->getElementMaterialById(0));
	material->setDensity(density);
	material->setE(youngModulus);
	material->setNu(poisson);
	return true;
}

void MIPC::MipcSimulator::doTimeGpuDenseSystem(int frame)
{
	int newton_iter = 0;
	qeal dt = _timeStep;

	cudaMemcpy(_devSysXn, _devSysX, _sysDim * sizeof(qeal), cudaMemcpyDeviceToDevice);
	cudaMemcpy(_devReducedXn, _devReducedX, _sysReducedDim * sizeof(qeal), cudaMemcpyDeviceToDevice);
	cudaMemcpy(_devReducedVn, _devReducedVelocity, _sysReducedDim * sizeof(qeal), cudaMemcpyDeviceToDevice);

	// compute fullspace/reduced external force
	cudaMemcpy(_devExternalForce, _sysGravityForce.data(), _sysGravityForce.size() * sizeof(qeal), cudaMemcpyHostToDevice);

	for (int mid = 0; mid < models.size(); mid++)
	{
		for (int idx = 0; idx < 3; idx++)
		{
			cublasDgemv(blasHandle, CUBLAS_OP_T, _hostTetPointProjectionRowsXYZ[mid], _hostTetPointProjectionColsXYZ[mid], &cublas_pos_one, _devTetPointProjectionXYZ[mid], _hostTetPointProjectionRowsXYZ[mid], _devExternalForce + _hostTetPointXOffset[mid] + idx, 3, &cublas_zero, _devReducedExternalForce + _hostFrameBufferOffset[mid] + idx, 3);
			cudaDeviceSynchronize();
		}
	}

	// compute predictive pos
	updatePredictivePos
	(
		_sysReducedDim,
		_devReducedDim,
		_devReducedXn,
		_devReducedVelocity,
		_devTimeStep,
		_devReducedXtilde
	);
	for (int mid = 0; mid < models.size(); mid++)
	{
		for (int idx = 0; idx < 3; idx++)
		{
			cublasDgemv(blasHandle, CUBLAS_OP_N, _hostTetPointProjectionRowsXYZ[mid], _hostTetPointProjectionColsXYZ[mid], &cublas_pos_one, _devTetPointProjectionXYZ[mid], _hostTetPointProjectionRowsXYZ[mid], _devReducedXtilde + _hostFrameBufferOffset[mid] + idx, 3, &cublas_zero, _devSysXtilde + _hostTetPointXOffset[mid] + idx, 3);
			cudaDeviceSynchronize();
		}
	}

	constructConstraintSet(_kappa, true);
	
	qeal Ep = computeEnergy(_devSysX, _devSysXtilde);
	// compute ipc constraint set
	do
	{
		computeElasticsHessianAndGradient(_sysReducedRhs.data(), _sysReducedMatrix.data());

		for (int i = 0; i < _activeCollisionEvents.size(); i++)
			_activeCollisionEvents[i]->getGradientAndHessian(_kappa, _sysReducedRhs, _sysReducedMatrix);

		for (int i = 0; i < _frictionCollisionEvents.size(); i++)
			_frictionCollisionEvents[i]->getFrictionGradientAndHessian(_kappa, _sysReducedRhs, _sysReducedMatrix);

		// solve
		cudaMemcpy(_devReducedMatrix, _sysReducedMatrix.data(), _sysReducedDim * _sysReducedDim * sizeof(qeal), cudaMemcpyHostToDevice);
		cudaMemcpy(_devReducedRhs, _sysReducedRhs.data(), _sysReducedDim * sizeof(qeal), cudaMemcpyHostToDevice);
		computeSystemUsingCusolverDenseChol(_devReducedMatrix, _devReducedRhs, _devReducedDir, _sysReducedDim);
		cudaMemcpy(_sysReducedDir.data(), _devReducedDir, _sysReducedDim * sizeof(qeal), cudaMemcpyDeviceToHost);

		for (int mid = 0; mid < models.size(); mid++)
		{
			for (int idx = 0; idx < 3; idx++)
			{
				cublasDgemv(blasHandle, CUBLAS_OP_N, _hostTetPointProjectionRowsXYZ[mid], _hostTetPointProjectionColsXYZ[mid], &cublas_pos_one, _devTetPointProjectionXYZ[mid], _hostTetPointProjectionRowsXYZ[mid], _devReducedDir + _hostFrameBufferOffset[mid] + idx, 3, &cublas_zero, _devDir + _hostTetPointXOffset[mid] + idx, 3);
				cudaDeviceSynchronize();
			}
		}
		cudaMemcpy(_sysDir.data(), _devDir, _sysDim * sizeof(qeal), cudaMemcpyDeviceToHost);

		qeal res = _sysDir.cwiseAbs().maxCoeff() / dt;	
		if (newton_iter > 0 && res <= tol)
		{
			std::cout << "Frame " << frame << " converges to " << res << " after " << newton_iter <<" iters."<< std::endl;
			break;
		}
		else
		{
			std::cout << "--  Frame: " << frame << " iter: " << newton_iter << "  res: " << res << std::endl;
		}

		// line search
		cudaMemcpy(_devMedialPointPosition, medialPointsBuffer.buffer.data(), 3 * totalMedialPoinsNum * sizeof(qeal), cudaMemcpyHostToDevice);
		computeMedialPointsMovingDir
		(
			totalMedialPoinsNum,
			_devMedialPointsNum,
			_devMedialOriPointPosition,
			_devReducedDir,
			_devMedialPointMovingDir
		);

		qeal toi = 1.0;	
		getToI(toi);

		qeal E; 	
		cudaMemcpy(_devSearchReducedX, _devReducedX, _sysReducedDim * sizeof(qeal), cudaMemcpyDeviceToDevice);
		do
		{
			updateLineSearchX
			(
				_sysReducedDim,
				_devReducedDim,
				_devSearchReducedX,
				_devReducedDir,
				toi,
				_devReducedX
			);
			cudaMemcpy(_sysReducedX.data(), _devReducedX, _sysReducedDim * sizeof(qeal), cudaMemcpyDeviceToHost);

			for (int mid = 0; mid < models.size(); mid++)
			{
				for (int idx = 0; idx < 3; idx++)
				{
					cublasDgemv(blasHandle, CUBLAS_OP_N, _hostTetPointProjectionRowsXYZ[mid], _hostTetPointProjectionColsXYZ[mid], &cublas_pos_one, _devTetPointProjectionXYZ[mid], _hostTetPointProjectionRowsXYZ[mid], _devReducedX + _hostFrameBufferOffset[mid] + idx, 3, &cublas_zero, _devSysX + _hostTetPointXOffset[mid] + idx, 3);
					cudaDeviceSynchronize();
				}
			}
			constructConstraintSet(_kappa, false);
			E = computeEnergy(_devSysX, _devSysXtilde);
			std::cout << "----  line search  Ep: " << E <<"  E: " <<E <<"  toi: " << toi << std::endl;
			toi *= 0.5;
		} while ((E - Ep) > MIN_VALUE);
		Ep = E;

	} while (++newton_iter);


	updatedVelocity
	(
		_sysReducedDim,
		_devReducedDim,
		_devReducedX,
		_devReducedXtilde,
		_devTimeStep,
		_devReducedVelocity
	);

	for (int mid = 0; mid < models.size(); mid++)
	{
		for (int idx = 0; idx < 3; idx++)
		{
			cublasDgemv(blasHandle, CUBLAS_OP_N, _hostTetPointProjectionRowsXYZ[mid], _hostTetPointProjectionColsXYZ[mid], &cublas_pos_one, _devTetPointProjectionXYZ[mid], _hostTetPointProjectionRowsXYZ[mid], _devReducedVelocity + _hostFrameBufferOffset[mid] + idx, 3, &cublas_zero, _devSysVelocity + _hostTetPointXOffset[mid] + idx, 3);
			cudaDeviceSynchronize();
		}
	}
	cudaMemcpy(_sysX.data(), _devSysX, _sysDim * sizeof(qeal), cudaMemcpyDeviceToHost);

}

void MIPC::MipcSimulator::doTimeGpuSparseSystem(int frame)
{

}

qeal MIPC::MipcSimulator::computeEnergy(qeal* devXn, qeal* devXtilde)
{
	// compute static force energy
	qeal e0;
	cublasDdot(blasHandle, _sysDim, _devExternalForce, 1, devXn, 1, &e0);
	e0 *= -1.0 * _timeStep * _timeStep;

	//compute inertial energy
	computeInertialDiffX
	(
		_sysDim,
		_devDim,
		devXn,
		devXtilde,
		_devTetPointsMass,
		_devDiffX,
		_devMassMulDiffX
	);
	qeal e1;
	cublasDdot(blasHandle, _sysDim, _devDiffX, 1, _devMassMulDiffX, 1, &e1);
	e1 *= 0.5;
	
	// compute elastics energy
	assembleTetELementX
	(
		totalTetElementNum,
		_devTotalTetElementNum,
		_devTetElementIndices,
		_devSysX,
		_devTetElementX
	);
	computeElementsEnergy
	(
		totalTetElementNum,
		_devTotalTetElementNum,
		_devTetElementX,
		_devTetElementDm,
		_devTetElementInvDm,
		_devTetElementAttri,
		_devTetElementVol,
		_devTimeStep,
		_devTetElementPotentialEnergy
	);
	qeal e2;
	cublasDdot(blasHandle, totalTetElementNum, _devTetElementPotentialEnergy, 1, _devTetElementSizeOne, 1, &e2);

	qeal e3;
	e3 = 0.0;
	for (int i = 0; i < _activeCollisionEvents.size(); i++)
		e3 += _activeCollisionEvents[i]->getEnergy(_kappa);
	e3 *= _timeStep * _timeStep;

	qeal e4;
	e4 = 0.0;
	for (int i = 0; i < _frictionCollisionEvents.size(); i++)
		e4 += _frictionCollisionEvents[i]->frictionEnergy();
	e4 *= _timeStep * _timeStep;
//	std::cout << e0 << " " << e1 << " " << e2 << " " << e3 <<" " << _activeCollisionEvents.size() << std::endl;
	return e0 + e1 + e2 + e3 + e4;
}

void MIPC::MipcSimulator::computeElasticsHessianAndGradient(qeal * elasticsDerivative, qeal * elasticsHessian)
{
	qeal timeStep2 = _timeStep * _timeStep;
	assembleTetELementX
	(
		totalTetElementNum,
		_devTotalTetElementNum,
		_devTetElementIndices,
		_devSysX,
		_devTetElementX
	);

	computeTetElementInternalForce
	(
		totalTetElementNum,
		_devTotalTetElementNum,
		_devTetElementX,
		_devTetElementDm,
		_devTetElementInvDm,
		_devTetElementdFPK,
		_devTetElementAttri,
		_devTetElementForce
	);

	computeTetElementStiffness
	(
		totalTetElementNum,
		_devTotalTetElementNum,
		_devTetElementX,
		_devTetElementDm,
		_devTetElementInvDm,
		_devTetElementdFdu,
		_devTetElementdFPK,// as eigen value of dPdF
		_devTetElementdPdF,// as eigen vector of dPdF
		_devTetElementAttri,
		_devTetElementStiffness
	);

	assembleTetPointsForceFromElementForce
	(
		totalTetPointsNum,
		_devTotalTetPointsNum,
		_devTetPointsSharedElementNum,
		_devTetPointsSharedElementOffset,
		_devTetPointsSharedElementList,
		_devTetElementForce,
		_devInternalForce
	);

	for (int mid = 0; mid < models.size(); mid++)
	{
		for (int idx = 0; idx < 3; idx++)
		{
			cublasDgemv(blasHandle, CUBLAS_OP_T, _hostTetPointProjectionRowsXYZ[mid], _hostTetPointProjectionColsXYZ[mid], &cublas_pos_one, _devTetPointProjectionXYZ[mid], _hostTetPointProjectionRowsXYZ[mid], _devInternalForce + _hostTetPointXOffset[mid] + idx, 3, &cublas_zero, _devReducedInternalForce + _hostFrameBufferOffset[mid] + idx, 3);
			cudaDeviceSynchronize();
		}
	}

	assembleReducedStiffness
	(
		_hostAssembleBlockNum,
		_devAssembleBlockNum,
		_devAssembleBlockIndex,
		_devStiffnessBlockSharedTetElementList,
		_devStiffnessBlockSharedTetElementNum,
		_devStiffnessBlockSharedTetElementOffset,
		_devPojectionStiffnessList,
		_devTetElementSharedFrameList,
		_devTetElementSharedFrameOffset,
		_devTetElementFrameProjectionBuffer,
		_devTetElementFrameProjectionNum,
		_devTetElementFrameProjectionOffset,
		_devTetElementStiffness,
		_devReducedDim,
		_devReducedStiffness
	);

	cudaMemcpy(_devReducedMatrix, _devReducedMassMatrix, _sysReducedDim * _sysReducedDim * sizeof(qeal), cudaMemcpyDeviceToDevice);
	cublasDaxpy(blasHandle, _sysReducedDim * _sysReducedDim, &timeStep2, _devReducedStiffness, 1, _devReducedMatrix, 1);

	computeReducedInertia
	(
		_sysReducedDim,
		_devReducedDim,
		_devReducedX,
		_devReducedXn,
		_devReducedVn,
		_devTimeStep,
		_devInertia
	);

	cublasDgemv(blasHandle, CUBLAS_OP_N, _sysReducedDim, _sysReducedDim, &cublas_pos_one, _devReducedMassMatrix, _sysReducedDim, _devInertia, 1, &cublas_zero, _devReducedRhs, 1);
	cudaDeviceSynchronize();

	cublasDaxpy(blasHandle, _sysReducedDim, &timeStep2, _devReducedInternalForce, 1, _devReducedRhs, 1);

	timeStep2 *= -1.0;
	cublasDaxpy(blasHandle, _sysReducedDim, &timeStep2, _devReducedExternalForce, 1, _devReducedRhs, 1);

	qeal alpha = -1.0;
	cublasDscal(blasHandle, _sysReducedDim, &alpha, _devReducedRhs, 1);

	cudaMemcpy(elasticsDerivative, _devReducedRhs, _sysReducedDim * sizeof(qeal), cudaMemcpyDeviceToHost);
	cudaMemcpy(elasticsHessian, _devReducedMatrix, _sysReducedDim * _sysReducedDim * sizeof(qeal), cudaMemcpyDeviceToHost);
}

void MIPC::MipcSimulator::constructConstraintSet(const qeal kappa, bool updateFriction)
{
	// update collision set
	for (int i = 0; i < _reducedFrameList.size(); i++)
		_reducedFrameList[i]->transform();
	_activeCollisionEvents.clear();
	if (updateFriction && enableFriction)
		_frictionCollisionEvents.clear();

	for (int i = 0; i < _overallCollisionEvents.size(); i++)
	{
		_overallCollisionEvents[i]->computeDistance();
		if (_overallCollisionEvents[i]->isActive())
		{
			_activeCollisionEvents.push_back(_overallCollisionEvents[i]);
			if (updateFriction && enableFriction)
			{
				_overallCollisionEvents[i]->computeLagTangentBasis(kappa);
				_frictionCollisionEvents.push_back(_overallCollisionEvents[i]);
			}
		}
	}
}

void MIPC::MipcSimulator::getToI(qeal& toi)
{
	if (_hostCollisionEventNum == 0)
	{
		toi = 1.0;
		return;
	}

	MPsCCD
	(
		_hostCollisionEventNum,
		_devCollisionEventNum,
		_devMedialPointPosition,
		_devMedialPointRadius,
		_devStaticMedialPointPosition,
		_devStaticMedialPointRadius,
		_devMedialPointMovingDir,
		_devCollisionEventList,
		_devCCD
	);

	int idx = 0;
	cublasIdamin(blasHandle, _hostCollisionEventNum,
		_devCCD, 1, &idx);
	cudaMemcpy(&toi, _devCCD + (idx - 1), sizeof(qeal), cudaMemcpyDeviceToHost);

	if (toi < 1.0)
		toi *= 0.8;
}

void MIPC::MipcSimulator::computeSystemUsingCusolverDenseChol(qeal* devSys, qeal* devRhs, qeal* devX, int dim)
{
	int bufferSize = 0;
	cusolverDnDpotrf_bufferSize(dnHandle, dnUplo, dim, devSys, dim, &bufferSize);

	double* devBuffer = NULL;

	int* devDnInfo;
	cudaMalloc((void**)&devDnInfo, sizeof(int));
	cudaMalloc(&devBuffer, bufferSize * sizeof(double));
	cudaMemset(devDnInfo, 0, sizeof(int));

	cusolverStatus_t satus = cusolverDnDpotrf(dnHandle, dnUplo, dim, devSys, dim, devBuffer, bufferSize, devDnInfo);

	int host_info;
	cudaMemcpy(&host_info, devDnInfo, sizeof(int), cudaMemcpyDeviceToHost);

	if (host_info != 0 || satus != CUSOLVER_STATUS_SUCCESS)
	{
		fprintf(stderr, "Error: Cholesky factorization failed\n");
		printf("%d\n", host_info);
		Eigen::LDLT<MatrixX> ldlt;
		ldlt.compute(_sysReducedMatrix);
		VectorX sdf = ldlt.solve(_sysReducedRhs);
		cudaMemcpy(_devReducedXtilde, sdf.data(), _sysReducedDim * sizeof(qeal), cudaMemcpyHostToDevice);

		cudaFree(devBuffer);
		cudaFree(devDnInfo);
		return;
	}
	cudaMemcpy(devX, devRhs, dim * sizeof(qeal), cudaMemcpyDeviceToDevice);
	cusolverDnDpotrs(dnHandle, dnUplo, dim, 1, devSys, dim, devX, dim, devDnInfo);
	cudaDeviceSynchronize();
	cudaFree(devBuffer);
	cudaFree(devDnInfo);
}

void MIPC::MipcSimulator::computeSystemUsingCusolverSparseChol(qeal* devSys, int* devSysRowPtr, int* devSysColInd, int nnz, qeal* devRhs, qeal* devX, int dim)
{
	qeal tol = 1.e-12;
	const int reorder = 0; /* no reordering */
	int singularity = 0; /* -1 if A is invertible under tol. */

	cusolverStatus_t  t = cusolverSpDcsrlsvchol(
		cusolverSpH, dim, nnz,
		spdescrA, devSys, devSysRowPtr, devSysColInd,
		devRhs, tol, reorder, devX, &singularity);

	cudaDeviceSynchronize();
	if (0 <= singularity)
	{
		printf("WARNING: the matrix is singular at row %d under tol (%E)\n", singularity, tol);
	}
}

void MIPC::MipcSimulator::initialization()
{
	FemSimulator::initialization();
	_sysReducedDim = 0;
	_sysReducedOffsetByModel.resize(models.size());

	_staticFramesNum = 0;
	_linearFramesNum = 0;
	_quadraticFramesNum = 0;
	for (size_t mid = 0; mid < models.size(); mid++)
	{
		MipcModel* m = getModel(mid);
		std::string frameFilename = m->dir + "frames.dofs";
		m->readFrameMatList(frameFilename);
		_sysReducedOffsetByModel[mid] = _sysReducedDim;
		_sysReducedDim += m->getReducedDim();
		for (size_t i = 0; i < m->getNonStaticFramesNum(); i++)
		{
			int mvid = m->getNonStaticFrameMedialId(i);
			int gmvid = m->getMedialPointOverallId(mvid);
			_sysNonStaticFrameOverallId.push_back(gmvid);
		}

		_staticFramesNum += m->getStaticFramesNum();
		_linearFramesNum += m->getLinearFramesNum();
		_quadraticFramesNum += m->getQuadraticFramesNum();
		_translationFramesNum += m->getTranslationFramesNum();
	}
	_nonStaticFramesNum = _linearFramesNum + _quadraticFramesNum + _translationFramesNum;

	_sysReducedMatrix.resize(_sysReducedDim, _sysReducedDim);
	_sysReducedRhs.resize(_sysReducedDim);
	_sysReducedX.resize(_sysReducedDim);
	_sysReducedX.setZero();
	_sysReducedXtilde.resize(_sysReducedDim);
	_sysReducedXtilde.setZero();
	_sysReducedDir.resize(_sysReducedDim);
	_sysReducedDir.setZero();

	_sysReducedVelocity.resize(_sysReducedDim);
	_sysReducedVelocity.setZero();
	_sysReducedAccrelation.resize(_sysReducedDim);
	_sysReducedAccrelation.setZero();

	_sysReducedXn.resize(_sysReducedDim);
	_sysReducedXn.setZero();
	_sysReducedVn.resize(_sysReducedDim);
	_sysReducedVn.setZero();
	_sysReducedAccn.resize(_sysReducedDim);
	_sysReducedAccn.setZero();

	_sysReducedExternalForce.resize(_sysReducedDim);
	_sysReducedExternalForce.setZero();
	_sysReducedInternalForce.resize(_sysReducedDim);
	_sysReducedInternalForce.setZero();

	_sysReducedMass.resize(_sysReducedDim, _sysReducedDim);
	_sysReducedStiffness.resize(_sysReducedDim, _sysReducedDim);

	int bufferOffset = 0;
	int frameOffset = 0;
	for (size_t mid = 0; mid < models.size(); mid++)
		getModel(mid)->createReducedFrame(frameOffset, bufferOffset, _sysReducedX.data(), _sysReducedXn.data(), _sysReducedXtilde.data(), _sysReducedVelocity.data(), _sysReducedVn.data(), _sysReducedAccrelation.data(), _sysReducedAccn.data(), _reducedFrameList);

	_sysReducedSparseProjection.resize(getSysDimension(), _sysReducedDim);

	std::vector<TripletX> triplet;
	for (size_t mid = 0; mid < models.size(); mid++)
	{
		MipcModel* m = getModel(mid);
		int rowBufferOffset = m->tetPoints.offset;
		int colBufferOffset = _sysReducedOffsetByModel[mid];

		m->getGlobalSparseProjectionMatrixTriplet(rowBufferOffset, colBufferOffset, triplet);
	}

	_sysReducedSparseProjection.setFromTriplets(triplet.begin(), triplet.end());
	_sysReducedSparseProjectionT = _sysReducedSparseProjection.transpose();
	_sysReducedMass = _sysReducedSparseProjectionT * _sysMassMatrix * _sysReducedSparseProjection;
	_sysMassMatrix.resize(0, 0);
	_sysMatrix.resize(0, 0);
	if (_mu > 0)enableFriction = true;
	tol = 1e-3 * _diagLen;

	genOverallCollisionEvents();
	initForGpu();
}

void MIPC::MipcSimulator::run(int frame)
{
	if (_sysMatType == SPARSE)
		doTimeGpuSparseSystem(frame);
	else doTimeGpuDenseSystem(frame);
}

void MIPC::MipcSimulator::postRun()
{
	_sysCurrentPosition = _sysOriginalPosition + _sysX;
	std::copy(_sysCurrentPosition.data(), _sysCurrentPosition.data() + _sysCurrentPosition.size(), tetPointsBuffer.buffer.data());
	alignAllMesh(tetPointsBuffer.buffer.data());
}

void MIPC::MipcSimulator::initForGpu()
{
	if (_runPlatform != RunPlatform::CUDA)
		return;
	gpuSize = 0;
	std::cout << "init gpu buffer" << std::endl;
	std::cout << "  -- general memory" << std::endl;
	initCudaGeneralSysMemory();
	std::cout << "  -- tet mesh memory" << std::endl;
	initCudaTetMeshMemory();
	std::cout << "  -- medial mesh memory" << std::endl;
	initCudaMedialMeshMemory();
	std::cout << "  -- reduced projection memory" << std::endl;
	initCudaReducedProjectionMemory();
	if (_sysMatType == SPARSE)
	{
		std::cout << "  -- sparse solver memory" << std::endl;
		initCudaSparseSysMemory();
	}
	std::cout << "  -- collision event memory" << std::endl;
	initCudaCollisionMemory();
	std::cout << "  -- cost gpu memory: " << gpuSize / 1024 / 1024 << " Mb" << std::endl;

}

void MIPC::MipcSimulator::initCudaTetMeshMemory()
{
	CUDA_CALL(cudaMalloc((void**)&_devTotalTetElementNum, sizeof(int))); gpuSize += sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTotalTetElementNum, &totalTetElementNum, sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTotalTetPointsNum, sizeof(int))); gpuSize += sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTotalTetPointsNum, &totalTetPointsNum, sizeof(int), cudaMemcpyHostToDevice));

	std::vector<int> hostTetElementIndices(4 * totalTetElementNum);
	int hostTetElementCount = 0;
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->tetElementNum; j++)
		{
			Vector4i indices = models[i]->getTetElement(j);
			hostTetElementIndices[4 * hostTetElementCount] = models[i]->getTetPointOverallId(indices.data()[0]);
			hostTetElementIndices[4 * hostTetElementCount + 1] = models[i]->getTetPointOverallId(indices.data()[1]);
			hostTetElementIndices[4 * hostTetElementCount + 2] = models[i]->getTetPointOverallId(indices.data()[2]);
			hostTetElementIndices[4 * hostTetElementCount + 3] = models[i]->getTetPointOverallId(indices.data()[3]);
			hostTetElementCount++;
		}
	}

	CUDA_CALL(cudaMalloc((void**)&_devTetElementIndices, 4 * totalTetElementNum * sizeof(int))); gpuSize += 4 * totalTetElementNum * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTetElementIndices, hostTetElementIndices.data(), 4 * totalTetElementNum * sizeof(int), cudaMemcpyHostToDevice));


	std::vector<int> hostTetPointsSharedElementNum(totalTetPointsNum);
	std::vector<int> hostTetPointsSharedElementOffset(totalTetPointsNum);
	std::vector<std::vector<int>> hostTetPointsSharedElementList(totalTetPointsNum);

	for (int i = 0; i < models.size(); i++)
	{
		BaseModel* m = models[i];
		for (int j = 0; j < m->tetElementNum; j++)
		{
			int geleId = m->getTetElementOverallId(j);
			Vector4i indices = m->getTetElement(j);
			for (int k = 0; k < 4; k++)
			{
				int gvid = m->getTetPointOverallId(indices.data()[k]);
				hostTetPointsSharedElementList[gvid].push_back(geleId);
				hostTetPointsSharedElementList[gvid].push_back(k);
			}
		}
	}

	std::vector<int> flatTetPointsSharedElementList;
	for (int i = 0; i < totalTetPointsNum; i++)
	{
		int num = hostTetPointsSharedElementList[i].size() / 2;
		hostTetPointsSharedElementNum[i] = num;
		int offset = flatTetPointsSharedElementList.size();
		hostTetPointsSharedElementOffset[i] = offset;
		for (int j = 0; j < num; j++)
		{
			flatTetPointsSharedElementList.push_back(hostTetPointsSharedElementList[i][2 * j]);
			flatTetPointsSharedElementList.push_back(hostTetPointsSharedElementList[i][2 * j + 1]);
		}
	}

	CUDA_CALL(cudaMalloc((void**)&_devTetPointsSharedElementNum, totalTetPointsNum * sizeof(int))); gpuSize += totalTetPointsNum * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTetPointsSharedElementNum, hostTetPointsSharedElementNum.data(), totalTetPointsNum * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)& _devTetPointsSharedElementOffset, totalTetPointsNum * sizeof(int))); gpuSize += totalTetPointsNum * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTetPointsSharedElementOffset, hostTetPointsSharedElementOffset.data(), totalTetPointsNum * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTetPointsSharedElementList, flatTetPointsSharedElementList.size() * sizeof(int))); gpuSize += flatTetPointsSharedElementList.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTetPointsSharedElementList, flatTetPointsSharedElementList.data(), flatTetPointsSharedElementList.size() * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTetElementX, 12 * totalTetElementNum * sizeof(qeal))); gpuSize += 12 * totalTetElementNum * sizeof(qeal);

	CUDA_CALL(cudaMalloc((void**)&_devTetElementPotentialEnergy, totalTetElementNum * sizeof(qeal))); gpuSize += totalTetElementNum * sizeof(qeal);

	CUDA_CALL(cudaMalloc((void**)&_devTetElementSizeOne, totalTetElementNum * sizeof(qeal))); gpuSize += totalTetElementNum * sizeof(qeal);
	VectorX hostTetElementSizeOne(totalTetElementNum);
	hostTetElementSizeOne.setOnes();
	CUDA_CALL(cudaMemcpy(_devTetElementSizeOne, hostTetElementSizeOne.data(), totalTetElementNum * sizeof(qeal), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTetElementForce, 12 * totalTetElementNum * sizeof(qeal))); gpuSize += 12 * totalTetElementNum * sizeof(qeal);

	CUDA_CALL(cudaMalloc((void**)&_devTetElementStiffness, 144 * totalTetElementNum * sizeof(qeal))); gpuSize += 144 * totalTetElementNum * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devTetElementdPdF, 81 * totalTetElementNum * sizeof(qeal))); gpuSize += 81 * totalTetElementNum * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devTetElementdFPK, 9 * totalTetElementNum * sizeof(qeal))); gpuSize += 9 * totalTetElementNum * sizeof(qeal);

	std::vector<qeal> tempBuffer0(9 * totalTetElementNum);
	std::vector<qeal> tempBuffer1(9 * totalTetElementNum);
	for (int mid = 0; mid < models.size(); mid++)
	{
		MipcModel* m = getModel(mid);
		for (int eid = 0; eid < m->tetElementNum; eid++)
		{
			int geid = m->getTetElementOverallId(eid);
			BaseTetElementParam* para = m->getTetMeshHandle()->getTetElementParam(eid);
			Matrix3 Dm, invDm;
			Dm = para->Dm;
			invDm = para->invDm;
			std::copy(Dm.data(), Dm.data() + 9, tempBuffer0.data() + 9 * geid);
			std::copy(invDm.data(), invDm.data() + 9, tempBuffer1.data() + 9 * geid);
		}
	}

	CUDA_CALL(cudaMalloc((void**)&_devTetElementDm, tempBuffer0.size() * sizeof(qeal))); gpuSize += tempBuffer0.size() * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devTetElementDm, tempBuffer0.data(), tempBuffer0.size() * sizeof(qeal), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTetElementInvDm, tempBuffer1.size() * sizeof(qeal))); gpuSize += tempBuffer1.size() * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devTetElementInvDm, tempBuffer1.data(), tempBuffer1.size() * sizeof(qeal), cudaMemcpyHostToDevice));

	tempBuffer0.resize(9 * 12 * totalTetElementNum);
	tempBuffer1.resize(3 * totalTetElementNum);
	for (int mid = 0; mid < models.size(); mid++)
	{
		MipcModel* m = getModel(mid);
		for (int eid = 0; eid < m->tetElementNum; eid++)
		{
			int geid = m->getTetElementOverallId(eid);
			BaseTetElementParam* para = m->getTetMeshHandle()->getTetElementParam(eid);
			MatrixX dFdu = m->getTetMeshHandle()->getTetElementParam(eid)->dFdu;
			std::copy(dFdu.data(), dFdu.data() + 108, tempBuffer0.data() + 108 * geid);

			qeal volume = para->volume;
			ENuMaterial * material = downcastENuMaterial(m->getTetMeshHandle()->getElementMaterial(eid));

			qeal lameMiu = material->getMu();
			qeal lameLamda = material->getLambda();
			tempBuffer1[3 * geid] = volume;
			tempBuffer1[3 * geid + 1] = lameMiu;
			tempBuffer1[3 * geid + 2] = lameLamda;
		}
	}

	CUDA_CALL(cudaMalloc((void**)&_devTetElementdFdu, tempBuffer0.size() * sizeof(qeal))); gpuSize += tempBuffer0.size() * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devTetElementdFdu, tempBuffer0.data(), tempBuffer0.size() * sizeof(qeal), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTetElementAttri, tempBuffer1.size() * sizeof(qeal))); gpuSize += tempBuffer1.size() * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devTetElementAttri, tempBuffer1.data(), tempBuffer1.size() * sizeof(qeal), cudaMemcpyHostToDevice));

	std::vector<qeal> tetElementVol(totalTetElementNum);
	for (size_t i = 0; i < models.size(); i++)
	{
		BaseModel* m = getModel(i);
		for (size_t eid = 0; eid < m->tetElementNum; eid++)
		{
			int geid = m->getTetElementOverallId(eid);
			qeal v = m->getTetMeshHandle()->getTetElementParam(eid)->volume;
			tetElementVol[geid] = v;
		}
	}

	CUDA_CALL(cudaMalloc((void**)&_devTetElementVol, totalTetElementNum * sizeof(qeal))); gpuSize += totalTetElementNum * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devTetElementVol, tetElementVol.data(), tetElementVol.size() * sizeof(qeal), cudaMemcpyHostToDevice));
}

void MIPC::MipcSimulator::initCudaMedialMeshMemory()
{
	CUDA_CALL(cudaMalloc((void**)&_devMedialPointsNum, sizeof(int))); gpuSize += sizeof(int);
	CUDA_CALL(cudaMemcpy(_devMedialPointsNum, &totalMedialPoinsNum, sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMalloc((void**)&_devMedialOriPointPosition, 3 * totalMedialPoinsNum * sizeof(qeal))); gpuSize += 3 * totalMedialPoinsNum * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devMedialOriPointPosition, medialPointsBuffer.buffer.data(), 3 * totalMedialPoinsNum * sizeof(qeal), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMalloc((void**)&_devMedialPointPosition, 3 * totalMedialPoinsNum * sizeof(qeal))); gpuSize += 3 * totalMedialPoinsNum * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devMedialPointPosition, medialPointsBuffer.buffer.data(), 3 * totalMedialPoinsNum * sizeof(qeal), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devMedialPointRadius, totalMedialPoinsNum * sizeof(qeal))); gpuSize += totalMedialPoinsNum * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devMedialPointRadius, medialRadiusBuffer.buffer.data(), totalMedialPoinsNum * sizeof(qeal), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devStaticMedialPointPosition, 3 * staticModelPool.totalMedialPoinsNum * sizeof(qeal))); gpuSize += 3 * staticModelPool.totalMedialPoinsNum * sizeof(qeal);

	CUDA_CALL(cudaMemcpy(_devStaticMedialPointPosition, staticModelPool.medialPointsBuffer.buffer.data(), 3 * staticModelPool.totalMedialPoinsNum * sizeof(qeal), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devStaticMedialPointRadius, staticModelPool.totalMedialPoinsNum * sizeof(qeal))); gpuSize += staticModelPool.totalMedialPoinsNum * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devStaticMedialPointRadius, staticModelPool.medialRadiusBuffer.buffer.data(), staticModelPool.totalMedialPoinsNum * sizeof(qeal), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devMedialPointMovingDir, 3 * totalMedialPoinsNum * sizeof(qeal))); gpuSize += 3 * totalMedialPoinsNum * sizeof(qeal);
}

void MIPC::MipcSimulator::initCudaGeneralSysMemory()
{
	CUDA_CALL(cudaMalloc((void**)&_devDim, sizeof(int))); gpuSize += sizeof(int);
	CUDA_CALL(cudaMemcpy(_devDim, &_sysDim, sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devReducedDim, sizeof(int))); gpuSize += sizeof(int);
	CUDA_CALL(cudaMemcpy(_devReducedDim, &_sysReducedDim, sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTimeStep, sizeof(qeal))); gpuSize += sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devTimeStep, &_timeStep, sizeof(qeal), cudaMemcpyHostToDevice));

	// gpu fullspace sys vector & matrix
	VectorX sysZeroVec(_sysDim);
	sysZeroVec.setZero();
	CUDA_CALL(cudaMalloc((void**)&_devSysX, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devSysXn, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devSysXtilde, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devDiffX, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devSysVelocity, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);
	cudaMemcpy(_devSysVelocity, sysZeroVec.data(), sizeof(qeal) *_sysDim, cudaMemcpyHostToDevice);
	CUDA_CALL(cudaMalloc((void**)&_devMassMulDiffX, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devExternalForce, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devInternalForce, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devDir, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devSearchX, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);

	VectorX massVector(_sysDim);
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->tetPointsNum; j++)
		{
			int pointId = models[i]->getTetPointOverallId(j);
			qeal density = models[i]->getTetMeshHandle()->getNodeMaterial(j)->getDensity();
			qeal volume = models[i]->getTetMeshHandle()->getTetNodeParam(j)->volume;

			qeal v = density * volume;
			massVector.data()[3 * pointId] = v;
			massVector.data()[3 * pointId + 1] = v;
			massVector.data()[3 * pointId + 2] = v;
		}
	}
	CUDA_CALL(cudaMalloc((void**)&_devTetPointsMass, _sysDim * sizeof(qeal))); gpuSize += _sysDim * sizeof(qeal);
	cudaMemcpy(_devTetPointsMass, massVector.data(), _sysDim * sizeof(qeal), cudaMemcpyHostToDevice);

	// gpu reduced sys vector & matrix
	VectorX ZeroReducedVector(_sysReducedDim);
	ZeroReducedVector.setZero();
	CUDA_CALL(cudaMalloc((void**)&_devZeroReducedVector, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devZeroReducedVector, ZeroReducedVector.data(), ZeroReducedVector.size() * sizeof(qeal), cudaMemcpyHostToDevice));

	//
	CUDA_CALL(cudaMalloc((void**)&_devReducedX, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devReducedXn, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devReducedVelocity, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	cudaMemcpy(_devReducedVelocity, ZeroReducedVector.data(), sizeof(qeal) *_sysReducedDim, cudaMemcpyHostToDevice);
	CUDA_CALL(cudaMalloc((void**)&_devReducedVn, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devInertia, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devReducedXtilde, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devReducedDir, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devSearchReducedX, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);


	_hostReducedSolvedResidual.resize(_sysReducedDim);
	CUDA_CALL(cudaMalloc((void**)&_devReducedSolvedResidual, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	//
	CUDA_CALL(cudaMalloc((void**)&_devReducedExternalForce, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devReducedInternalForce, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devReducedRhs, _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * sizeof(qeal);


	CUDA_CALL(cudaMalloc((void**)&_devReducedMassMatrix, _sysReducedMass.size() * sizeof(qeal))); gpuSize += _sysReducedMass.size() * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(_devReducedMassMatrix, _sysReducedMass.data(), _sysReducedMass.size() * sizeof(qeal), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devReducedStiffness, _sysReducedDim * _sysReducedDim * sizeof(qeal))); gpuSize += _sysReducedDim * _sysReducedDim * sizeof(qeal);
	CUDA_CALL(cudaMalloc((void**)&_devReducedMatrix, _sysReducedMatrix.size() * sizeof(qeal))); gpuSize += _sysReducedMatrix.size() * sizeof(qeal);
}

void MIPC::MipcSimulator::initCudaReducedProjectionMemory()
{
	int bufferOffset = 0;
	CUDA_CALL(cudaMalloc((void**)&_devNonStaticFramesNum, sizeof(int))); gpuSize += sizeof(int);
	CUDA_CALL(cudaMemcpy(_devNonStaticFramesNum, &_nonStaticFramesNum, sizeof(int), cudaMemcpyHostToDevice));

	_hostTetPointXOffset.resize(models.size());
	_hostTetElementXOffset.resize(models.size());
	_hostFrameBufferOffset.resize(models.size());
	_hostFrameBufferDim.resize(models.size());
	for (int i = 0; i < models.size(); i++)
	{
		MipcModel* m = getModel(i);
		_hostTetPointXOffset[i] = m->tetPoints.offset;
		_hostTetElementXOffset[i] = m->tetElementIndices.offset * 3;
		ReducedFrame* frame = m->getReducedFrame(0);
		_hostFrameBufferOffset[i] = frame->getOffset();
		_hostFrameBufferDim[i] = m->getNonStaticFramesNum() * 12;
	}
	//
	_devTetPointProjectionXYZ.resize(models.size());
	_hostTetPointProjectionRowsXYZ.resize(models.size());
	_hostTetPointProjectionColsXYZ.resize(models.size());

	std::vector<MatrixX> hostTetElementProjection(models.size());
	_devTetElementProjectionXYZ.resize(models.size());
	_hostTetElementProjectionRowsXYZ.resize(models.size());
	_hostTetElementProjectionColsXYZ.resize(models.size());

	for (int i = 0; i < models.size(); i++)
	{
		MipcModel* m = getModel(i);
		int sub_rows_size;
		int sub_cols_size;
		MatrixX reducedProjection = m->getReducedProjection();
		MatrixX sysProjectMatrixX;

		sub_rows_size = reducedProjection.rows() / 3;
		sub_cols_size = reducedProjection.cols() / 3;
		sysProjectMatrixX.resize(sub_rows_size, sub_cols_size);
		sysProjectMatrixX.setZero();

		for (uint32_t i = 0; i < sub_rows_size; i++)
		{
			for (uint32_t j = 0; j < sub_cols_size; j++)
			{
				qeal val_x = reducedProjection(3 * i, 3 * j);
				sysProjectMatrixX(i, j) = val_x;
			}
		}

		CUDA_CALL(cudaMalloc((void**)&_devTetPointProjectionXYZ[i], sysProjectMatrixX.size() * sizeof(qeal))); gpuSize += sysProjectMatrixX.size() * sizeof(qeal);
		CUDA_CALL(cudaMemcpy(_devTetPointProjectionXYZ[i], sysProjectMatrixX.data(), sysProjectMatrixX.size() * sizeof(qeal), cudaMemcpyHostToDevice));

		_hostTetPointProjectionRowsXYZ[i] = sub_rows_size;
		_hostTetPointProjectionColsXYZ[i] = sub_cols_size;
	}

	bufferOffset = 0;
	int bufferCount = 0;
	std::vector<int>_hostTetElementFrameProjectionBufferNum;
	std::vector<int>_hostTetElementFrameProjectionBufferOffset;
	_hostTetElementFrameProjectionBufferNum.resize(totalTetElementNum, 0);
	_hostTetElementFrameProjectionBufferOffset.resize(totalTetElementNum, 0);
	std::vector<qeal> hostTetElementFrameWeightList;

	for (int i = 0; i < models.size(); i++)
	{
		MipcModel* m = getModel(i);
		for (int eleId = 0; eleId < m->tetElementNum; eleId++)
		{
			int geleId = m->getTetElementOverallId(eleId);
			std::set<int>::iterator it = m->_tetElementShareFramesList[eleId].begin();

			Vector4i ele = m->getTetElement(eleId);
			Vector3 v0 = m->getTetPoint(ele[0]);
			Vector3 v1 = m->getTetPoint(ele[1]);
			Vector3 v2 = m->getTetPoint(ele[2]);
			Vector3 v3 = m->getTetPoint(ele[3]);

			int num = 0;
			for (; it != m->_tetElementShareFramesList[eleId].end(); ++it)
			{
				ReducedFrame* frame = m->getReducedFrame(*it);
				int gFrameId = frame->getFrameId();

				qeal w0, w1, w2, w3;
				w0 = m->_harmonicWeight.data()[(*it) * m->_harmonicWeight.rows() + ele[0]];
				w1 = m->_harmonicWeight.data()[(*it) * m->_harmonicWeight.rows() + ele[1]];
				w2 = m->_harmonicWeight.data()[(*it) * m->_harmonicWeight.rows() + ele[2]];
				w3 = m->_harmonicWeight.data()[(*it) * m->_harmonicWeight.rows() + ele[3]];
				hostTetElementFrameWeightList.push_back(w0);
				hostTetElementFrameWeightList.push_back(w1);
				hostTetElementFrameWeightList.push_back(w2);
				hostTetElementFrameWeightList.push_back(w3);

				bufferCount += 16;

				num++;
			}
			_hostTetElementFrameProjectionBufferNum[geleId] = num;
			_hostTetElementFrameProjectionBufferOffset[geleId] = bufferOffset;


			bufferOffset += num;
		}
	}

	CUDA_CALL(cudaMalloc((void**)&_devTetElementFrameProjectionNum, _hostTetElementFrameProjectionBufferNum.size() * sizeof(int))); gpuSize += _hostTetElementFrameProjectionBufferNum.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTetElementFrameProjectionNum, _hostTetElementFrameProjectionBufferNum.data(), _hostTetElementFrameProjectionBufferNum.size() * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTetElementFrameProjectionOffset, _hostTetElementFrameProjectionBufferOffset.size() * sizeof(int))); gpuSize += _hostTetElementFrameProjectionBufferOffset.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTetElementFrameProjectionOffset, _hostTetElementFrameProjectionBufferOffset.data(), _hostTetElementFrameProjectionBufferOffset.size() * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTetElementFrameProjectionBuffer, bufferCount * sizeof(qeal))); gpuSize += bufferCount * sizeof(qeal);

	qeal* devTetPoints;
	CUDA_CALL(cudaMalloc((void**)&devTetPoints, 3 * totalTetPointsNum * sizeof(qeal))); gpuSize += 3 * totalTetPointsNum * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(devTetPoints, tetPointsBuffer.buffer.data(), 3 * totalTetPointsNum * sizeof(qeal), cudaMemcpyHostToDevice));

	qeal* devTetElementFrameWeightList;
	CUDA_CALL(cudaMalloc((void**)&devTetElementFrameWeightList, hostTetElementFrameWeightList.size() * sizeof(qeal))); gpuSize += hostTetElementFrameWeightList.size() * sizeof(qeal);
	CUDA_CALL(cudaMemcpy(devTetElementFrameWeightList, hostTetElementFrameWeightList.data(), hostTetElementFrameWeightList.size() * sizeof(qeal), cudaMemcpyHostToDevice));

	fillTetElementFrameProjectionBuffer
	(
		totalTetElementNum,
		_devTotalTetElementNum,
		_devTetElementIndices,
		devTetPoints,
		devTetElementFrameWeightList,
		_devTetElementFrameProjectionNum,
		_devTetElementFrameProjectionOffset,
		_devTetElementFrameProjectionBuffer
	);
	cudaFree(devTetPoints);
	cudaFree(devTetElementFrameWeightList);

	//
	_hostPojectionStiffnessNum = 0;
	std::vector<int> hostPojectionStiffnessList;
	std::vector<std::vector<int>> tetElementSharedFrameList(totalTetElementNum);
	std::vector<int> hostTetElementSharedFrameList;
	std::vector<int> hostTetElementSharedFrameNum(totalTetElementNum);
	std::vector<int> hostTetElementSharedFrameOffset(totalTetElementNum);

	for (int mid = 0; mid < models.size(); mid++)
	{
		MipcModel* m = getModel(mid);

		for (int eleId = 0; eleId < m->tetElementNum; eleId++)
		{
			int geleId = m->getTetElementOverallId(eleId);
			std::set<int>::iterator it = m->_tetElementShareFramesList[eleId].begin();
			for (; it != m->_tetElementShareFramesList[eleId].end(); ++it)
			{
				int frameIndex = *it;
				int frameId = m->getReducedFrame(*it)->getFrameId();
				tetElementSharedFrameList[geleId].push_back(frameId);
			}

			hostTetElementSharedFrameNum[geleId] = tetElementSharedFrameList[geleId].size();
			hostTetElementSharedFrameOffset[geleId] = hostTetElementSharedFrameList.size();

			for (int i = 0; i < tetElementSharedFrameList[geleId].size(); i++)
			{
				int iFrameId = tetElementSharedFrameList[geleId][i];

				hostTetElementSharedFrameList.push_back(iFrameId);

				hostPojectionStiffnessList.push_back(geleId);
				hostPojectionStiffnessList.push_back(i);
				hostPojectionStiffnessList.push_back(i);

				for (int j = i + 1; j < tetElementSharedFrameList[geleId].size(); j++)
				{
					int jFrameId = tetElementSharedFrameList[geleId][j];
					if (iFrameId < jFrameId)
					{
						hostPojectionStiffnessList.push_back(geleId);
						hostPojectionStiffnessList.push_back(i);
						hostPojectionStiffnessList.push_back(j);
					}
					else
					{
						hostPojectionStiffnessList.push_back(geleId);
						hostPojectionStiffnessList.push_back(j);
						hostPojectionStiffnessList.push_back(i);
					}
				}
			}
		}
	}
	_hostPojectionStiffnessNum = hostPojectionStiffnessList.size() / 3;

	CUDA_CALL(cudaMalloc((void**)&_devPojectionStiffnessNum, sizeof(int))); gpuSize += sizeof(int);
	CUDA_CALL(cudaMemcpy(_devPojectionStiffnessNum, &_hostPojectionStiffnessNum, sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devPojectionStiffnessList, hostPojectionStiffnessList.size() * sizeof(int))); gpuSize += hostPojectionStiffnessList.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devPojectionStiffnessList, hostPojectionStiffnessList.data(), hostPojectionStiffnessList.size() * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTetElementSharedFrameList, hostTetElementSharedFrameList.size() * sizeof(int))); gpuSize += hostTetElementSharedFrameList.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTetElementSharedFrameList, hostTetElementSharedFrameList.data(), hostTetElementSharedFrameList.size() * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTetElementSharedFrameNum, hostTetElementSharedFrameNum.size() * sizeof(int))); gpuSize += hostTetElementSharedFrameNum.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTetElementSharedFrameNum, hostTetElementSharedFrameNum.data(), hostTetElementSharedFrameNum.size() * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devTetElementSharedFrameOffset, hostTetElementSharedFrameOffset.size() * sizeof(int))); gpuSize += hostTetElementSharedFrameOffset.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devTetElementSharedFrameOffset, hostTetElementSharedFrameOffset.data(), hostTetElementSharedFrameOffset.size() * sizeof(int), cudaMemcpyHostToDevice));

	std::vector<std::vector<int>> stiffnessBlockSharedTetElementList(_nonStaticFramesNum * _nonStaticFramesNum);

	for (uint32_t i = 0; i < _hostPojectionStiffnessNum; i++)
	{
		int eid = hostPojectionStiffnessList[3 * i]; // geid
		int ci = hostPojectionStiffnessList[3 * i + 1];
		int cj = hostPojectionStiffnessList[3 * i + 2];

		int iFrameId = tetElementSharedFrameList[eid][ci];
		int jFrameId = tetElementSharedFrameList[eid][cj];

		stiffnessBlockSharedTetElementList[jFrameId * _nonStaticFramesNum + iFrameId].push_back(i);
	}
	std::vector<int> hostAssembleBlockIndex;
	std::vector<int> hoststiffnessBlockSharedTetElementList;
	std::vector<int> hostStiffnessBlockSharedTetElementNum;
	std::vector<int> hostStiffnessBlockSharedTetElementOffset;

	for (int i = 0; i < stiffnessBlockSharedTetElementList.size(); i++)
	{
		if (stiffnessBlockSharedTetElementList[i].size() == 0)
			continue;

		int iFrameId = i % _nonStaticFramesNum;
		int jFrameId = (i - iFrameId) / _nonStaticFramesNum;
		hostAssembleBlockIndex.push_back(iFrameId);
		hostAssembleBlockIndex.push_back(jFrameId);

		int num = stiffnessBlockSharedTetElementList[i].size();
		hostStiffnessBlockSharedTetElementNum.push_back(num);
		hostStiffnessBlockSharedTetElementOffset.push_back(hoststiffnessBlockSharedTetElementList.size());

		for (int j = 0; j < num; j++)
		{
			int idx = stiffnessBlockSharedTetElementList[i][j];
			hoststiffnessBlockSharedTetElementList.push_back(idx);
		}
	}

	_hostAssembleBlockNum = hostAssembleBlockIndex.size() / 2;
	CUDA_CALL(cudaMalloc((void**)&_devAssembleBlockNum, sizeof(int))); gpuSize += sizeof(int);
	CUDA_CALL(cudaMemcpy(_devAssembleBlockNum, &_hostAssembleBlockNum, sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devAssembleBlockIndex, hostAssembleBlockIndex.size() * sizeof(int))); gpuSize += hostAssembleBlockIndex.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devAssembleBlockIndex, hostAssembleBlockIndex.data(), hostAssembleBlockIndex.size() * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devStiffnessBlockSharedTetElementList, hoststiffnessBlockSharedTetElementList.size() * sizeof(int))); gpuSize += hoststiffnessBlockSharedTetElementList.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devStiffnessBlockSharedTetElementList, hoststiffnessBlockSharedTetElementList.data(), hoststiffnessBlockSharedTetElementList.size() * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devStiffnessBlockSharedTetElementNum, hostStiffnessBlockSharedTetElementNum.size() * sizeof(int))); gpuSize += hostStiffnessBlockSharedTetElementNum.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devStiffnessBlockSharedTetElementNum, hostStiffnessBlockSharedTetElementNum.data(), hostStiffnessBlockSharedTetElementNum.size() * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devStiffnessBlockSharedTetElementOffset, hostStiffnessBlockSharedTetElementOffset.size() * sizeof(int))); gpuSize += hostStiffnessBlockSharedTetElementOffset.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devStiffnessBlockSharedTetElementOffset, hostStiffnessBlockSharedTetElementOffset.data(), hostStiffnessBlockSharedTetElementOffset.size() * sizeof(int), cudaMemcpyHostToDevice));
}

void MIPC::MipcSimulator::initCudaSparseSysMemory()
{
	// implement for sparse system
}

void MIPC::MipcSimulator::initCudaCollisionMemory()
{
	_hostCollisionEventNum = _overallCollisionEvents.size();
	CUDA_CALL(cudaMalloc((void**)&_devCollisionEventNum, sizeof(int))); gpuSize += sizeof(int);
	CUDA_CALL(cudaMemcpy(_devCollisionEventNum, &_hostCollisionEventNum, sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devCollisionEventList, _hostCollisionEventList.size() * sizeof(int))); gpuSize += _hostCollisionEventList.size() * sizeof(int);
	CUDA_CALL(cudaMemcpy(_devCollisionEventList, _hostCollisionEventList.data(), _hostCollisionEventList.size() * sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMalloc((void**)&_devCCD, _hostCollisionEventNum * sizeof(qeal))); gpuSize += _hostCollisionEventNum * sizeof(qeal);
}

void MIPC::MipcSimulator::genOverallCollisionEvents()
{
	// create fake reduced frame for static model;
	staticReducedFrameList.resize(staticModelPool.totalMedialPoinsNum, nullptr);
	for (int i = 0; i < staticModelPool.totalMedialPoinsNum; i++)
		staticReducedFrameList[i] = new ReducedFrame(staticModelPool.medialPointsBuffer.buffer.data() + (3 * i), -i - 1);

	//create collsion medial sphere
	_collisionMedialSpheres.resize(totalMedialPoinsNum);
	for (int i = 0; i < totalMedialPoinsNum; i++)
	{
		MedialSphereFrame* frame = _reducedFrameList[i];
		qeal* r = medialRadiusBuffer.buffer.data() + i;
		_collisionMedialSpheres[i] = new CollideMedialSphere(frame, r);
	}
	_collisionStaticMedialSpheres.resize(staticModelPool.totalMedialPoinsNum);
	for (int i = 0; i < staticModelPool.totalMedialPoinsNum; i++)
	{
		MedialSphereFrame* frame = staticReducedFrameList[i];
		qeal* r = staticModelPool.medialRadiusBuffer.buffer.data() + i;
		_collisionStaticMedialSpheres[i] = new CollideMedialSphere(frame, r);
	}

	for (int i = 0; i < models.size(); i++)
	{
		BaseMedialMesh* mi = models[i]->getMedialMeshHandle()->getMesh();
		genOverallInterCollisionEvents(i, mi, _overallCollisionEvents);

		for (int j = i + 1; j < models.size(); j++)
		{
			BaseMedialMesh* mj = models[j]->getMedialMeshHandle()->getMesh();
			genOverallIntraCollisionEvents(i, mi, j, mj, _overallCollisionEvents);
		}

		for (int j = 0; j < staticModels.size(); j++)
		{
			BaseMedialMesh* mj = staticModels[j]->getMedialMeshHandle()->getMesh();
			genOverallStaticCollisionEvents(i, mi, j, mj, _overallCollisionEvents);
		}

	}

}

void MIPC::MipcSimulator::genOverallInterCollisionEvents(int mid, BaseMedialMesh* m, std::vector<MipcConstraint*>& collisionEventsList)
{
	if (m == nullptr)
		return;
	for (int i = 0; i < m->medialPointsNum; i++)
	{
		int gnid = m->getMedialPointOverallId(i);
		CollideMedialSphere* s = _collisionMedialSpheres[gnid];
		for (int j = 0; j < m->medialSlabsNum; j++)
		{
			Vector3i slab = m->getMedialSlab(j);
			if (slab[2] < 0)
				continue;
			if (i == slab.data()[0] || i == slab.data()[1] || i == slab.data()[2])
				continue;
			bool isNeighbor = false;
			for (int ns = 0; ns < m->medialPointsNeighborList[i].span; ns++)
			{
				int ns_id = m->medialPointsNeighborList[i].buffer[ns];
				if (ns_id == slab.data()[0] || ns_id == slab.data()[1] || ns_id == slab.data()[2])
				{
					isNeighbor = true;
					break;
				}
			}
			if (isNeighbor)
				continue;

			int n_gnid0 = m->getMedialPointOverallId(slab.data()[0]);
			int n_gnid1 = m->getMedialPointOverallId(slab.data()[1]);
			int n_gnid2 = m->getMedialPointOverallId(slab.data()[2]);
			CollideMedialSphere* ns0 = _collisionMedialSpheres[n_gnid0];
			CollideMedialSphere* ns1 = _collisionMedialSpheres[n_gnid1];
			CollideMedialSphere* ns2 = _collisionMedialSpheres[n_gnid2];

			//
			CollisionType type = DefromableWithDefromable;
			MipcSlabSphereConstraint* event = new MipcSlabSphereConstraint(collisionEventsList.size(), ns0, ns1, ns2, s, _dHat, _mu, _ev * _timeStep, type);
			if (event->isActive() || event->distance < 0)
			{
				free(event);
				continue;
			}

			collisionEventsList.push_back(event);
			_hostCollisionEventList.push_back(COLLISION_SS);
			_hostCollisionEventList.push_back(n_gnid0);
			_hostCollisionEventList.push_back(n_gnid1);
			_hostCollisionEventList.push_back(n_gnid2);
			_hostCollisionEventList.push_back(gnid);
		}
	}

	for (int i = 0; i < m->edgeList.size(); i++)
	{
		Vector2i mc = m->edgeList[i];
		CollideMedialSphere* cs0 = _collisionMedialSpheres[m->getMedialPointOverallId(mc.data()[0])];
		CollideMedialSphere* cs1 = _collisionMedialSpheres[m->getMedialPointOverallId(mc.data()[1])];
		for (int j = i + 1; j < m->edgeList.size(); j++)
		{
			Vector2i n_mc = m->edgeList[j];
			if (mc.data()[0] == n_mc.data()[0] || mc.data()[1] == n_mc.data()[1] || mc.data()[0] == n_mc.data()[1] || mc.data()[1] == n_mc.data()[0]) continue;
			CollideMedialSphere* cs3 = _collisionMedialSpheres[m->getMedialPointOverallId(n_mc.data()[0])];
			CollideMedialSphere* cs4 = _collisionMedialSpheres[m->getMedialPointOverallId(n_mc.data()[1])];

			CollisionType type = DefromableWithDefromable;
			MipcConeConeConstraint* event = new MipcConeConeConstraint(collisionEventsList.size(), cs0, cs1, cs3, cs4, _dHat, _mu, _ev * _timeStep, type);
			if (event->isActive() || event->distance < 0)
			{
				free(event);
				continue;
			}

			collisionEventsList.push_back(event);
			_hostCollisionEventList.push_back(COLLISION_CC);
			_hostCollisionEventList.push_back(m->getMedialPointOverallId(mc.data()[0]));
			_hostCollisionEventList.push_back(m->getMedialPointOverallId(mc.data()[1]));
			_hostCollisionEventList.push_back(m->getMedialPointOverallId(n_mc.data()[0]));
			_hostCollisionEventList.push_back(m->getMedialPointOverallId(n_mc.data()[1]));


		}
	}
}

void MIPC::MipcSimulator::genOverallIntraCollisionEvents(int mid1, BaseMedialMesh* m1, int mid2, BaseMedialMesh* m2, std::vector<MipcConstraint*>& collisionEventsList)
{
	if (m1 == nullptr || m2 == nullptr)
		return;
	for (int i = 0; i < m1->medialPointsNum; i++)
	{
		CollideMedialSphere* s = _collisionMedialSpheres[m1->getMedialPointOverallId(i)];
		int i1 = m1->getMedialPointOverallId(i);
		for (int j = 0; j < m2->medialSlabsNum; j++)
		{
			Vector3i slab = m2->getMedialSlab(j);
			CollideMedialSphere* ns0 = _collisionMedialSpheres[m2->getMedialPointOverallId(slab.data()[0])];
			int id = m2->getMedialPointOverallId(slab.data()[0]);
			CollideMedialSphere* ns1 = _collisionMedialSpheres[m2->getMedialPointOverallId(slab.data()[1])];
			CollideMedialSphere* ns2 = _collisionMedialSpheres[m2->getMedialPointOverallId(slab.data()[2])];

			CollisionType type = DefromableWithDefromable;
			MipcSlabSphereConstraint* event = new MipcSlabSphereConstraint(collisionEventsList.size(), ns0, ns1, ns2, s, _dHat, _mu, _ev * _timeStep, type);
			if (event->isActive() || event->distance < 0)
			{
				free(event);
				continue;
			}
			collisionEventsList.push_back(event);
			_hostCollisionEventList.push_back(COLLISION_SS);
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(slab.data()[0]));
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(slab.data()[1]));
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(slab.data()[2]));
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(i));
		}
	}

	for (int i = 0; i < m2->medialPointsNum; i++)
	{
		CollideMedialSphere* s = _collisionMedialSpheres[m2->getMedialPointOverallId(i)];
		for (int j = 0; j < m1->medialSlabsNum; j++)
		{
			Vector3i slab = m1->getMedialSlab(j);
			CollideMedialSphere* ns0 = _collisionMedialSpheres[m1->getMedialPointOverallId(slab.data()[0])];
			CollideMedialSphere* ns1 = _collisionMedialSpheres[m1->getMedialPointOverallId(slab.data()[1])];
			CollideMedialSphere* ns2 = _collisionMedialSpheres[m1->getMedialPointOverallId(slab.data()[2])];

			CollisionType type = DefromableWithDefromable;
			MipcSlabSphereConstraint* event = new MipcSlabSphereConstraint(collisionEventsList.size(), ns0, ns1, ns2, s, _dHat, _mu, _ev * _timeStep, type);
			if (event->isActive() || event->distance < 0)
			{
				free(event);
				continue;
			}

			collisionEventsList.push_back(event);

			_hostCollisionEventList.push_back(COLLISION_SS);
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(slab.data()[0]));
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(slab.data()[1]));
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(slab.data()[2]));
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(i));
		}
	}

	for (int i = 0; i < m1->edgeList.size(); i++)
	{
		Vector2i mc = m1->edgeList[i];
		CollideMedialSphere* cs0 = _collisionMedialSpheres[m1->getMedialPointOverallId(mc.data()[0])];
		CollideMedialSphere* cs1 = _collisionMedialSpheres[m1->getMedialPointOverallId(mc.data()[1])];
		for (int j = 0; j < m2->edgeList.size(); j++)
		{
			Vector2i n_mc = m2->edgeList[j];
			CollideMedialSphere* cs3 = _collisionMedialSpheres[m2->getMedialPointOverallId(n_mc.data()[0])];
			CollideMedialSphere* cs4 = _collisionMedialSpheres[m2->getMedialPointOverallId(n_mc.data()[1])];

			CollisionType type = DefromableWithDefromable;
			MipcConeConeConstraint* event = new MipcConeConeConstraint(collisionEventsList.size(), cs0, cs1, cs3, cs4, _dHat, _mu, _ev * _timeStep, type);
			if (event->isActive() || event->distance < 0)
			{
				free(event);
				continue;
			}

			collisionEventsList.push_back(event);
			_hostCollisionEventList.push_back(COLLISION_CC);
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(mc.data()[0]));
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(mc.data()[1]));
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(n_mc.data()[0]));
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(n_mc.data()[1]));
		}
	}
}

void MIPC::MipcSimulator::genOverallStaticCollisionEvents(int mid1, BaseMedialMesh* m1, int mid2, BaseMedialMesh* m2, std::vector<MipcConstraint*>& collisionEventsList)
{
	if (m1 == nullptr || m2 == nullptr)
		return;

	for (int i = 0; i < m1->medialPointsNum; i++)
	{
		CollideMedialSphere* s = _collisionMedialSpheres[m1->getMedialPointOverallId(i)];
		for (int j = 0; j < m2->medialSlabsNum; j++)
		{
			Vector3i slab = m2->getMedialSlab(j);
			CollideMedialSphere* ns0 = _collisionStaticMedialSpheres[m2->getMedialPointOverallId(slab.data()[0])];
			CollideMedialSphere* ns1 = _collisionStaticMedialSpheres[m2->getMedialPointOverallId(slab.data()[1])];
			CollideMedialSphere* ns2 = _collisionStaticMedialSpheres[m2->getMedialPointOverallId(slab.data()[2])];

			CollisionType type = DefromableWithStatic;
			MipcSlabSphereConstraint* event = new MipcSlabSphereConstraint(collisionEventsList.size(), ns0, ns1, ns2, s, _dHat, _mu, _ev * _timeStep, type);

			if (event->isActive() || event->distance < 0)
			{
				free(event);
				continue;
			}

			collisionEventsList.push_back(event); 
			_hostCollisionEventList.push_back(COLLISION_STATIC_WITH_DEFORMABLE_SS);
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(slab.data()[0]));
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(slab.data()[1]));
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(slab.data()[2]));
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(i));
		}
	}

	for (int i = 0; i < m2->medialPointsNum; i++)
	{
		CollideMedialSphere* s = _collisionStaticMedialSpheres[m2->getMedialPointOverallId(i)];
		for (int j = 0; j < m1->medialSlabsNum; j++)
		{
			Vector3i slab = m1->getMedialSlab(j);
			CollideMedialSphere* ns0 = _collisionMedialSpheres[m1->getMedialPointOverallId(slab.data()[0])];
			CollideMedialSphere* ns1 = _collisionMedialSpheres[m1->getMedialPointOverallId(slab.data()[1])];
			CollideMedialSphere* ns2 = _collisionMedialSpheres[m1->getMedialPointOverallId(slab.data()[2])];

			CollisionType type = DefromableWithStatic;
			MipcSlabSphereConstraint* event = new MipcSlabSphereConstraint(collisionEventsList.size(), ns0, ns1, ns2, s, _dHat, _mu, _ev * _timeStep, type);

			if (event->isActive() || event->distance < 0)
			{
				free(event);
				continue;
			}
			collisionEventsList.push_back(event);
			_hostCollisionEventList.push_back(COLLISION_DEFORMABLE_WITH_STATIC_SS);
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(slab.data()[0]));
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(slab.data()[1]));
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(slab.data()[2]));
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(i));
		}
	}

	//edge && edge
	for (int i = 0; i < m1->edgeList.size(); i++)
	{
		Vector2i mc = m1->edgeList[i];
		CollideMedialSphere* cs0 = _collisionMedialSpheres[m1->getMedialPointOverallId(mc.data()[0])];
		CollideMedialSphere* cs1 = _collisionMedialSpheres[m1->getMedialPointOverallId(mc.data()[1])];
		for (int j = 0; j < m2->edgeList.size(); j++)
		{
			Vector2i n_mc = m2->edgeList[j];
			CollideMedialSphere* cs3 = _collisionStaticMedialSpheres[m2->getMedialPointOverallId(n_mc.data()[0])];
			CollideMedialSphere* cs4 = _collisionStaticMedialSpheres[m2->getMedialPointOverallId(n_mc.data()[1])];

			CollisionType type = DefromableWithStatic;
			MipcConeConeConstraint* event = new MipcConeConeConstraint(collisionEventsList.size(), cs0, cs1, cs3, cs4, _dHat, _mu, _ev * _timeStep, type);

			if (event->isActive() || event->distance < 0)
			{
				free(event);
				continue;
			}
			collisionEventsList.push_back(event);
			_hostCollisionEventList.push_back(COLLISION_DEFORMABLE_WITH_STATIC_CC);
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(mc.data()[0]));
			_hostCollisionEventList.push_back(m1->getMedialPointOverallId(mc.data()[1]));
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(n_mc.data()[0]));
			_hostCollisionEventList.push_back(m2->getMedialPointOverallId(n_mc.data()[1]));
		}
	}
}