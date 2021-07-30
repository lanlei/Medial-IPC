#ifndef Mipc_SIMULATOR_H
#define Mipc_SIMULATOR_H

#include "Simulator\FiniteElementMethod\FemSimulator.h"
#include "MpsCCD.cuh"
#include "GpuFunc.cuh"
#include "MipcModel.h"
#include "MipcConstraint.h"


namespace MIPC
{	
	class MipcSimulator : public FemSimulator
	{
		enum SysMatType
		{
			DENSE = 0,
			SPARSE = 1
		};

		enum CollisionType
		{
			DefromableWithStatic = 0,
			DefromableWithDefromable = 1,
			FrictionCase = 2
		};

	public:
		MipcSimulator(std::string simName, RunPlatform runPlatform = RunPlatform::CPU) :FemSimulator(simName, runPlatform)
		{
			solverConfig = new ImplicitNewMarkSolverConfig(_timeStep);
			_staticFramesNum = 0;
			_linearFramesNum = 0;
			_quadraticFramesNum = 0;
			_translationFramesNum = 0;
			_nonStaticFramesNum = 0;

			_kappa = 500.0;
			_dHat = 1 / 1000.0;
			_mu = 0.0;
			_ev = 1e-2;
			_sysMatType = DENSE;
		}
		virtual bool addModelFromConfigFile(const std::string filename, TiXmlElement* item);
		MipcModel* getModel(const int id) { return dynamic_cast<MipcModel*> (models[id]); }

		virtual void doTimeGpuDenseSystem(int frame = 0);
		virtual void doTimeGpuSparseSystem(int frame = 0);
		virtual qeal computeEnergy(qeal* devXn, qeal* devXtilde);
		virtual void computeElasticsHessianAndGradient(qeal * elasticsDerivative, qeal * elasticsHessian);
		virtual void constructConstraintSet(const qeal kappa, bool updateFriction);
		virtual void getToI(qeal& toi);

		virtual void computeSystemUsingCusolverDenseChol(qeal* devSys, qeal* devRhs, qeal* devX, int dim);
		virtual void computeSystemUsingCusolverSparseChol(qeal* devSys, int* devSysRowPtr, int* devSysColInd, int nnz, qeal* devRhs, qeal* devX, int dim);

		virtual void initialization();
		virtual void run(int frame);
		virtual void postRun();

		int getReducedDimension() { return _sysReducedDim; }
	protected:
		virtual void initForGpu();
		virtual void initCudaTetMeshMemory();
		virtual void initCudaMedialMeshMemory();
		virtual void initCudaGeneralSysMemory();
		virtual void initCudaReducedProjectionMemory();
		virtual void initCudaSparseSysMemory();
		virtual void initCudaCollisionMemory();

		virtual void genOverallCollisionEvents();
		virtual  void genOverallInterCollisionEvents(int mid, BaseMedialMesh* m, std::vector<MipcConstraint*>& collisionEventsList);
		virtual  void genOverallIntraCollisionEvents(int mid1, BaseMedialMesh* m1, int mid2, BaseMedialMesh* m2, std::vector<MipcConstraint*>& collisionEventsList);
		virtual  void genOverallStaticCollisionEvents(int mid1, BaseMedialMesh* m1, int mid2, BaseMedialMesh* m2, std::vector<MipcConstraint*>& collisionEventsList);

		bool enableFriction = false;
		int _sysReducedDim;
		std::vector<int> _sysNonStaticFrameOverallId;
		std::vector<int> _sysReducedOffsetByModel;
		MatrixX _sysReducedMatrix;
		VectorX _sysReducedRhs;
		VectorX _sysReducedX;
		VectorX _sysReducedXtilde;
		VectorX _sysReducedVelocity;
		VectorX _sysReducedAccrelation;
		VectorX _sysReducedExternalForce;
		VectorX _sysReducedInternalForce;
		VectorX _sysReducedDir;

		VectorX _sysReducedXn;
		VectorX _sysReducedVn;
		VectorX _sysReducedAccn;

		MatrixX _sysReducedMass;
		MatrixX _sysReducedStiffness;
		MatrixX _sysReducedDamping;

		MatrixX _sysReducedProjection;
		MatrixX _sysReducedProjectionT;

		SparseMatrix _sysReducedSparseProjection;
		SparseMatrix _sysReducedSparseProjectionT;

		MatrixX _sysReducedMedialProjection;
		MatrixX _sysReducedMedialProjectionT;

		Eigen::LLT<MatrixX> _llt;

		int _staticFramesNum;
		int _linearFramesNum;
		int _quadraticFramesNum;
		int _translationFramesNum;
		int _nonStaticFramesNum;
		std::vector<ReducedFrame*> _reducedFrameList;

		// IPC
		qeal _kappa; qeal _dHat;
		qeal _mu, _ev;
		qeal tol;

		std::vector<ReducedFrame*> staticReducedFrameList;
		std::vector<CollideMedialSphere*> _collisionMedialSpheres;
		std::vector<CollideMedialSphere*> _collisionStaticMedialSpheres;
		std::vector<MipcConstraint*> _overallCollisionEvents;
		std::vector<MipcConstraint*> _activeCollisionEvents;
		std::vector<MipcConstraint*> _frictionCollisionEvents;

		//Gpu
		long long int gpuSize;
		SysMatType _sysMatType;
		int* _devDim;
		int* _devReducedDim;
		qeal* _devTimeStep;

		// gpu fullspace sys vector & matrix
		qeal* _devSysX;
		qeal* _devSysXn;
		qeal* _devSysXtilde;
		qeal* _devDir;
		qeal* _devSearchX;

		qeal* _devDiffX;
		qeal* _devSysVelocity;
		qeal* _devMassMulDiffX;
		qeal* _devExternalForce;
		qeal* _devInternalForce;
		qeal* _devTetPointsMass;

		// gpu reduced sys vector & matrix
		qeal* _devZeroReducedMatrix;
		qeal* _devZeroReducedVector;
		qeal* _devReducedX;
		qeal* _devReducedXn;
		qeal* _devReducedVelocity;
		qeal* _devReducedVn;
		qeal* _devInertia;
		qeal* _devReducedXtilde;
		qeal* _devReducedSolvedResidual;
		VectorX _hostReducedSolvedResidual;

		qeal* _devReducedExternalForce;
		qeal* _devReducedInternalForce;
		qeal* _devReducedRhs;

		qeal* _devReducedMassMatrix;
		qeal* _devReducedStiffness;
		qeal* _devReducedMatrix;
		qeal* _devReducedDir;
		qeal* _devSearchReducedX;

		// gpu sparse sys matrix
		// reduced sparse mass matrix
		SparseMatrix _hostReducedSparseMass;
		std::vector<int> _hostReducedSparseMassCsrRowPtr;
		int* _devReducedSparseMassCsrRowPtr;
		std::vector<int> _hostReducedSparseMassCsrColInd;
		int* _devReducedSparseMassCsrColInd;
		std::vector<qeal> _hostReducedSparseMassCsrVal;
		qeal* _devReducedSparseMassCsrVal;
		int _hostReducedSparseMassCsrNonZero;
		int* _devReducedSparseMassCsrNonZero;
		// reduced sparse stiffness matrix
		SparseMatrix _hostReducedSparseStiffness;
		std::vector<int> _hostReducedSparseStiffnessCsrRowPtr;
		int* _devReducedSparseStiffnessCsrRowPtr;
		std::vector<int> _hostReducedSparseStiffnessCsrColInd;
		int* _devReducedSparseStiffnessCsrColInd;
		std::vector<qeal> _hostReducedSparseStiffnessCsrVal;
		qeal* _devReducedSparseStiffnessCsrVal;
		int _hostReducedSparseStiffnessCsrNonZero;
		int* _devReducedSparseStiffnessCsrNonZero;
		int* _devSparseStiffnessBufferIndex;
		// reduced sparse elastic hessina matrix
		SparseMatrix _hostReducedSparseElasticsHessina;
		std::vector<int> _hostReducedSparseElasticsHessinaCsrRowPtr;
		int* _devReducedSparseElasticsHessinaCsrRowPtr;
		std::vector<int> _hostReducedSparseElasticsHessinaCsrColInd;
		int* _devReducedSparseElasticsHessinaCsrColInd;
		std::vector<qeal> _hostReducedSparseElasticsHessinaCsrVal;
		qeal* _devReducedSparseElasticsHessinaCsrVal;
		int _hostReducedSparseElasticsHessinaCsrNonZero;
		int* _devReducedSparseElasticsHessinaCsrNonZero;
		// reduced sparse collision & friction hessina matrix
		SparseMatrix _hostReducedSparseCFHessina;
		std::vector<TripletX> _hostCFHessinaTriplet;
		std::vector<int> _hostReducedSparseCFHessinaCsrRowPtr;
		int* _devReducedSparseCFHessinaCsrRowPtr;
		std::vector<int> _hostReducedSparseCFHessinaCsrColInd;
		int* _devReducedSparseCFHessinaCsrColInd;
		std::vector<qeal> _hostReducedSparseCFHessinaCsrVal;
		qeal* _devReducedSparseCFHessinaCsrVal;
		int _hostReducedSparseCFHessinaCsrNonZero;
		int* _devReducedSparseCFHessinaCsrNonZero;
		// reduced sparse system hessina matrix
		SparseMatrix _hostReducedSparseSysMatrix;
		std::vector<int> _hostReducedSparseSysMatrixCsrRowPtr;
		int* _devReducedSparseSysMatrixCsrRowPtr;
		std::vector<int> _hostReducedSparseSysMatrixCsrColInd;
		int* _devReducedSparseSysMatrixCsrColInd;
		std::vector<qeal> _hostReducedSparseSysMatrixCsrVal;
		qeal* _devReducedSparseSysMatrixCsrVal;
		int _hostReducedSparseSysMatrixCsrNonZero;
		int* _devReducedSparseSysMatrixCsrNonZero;

		// stiffness & projection & assemble
		int* _devTotalTetElementNum;
		int* _devNonStaticFramesNum;
		int* _devTotalTetPointsNum;

		qeal* _devTetElementX;
		int* _devTetElementIndices;

		int* _devTetPointsSharedElementNum;
		int* _devTetPointsSharedElementOffset;
		int* _devTetPointsSharedElementList;


		qeal* _devTetElementPotentialEnergy;
		qeal* _devTetElementSizeOne;
		qeal* _devTetElementForce;
		qeal* _devTetElementStiffness;
		qeal* _devTetElementDm;
		qeal* _devTetElementInvDm;
		qeal* _devTetElementdFdu;
		qeal* _devTetElementdPdF;
		qeal *_devTetElementdFPK;
		qeal* _devTetElementAttri;
		qeal* _devTetElementVol;

		std::vector<int> _hostTetPointXOffset;
		std::vector<int> _hostTetElementXOffset;
		std::vector<int> _hostFrameBufferOffset;
		std::vector<int> _hostFrameBufferDim;

		std::vector<qeal*> _devTetPointProjectionXYZ;
		std::vector<int> _hostTetPointProjectionRowsXYZ;
		std::vector<int> _hostTetPointProjectionColsXYZ;

		std::vector<qeal*> _devTetElementProjectionXYZ;
		std::vector<int> _hostTetElementProjectionRowsXYZ;
		std::vector<int> _hostTetElementProjectionColsXYZ;

		qeal* _devTetElementFrameProjectionBuffer;
		int* _devTetElementFrameProjectionNum;
		int* _devTetElementFrameProjectionOffset;

		int _hostPojectionStiffnessNum;
		int* _devPojectionStiffnessNum;
		int* _devPojectionStiffnessList;

		int* _devTetElementSharedFrameList;
		int* _devTetElementSharedFrameNum;
		int* _devTetElementSharedFrameOffset;

		int _hostAssembleBlockNum;
		int* _devAssembleBlockNum;
		int* _devAssembleBlockIndex;
		int* _devStiffnessBlockSharedTetElementList;
		int* _devStiffnessBlockSharedTetElementNum;
		int* _devStiffnessBlockSharedTetElementOffset;

		int _hostTetPointsStiffnessNum;
		int* _devTetPointsStiffnessNum;
		qeal* _devTetPointStiffness;
		int* _devTetPointStiffnessEleNum;
		int* _devTetPointStiffnessEleOffset;
		int* _devTetPointStiffnessEleList;

		std::vector<int>  tetPointUpperStiffnessList;
		std::vector<std::set<int>> relTetElementSet;
		// gpu medial mesh
		int* _devMedialPointsNum;
		qeal* _devMedialOriPointPosition;
		qeal* _devMedialPointPosition;
		qeal* _devMedialPointRadius;
		qeal* _devStaticMedialPointPosition;
		qeal* _devStaticMedialPointRadius;
		qeal* _devMedialPointMovingDir;
		// gpu ccd
		int _hostCollisionEventNum;
		int* _devCollisionEventNum;
		std::vector<int> _hostCollisionEventList;
		int* _devCollisionEventList;
		qeal* _devCCD;

	};


}





#endif