#pragma once
#ifndef BASE_SIMULATOR_H
#define BASE_SIMULATOR_H

#include "Commom\tinyxml\tinyxml.h"
#include "Model\BaseModel.h"
#include "Commom\EigenMatrixIO.h"
#include <unordered_set>
#include "Simulator\Cuda\CudaHeader.cuh"
#include "Simulator\Cuda\CudaHandle.h"

enum RunPlatform
{
	CPU,
	OPEMMP,
	CUDA
};

class StaticModelPool :public BaseSurfaceMeshBufferPool, public BaseTetMeshBufferPool, public BaseMedialMeshBufferPool
{
public:
	StaticModelPool(){}
};

class BaseSimulator: public BaseSurfaceMeshBufferPool, public BaseTetMeshBufferPool, public BaseMedialMeshBufferPool
{
public:
	BaseSimulator(std::string simName = "base_simulator", RunPlatform runPlatform = RunPlatform::CPU):
		_simulatorName(simName),
		_timeStep(0.01),
		_gravity(-9.8),
		_useGravity(true),
		_runPlatform(runPlatform)
	{
		cudaFree(0);
		setCublasAndCuSparse();	
	}
	~BaseSimulator() {
		freeCublasAndCusparse();
	}

	std::string getSimulatorName() { return _simulatorName; }

	bool readSimulatorFromConfigFile(const std::string filename, TiXmlElement* item);
	virtual void readExtraAttributeFromConfigFile(TiXmlElement* item);

	virtual bool addModelFromFile(std::vector<std::string> files, BaseModel* m);
	virtual bool addModelFromConfigFile(const std::string filename, TiXmlElement* item);
	virtual bool addModelFromConfigFile(const std::string filename, TiXmlElement* item, BaseModel* model);
	virtual bool getModelFilename(const std::string filename, TiXmlElement* item, std::vector<std::string>& files);

	virtual bool addStaticModelFromFile(std::vector<std::string> files, BaseModel* m);
	virtual bool addStaticModelFromConfigFile(const std::string filename, TiXmlElement* item);
	virtual bool addStaticModelFromConfigFile(const std::string filename, TiXmlElement* item, BaseModel* model);

	virtual void saveFile();
	virtual void saveSimulator();
	virtual void initialization() {};
	virtual void render(QOpenGLShaderProgram* program, QOpenGLFunctions* f, bool drawEdge = false);
	virtual void renderExtraElementOnCPU(QOpenGLFunctions * f);
	virtual void animate(int frame = 0) { run(frame); postRun(); }
	virtual void run(int frame = 0) {}
	virtual void postRun() {};
	virtual void reset(){}

	virtual void debug(){}

	int getModelsNum() {return models.size();}

	virtual void handleSurfacePointsSelectedEvent(){}
	virtual void handleTetPointsSelectedEvent(){}

	virtual void handleMouseForce(int nid, qeal& x, qeal& y, qeal& z) { std::cout << " handleMouseForce " << std::endl; }

	virtual void getSurfacePoint(const int id, qeal& x, qeal& y, qeal& z) { x = pointsBuffer.buffer[3 * id];  y = pointsBuffer.buffer[3 * id + 1]; z = pointsBuffer.buffer[3 * id + 2];}
	virtual void getTetPoint(const int id, qeal& x, qeal& y, qeal& z) { x = tetPointsBuffer.buffer[3 * id];  y = tetPointsBuffer.buffer[3 * id + 1]; z = tetPointsBuffer.buffer[3 * id + 2]; }
	//


	BaseModel* getStaticModel(const int id) { return staticModels[id]; }

	std::vector<BaseModel*> models;
	std::vector<BaseModel*> staticModels;

	virtual qeal getTimeStep() { return _timeStep; }
	virtual void setTimeStep(qeal t) { _timeStep = t; }

	virtual qeal getGraviry() { return _gravity; }
	virtual void setGraviry(qeal g) { _gravity = g; /* to do: renew gravity force*/}
	virtual bool isUseGravity() { return _useGravity; }
	virtual void enableGravity(bool enable) { _useGravity = enable; }

	std::unordered_set<int> selectSurfacePoints;
	std::unordered_set<int> selectTetPoints;

	//
	virtual int getSysDimension() { return _sysDim; }
	virtual void handleTetPointsConstraint();
	virtual SparseMatrix computeMassMatrix();
	virtual SparseMatrix computeInverseMassMatrix();
	virtual void alignAllMesh(qeal* tetPointsPtr);
protected:
	std::string _simulatorName;
	qeal _timeStep;
	qeal _gravity;
	bool _useGravity;

	RunPlatform _runPlatform;

	int _sysDim;

	//
	StaticModelPool staticModelPool;
	

};


#endif
