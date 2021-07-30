#pragma once
#ifndef Mipc_MODEL_H
#define Mipc_MODEL_H
#include "Simulator\FiniteElementMethod\FemModel.h"
#include "Simulator\FiniteElementMethod\Reduced\ReducedFrame.h"


namespace MIPC
{
	using namespace FiniteElementMethod;
	class MipcModel : public FemModel
	{
	public:
		typedef enum {
			DENSE = 0,
			SPARSE = 1
		} ReducedType;
		MipcModel(ReducedType type = DENSE) :_type(type), _linearFramesNum(0), _quadFramesNum(0), _reducedDim(0), setWeightLocality(false), weightLocalityThreshold(0) {}

		ReducedType getReducedType() { return _type; }
		int getNonStaticFramesNum() { return _matFrameId.size(); }
		int getStaticFramesNum() { return _staticFramesNum; }
		int getLinearFramesNum() { return _linearFramesNum; }
		int getQuadraticFramesNum() { return _quadFramesNum; }
		int getTranslationFramesNum() { return _translaFramesNum; }
		int getReducedDim() { return _reducedDim; }

		int getNonStaticFrameMedialId(int id) { return _matFrameId[id]; }
		ReducedFrame* getReducedFrame(int id) { return _reducedFrames[id]; }
		ReducedFrameType getReducedFrameType(int id) { return _reducedTypes[id]; }

		virtual bool readFrameMatList(std::string filename);
		virtual bool writeFrameMatList(std::string filename);

		void createReducedFrame(int& frameOffset, int& bufferOffset, qeal* X, qeal* preX, qeal* Xtilde, qeal* Vel, qeal* preVel, qeal* Acc, qeal* preAcc, std::vector<ReducedFrame*>& frameList);

		virtual void alignAllMesh(qeal* tetPointsPtr);
		virtual void computeElementSetFromWeightMatrix(MatrixX& weight);

		virtual void computeLaplacianMatrix(std::vector<TripletX>& matValue);
		virtual MatrixX computeHarmonicWeight();

		virtual void getGlobalSparseProjectionMatrixTriplet(int rOffset, int cOffset, std::vector<TripletX>& triplet);
		MatrixX getReducedProjection();


		friend class MipcSimulator;
		bool setWeightLocality;
		qeal weightLocalityThreshold;
	
	//protected:
		ReducedType _type;
		int _staticFramesNum, _linearFramesNum, _quadFramesNum, _translaFramesNum;
		std::vector<int> _matFrameId;
		std::vector<int> _matFrameInverseId;
		std::vector<ReducedFrame*> _reducedFrames;
		std::vector<ReducedFrameType> _reducedTypes;
		int _reducedDim;
		std::vector<std::set<int>> _tetPointShareFramesList;
		std::vector<std::set<int>> _tetElementShareFramesList;
		std::vector<std::set<int>> _frameShareTetPointsList;
		std::vector<std::set<int>> _frameShareTetElementList;
		MatrixX _reducedProjection;

		MatrixX _harmonicWeight;
		SparseMatrix _reducedSpProjection;
	};


}

#endif