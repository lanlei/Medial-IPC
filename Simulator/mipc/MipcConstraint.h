#pragma once
#ifndef MIPC_Constraint_H
#define MIPC_Constraint_H
#define  CHECK_MEDIAL_PRIMITIVE_VALID
#include <stdio.h>
#include <iostream>
#include <queue>
#include "Simulator\FiniteElementMethod\Reduced\ReducedFrame.h" 
#include "Commom\SPDProjectFunction.h"
#include "Simulator\CollisionDetection\CollisionDetectionMedialMesh.h"

namespace MIPC
{
	using namespace CDMM;

	typedef FiniteElementMethod::ReducedFrame MedialSphereFrame;
	typedef FiniteElementMethod::LinearReducedFrame MedialSphereLinearFrame;
	typedef FiniteElementMethod::ReducedFrameType FrameType;

	class CollideMedialSphere
	{
	public:
		CollideMedialSphere(MedialSphereFrame* c, qeal* r) : center(c), radius(r)
		{}
		MedialSphereFrame* center;
		qeal* radius;
	};

	class MipcConstraint
	{
	public:
		enum CollisionType
		{
			CC = 0,
			SS = 1
		};

		int info;
		int index;
		CollisionType collisionType;
		std::vector<CollideMedialSphere* > spheres;
		qeal distance;
		qeal dHat;
		qeal dHat2;

		//
		Vector3 sC1, sC2, sC3;
		Vector3 sV1, sV2, sV3;
		qeal sR1, sR2, sR3;
		qeal A, B, C, D, E, F;
		qeal delta; // 4AC - B^2
		qeal alpha, beta;//

		// friction
		Vector3 relU;
		qeal lagLamda;
		qeal mu;
		Eigen::Matrix<qeal, 3, 2> lagBasis;
		Vector3 cloestPoints[2];
		qeal epsvh;
		qeal lagAlpha, lagBbeta;//

		MipcConstraint(int id, CollideMedialSphere* s0, CollideMedialSphere* s1, CollideMedialSphere* s2, CollideMedialSphere* s3, qeal disHat = 1.0 / 1000.0, qeal fricMu = 0.0, qeal fricEpsvh = 1e-4, int debug_info = 0):
			index(id), dHat(disHat), mu(fricMu), epsvh(fricEpsvh), info(debug_info)
		{
			spheres.push_back(s0);
			spheres.push_back(s1);
			spheres.push_back(s2);
			spheres.push_back(s3);
			dHat2 = disHat * disHat;
			relU.setZero();
			lagLamda = 0;
			lagBasis.setZero();
		}

		CollisionType getCollisionType() { return collisionType; }
		virtual qeal getDistanceHat() { return dHat; }
		virtual void setDistanceHat(qeal disHat) {dHat = disHat; dHat2= disHat * disHat;}
		virtual qeal getDistance() { return distance; }
		virtual qeal getEnergy(qeal kappa)
		{
			if (distance > dHat2)
				return 0.0;
			qeal energy = -kappa * (distance - dHat2) * (distance - dHat2) * log(distance / dHat2);
			return energy;
		}
		virtual qeal frictionEnergy() = 0;
		virtual qeal getBarrierGradient()
		{
			if (distance > dHat2)
				return 0.0;
			return (dHat2 - distance) * (2 * log(distance / dHat2) - (dHat2 / distance) + 1.0);

		}
		virtual qeal getBarrierHessian()
		{
			if (distance > dHat2)
				return 0.0;
			double dhat_d = dHat2 / distance;
			return (dhat_d + 2) * dhat_d - 2 * log(distance / dHat2) - 3;
		}
		bool isActive() { return distance <= dHat2 && distance > 0.0; }

		virtual qeal computeDistance() = 0;
		virtual void getTanBasis(Eigen::Matrix<qeal, 3, 2>& lagBasis) = 0;
		virtual void computeLagTangentBasis(const qeal kappa) = 0;

		virtual void getGradientAndHessian(qeal kappa, VectorX& gradient, MatrixX& hessian) = 0;
		virtual void getFrictionGradientAndHessian(qeal kappa, VectorX& gradient, MatrixX& hessian) = 0;

		virtual void getGradientAndHessian(qeal kappa, VectorX& gradient, std::vector<TripletX>& triplet) = 0;
		virtual void getFrictionGradientAndHessian(qeal kappa, VectorX& gradient, std::vector<TripletX>& triplet) = 0;

		void fillOverallGradient(qeal S, VectorX& dbdx, VectorX& gradient);
		void fillOverallHessian(qeal S, MatrixX& dbdx2, MatrixX& hessian);
		void fillOverallHessian(qeal S, MatrixX & dbdx2, std::vector<TripletX>& triplet);
	};

	class MipcConeConeConstraint : public MipcConstraint
	{
	public:
		enum DistanceMode
		{
			TWO_ENDPOINTS = 0,
			ALPHA_ZERO = 1,
			ALPHA_ONE = 2,
			BETA_ZERO = 3,
			BETA_ONE = 4,
			ALPHA_BETA = 5
		};
		DistanceMode distanceMode;

		MipcConeConeConstraint(int id, CollideMedialSphere* s0, CollideMedialSphere* s1, CollideMedialSphere* s2, CollideMedialSphere* s3, qeal disHat = 1.0 / 1000.0, qeal fricMu = 0.0, qeal fricEpsvh = 1e-4, int debug_info = 0) :MipcConstraint(index, s0, s1, s2, s3, disHat, fricMu, fricEpsvh, debug_info)
		{
			collisionType = CollisionType::CC;
			computeDistance();
		}
		virtual qeal frictionEnergy();
		virtual qeal computeDistance();
		virtual void getTanBasis(Eigen::Matrix<qeal, 3, 2>& lagBasis);
		virtual void computeLagTangentBasis(const qeal kappa);

		virtual void getGradientAndHessian(qeal kappa, VectorX& gradient, MatrixX& hessian);
		virtual void getFrictionGradientAndHessian(qeal kappa, VectorX& gradient, MatrixX& hessian);

		virtual void getGradientAndHessian(qeal kappa, VectorX& gradient, std::vector<TripletX>& triplet);
		virtual void getFrictionGradientAndHessian(qeal kappa, VectorX& gradient, std::vector<TripletX>& triplet);

		inline void diff_F_x(VectorX& diff);
		inline void endPointsHessina(MatrixX& hessina);
		inline void alphaIsZeroHessina(MatrixX& hessina);
		inline void alphaIsOneHessina(MatrixX& hessina);
		inline void betaIsZeroHessina(MatrixX& hessina);
		inline void betaIsOneHessina(MatrixX& hessina);
		inline void alphaBetaHessina(MatrixX& hessina);

	};

	class MipcSlabSphereConstraint : public MipcConstraint
	{
	public:
		enum DistanceMode
		{
			TWO_ENDPOINTS = 0,
			ALPHA_ZERO = 1,
			BETA_ZERO = 2,
			ALPHA_BETA_ONE = 3,
			ALPHA_BETA = 4
		};
		DistanceMode distanceMode;

		MipcSlabSphereConstraint(int id, CollideMedialSphere* s0, CollideMedialSphere* s1, CollideMedialSphere* s2, CollideMedialSphere* s3, qeal disHat = 1.0 / 1000.0, qeal fricMu = 0.0, qeal fricEpsvh = 1e-4, int debug_info = 0) :MipcConstraint(index, s0, s1, s2, s3, disHat, fricMu, fricEpsvh, debug_info)
		{
			collisionType = CollisionType::SS;
			computeDistance();
		}

		virtual qeal frictionEnergy();

		virtual qeal computeDistance();
		virtual void getTanBasis(Eigen::Matrix<qeal, 3, 2>& lagBasis);
		virtual void computeLagTangentBasis(const qeal kappa);

		virtual void getGradientAndHessian(qeal kappa, VectorX& gradient, MatrixX& hessian);
		virtual void getFrictionGradientAndHessian(qeal kappa, VectorX& gradient, MatrixX& hessian);
		
		virtual void getGradientAndHessian(qeal kappa, VectorX& gradient, std::vector<TripletX>& triplet);
		virtual void getFrictionGradientAndHessian(qeal kappa, VectorX& gradient, std::vector<TripletX>& triplet);

		inline void diff_F_x(VectorX& diff);
		inline void endPointsHessina(MatrixX& hessina);
		inline void alphaIsZeroHessina(MatrixX& hessina);
		inline void betaIsZeroHessina(MatrixX& hessina);
		inline void alphaBetaPlusOneHessina(MatrixX& hessina);
		inline void alphaBetaHessina(MatrixX& hessina);
	};


}



















#endif