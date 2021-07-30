#include "MipcConstraint.h"

namespace MIPC
{
	void f0_SF(qeal x2, qeal epsvh, qeal& f0)
	{
		if (x2 >= epsvh * epsvh) {
			f0 = std::sqrt(x2);
		}
		else {
			f0 = x2 * (-std::sqrt(x2) / 3.0 + epsvh) / (epsvh * epsvh) + epsvh / 3.0;
		}
	}

	void f1_SF_Div_RelDXNorm(qeal x2, qeal epsvh, qeal& result)
	{
		if (x2 >= epsvh * epsvh) {
			result = 1 / std::sqrt(x2);
		}
		else {
			result = (-std::sqrt(x2) + 2.0 * epsvh) / (epsvh * epsvh);
		}
	}

	void f2_SF_Term(qeal x2, qeal epsvh, qeal& f2_term)
	{
		f2_term = -1 / (epsvh * epsvh);
	}

	void Point_Triangle_Tangent_Basis(
		const Eigen::Matrix<qeal, 3, 1>& v0,
		const Eigen::Matrix<qeal, 3, 1>& v1,
		const Eigen::Matrix<qeal, 3, 1>& v2,
		const Eigen::Matrix<qeal, 3, 1>& v3,
		Eigen::Matrix<qeal, 3, 2>& basis)
	{
		Eigen::Matrix<qeal, 3, 1> v12 = v2 - v1;
		basis.col(0) = v12.normalized();
		basis.col(1) = v12.cross(v3 - v1).cross(v12).normalized();
	}

	void Edge_Edge_Tangent_Basis(
		const Eigen::Matrix<qeal, 3, 1>& v0,
		const Eigen::Matrix<qeal, 3, 1>& v1,
		const Eigen::Matrix<qeal, 3, 1>& v2,
		const Eigen::Matrix<qeal, 3, 1>& v3,
		Eigen::Matrix<qeal, 3, 2>& basis)
	{
		Eigen::Matrix<qeal, 3, 1> v01 = v1 - v0;
		basis.col(0) = v01.normalized();
		basis.col(1) = v01.cross(v3 - v2).cross(v01).normalized();
	}

	void Point_Edge_Tangent_Basis(
		const Eigen::Matrix<qeal, 3, 1>& v0,
		const Eigen::Matrix<qeal, 3, 1>& v1,
		const Eigen::Matrix<qeal, 3, 1>& v2,
		Eigen::Matrix<qeal, 3, 2>& basis)
	{
		Eigen::Matrix<qeal, 3, 1> v12 = v2 - v1;
		basis.col(0) = v12.normalized();
		basis.col(1) = v12.cross(v0 - v1).normalized();
	}

	void Point_Point_Tangent_Basis(
		const Eigen::Matrix<qeal, 3, 1>& v0,
		const Eigen::Matrix<qeal, 3, 1>& v1,
		Eigen::Matrix<qeal, 3, 2>& basis)
	{
		Eigen::Matrix<qeal, 1, 3> v01 = (v1 - v0).transpose();
		Eigen::Matrix<qeal, 1, 3> xCross = Eigen::Matrix<qeal, 1, 3>::UnitX().cross(v01);
		Eigen::Matrix<qeal, 1, 3> yCross = Eigen::Matrix<qeal, 1, 3>::UnitY().cross(v01);
		if (xCross.squaredNorm() > yCross.squaredNorm()) {
			basis.col(0) = xCross.normalized().transpose();
			basis.col(1) = v01.cross(xCross).normalized().transpose();
		}
		else {
			basis.col(0) = yCross.normalized().transpose();
			basis.col(1) = v01.cross(yCross).normalized().transpose();
		}
	}

	void MipcConstraint::fillOverallGradient(qeal S, VectorX& dbdx, VectorX& gradient)
	{
		for (int i = 0; i < 4; i++)
		{
			FrameType frameType = spheres[i]->center->getFrameType();
			if (frameType == FrameType::STATIC)
				continue;
			int frameId = spheres[i]->center->getFrameId();
			int offset = spheres[i]->center->getOffset();

			if (frameType == FrameType::LINEAR)
			{
				qeal p[3];
				spheres[i]->center->getOriginalP(p);
				gradient[offset] += p[0] * S * dbdx.data()[3 * i];
				gradient[offset + 1] += p[0] * S * dbdx.data()[3 * i + 1];
				gradient[offset + 2] += p[0] * S * dbdx.data()[3 * i + 2];

				gradient[offset + 3] += p[1] * S * dbdx.data()[3 * i];
				gradient[offset + 4] += p[1] * S * dbdx.data()[3 * i + 1];
				gradient[offset + 5] += p[1] * S * dbdx.data()[3 * i + 2];

				gradient[offset + 6] += p[2] * S * dbdx.data()[3 * i];
				gradient[offset + 7] += p[2] * S * dbdx.data()[3 * i + 1];
				gradient[offset + 8] += p[2] * S * dbdx.data()[3 * i + 2];

				gradient[offset + 9] += S * dbdx.data()[3 * i];
				gradient[offset + 10] += S * dbdx.data()[3 * i + 1];
				gradient[offset + 11] += S * dbdx.data()[3 * i + 2];
			}
		}
	}

	void MipcConstraint::fillOverallHessian(qeal S, MatrixX& dbdx2, MatrixX& hessian)
	{
		int size = hessian.rows();
		for (int c = 0; c < 4; c++)
		{
			FrameType c_frameType = spheres[c]->center->getFrameType();
			if (c_frameType == FrameType::STATIC)
				continue;
			int c_frameId = spheres[c]->center->getFrameId();
			int c_offset = spheres[c]->center->getOffset();

			qeal c_p[3];
			spheres[c]->center->getOriginalP(c_p);
			for (int a = 0; a < 4; a++)
			{
				FrameType a_frameType = spheres[a]->center->getFrameType();
				if (a_frameType == FrameType::STATIC)
					continue;
				int a_frameId = spheres[a]->center->getFrameId();
				int a_offset = spheres[a]->center->getOffset();
				qeal a_p[3];
				spheres[a]->center->getOriginalP(a_p);

				qeal x0 = dbdx2.data()[12 * 3 * a + 3 * c + 0];
				qeal x1 = dbdx2.data()[12 * 3 * a + 3 * c + 1];
				qeal x2 = dbdx2.data()[12 * 3 * a + 3 * c + 2];

				qeal x3 = dbdx2.data()[12 * (3 * a + 1) + 3 * c + 0];
				qeal x4 = dbdx2.data()[12 * (3 * a + 1) + 3 * c + 1];
				qeal x5 = dbdx2.data()[12 * (3 * a + 1) + 3 * c + 2];

				qeal x6 = dbdx2.data()[12 * (3 * a + 2) + 3 * c + 0];
				qeal x7 = dbdx2.data()[12 * (3 * a + 2) + 3 * c + 1];
				qeal x8 = dbdx2.data()[12 * (3 * a + 2) + 3 * c + 2];

				// hard code
				if (c_frameType == FrameType::LINEAR)
				{
					if (a_frameType == FrameType::LINEAR)
					{
						// 24 x 24
						hessian.data()[(a_offset + 0) * size + c_offset + 0] += c_p[0] * x0 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 1] += c_p[0] * x1 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 2] += c_p[0] * x2 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 3] += c_p[1] * x0 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 4] += c_p[1] * x1 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 5] += c_p[1] * x2 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 6] += c_p[2] * x0 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 7] += c_p[2] * x1 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 8] += c_p[2] * x2 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 9] += 1.0 * x0 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 10] += 1.0 * x1 * a_p[0];
						hessian.data()[(a_offset + 0) * size + c_offset + 11] += 1.0 * x2 * a_p[0];

						hessian.data()[(a_offset + 1) * size + c_offset + 0] += c_p[0] * x3 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 1] += c_p[0] * x4 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 2] += c_p[0] * x5 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 3] += c_p[1] * x3 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 4] += c_p[1] * x4 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 5] += c_p[1] * x5 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 6] += c_p[2] * x3 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 7] += c_p[2] * x4 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 8] += c_p[2] * x5 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 9] += 1.0 * x3 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 10] += 1.0 * x4 * a_p[0];
						hessian.data()[(a_offset + 1) * size + c_offset + 11] += 1.0 * x5 * a_p[0];

						hessian.data()[(a_offset + 2) * size + c_offset + 0] += c_p[0] * x6 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 1] += c_p[0] * x7 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 2] += c_p[0] * x8 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 3] += c_p[1] * x6 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 4] += c_p[1] * x7 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 5] += c_p[1] * x8 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 6] += c_p[2] * x6 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 7] += c_p[2] * x7 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 8] += c_p[2] * x8 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 9] += 1.0 * x6 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 10] += 1.0 * x7 * a_p[0];
						hessian.data()[(a_offset + 2) * size + c_offset + 11] += 1.0 * x8 * a_p[0];

						//
						hessian.data()[(a_offset + 3) * size + c_offset + 0] += c_p[0] * x0 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 1] += c_p[0] * x1 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 2] += c_p[0] * x2 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 3] += c_p[1] * x0 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 4] += c_p[1] * x1 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 5] += c_p[1] * x2 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 6] += c_p[2] * x0 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 7] += c_p[2] * x1 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 8] += c_p[2] * x2 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 9] += 1.0 * x0 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 10] += 1.0 * x1 * a_p[1];
						hessian.data()[(a_offset + 3) * size + c_offset + 11] += 1.0 * x2 * a_p[1];

						hessian.data()[(a_offset + 4) * size + c_offset + 0] += c_p[0] * x3 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 1] += c_p[0] * x4 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 2] += c_p[0] * x5 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 3] += c_p[1] * x3 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 4] += c_p[1] * x4 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 5] += c_p[1] * x5 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 6] += c_p[2] * x3 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 7] += c_p[2] * x4 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 8] += c_p[2] * x5 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 9] += 1.0 * x3 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 10] += 1.0 * x4 * a_p[1];
						hessian.data()[(a_offset + 4) * size + c_offset + 11] += 1.0 * x5 * a_p[1];

						hessian.data()[(a_offset + 5) * size + c_offset + 0] += c_p[0] * x6 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 1] += c_p[0] * x7 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 2] += c_p[0] * x8 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 3] += c_p[1] * x6 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 4] += c_p[1] * x7 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 5] += c_p[1] * x8 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 6] += c_p[2] * x6 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 7] += c_p[2] * x7 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 8] += c_p[2] * x8 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 9] += 1.0 * x6 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 10] += 1.0 * x7 * a_p[1];
						hessian.data()[(a_offset + 5) * size + c_offset + 11] += 1.0 * x8 * a_p[1];

						//
						hessian.data()[(a_offset + 6) * size + c_offset + 0] += c_p[0] * x0 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 1] += c_p[0] * x1 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 2] += c_p[0] * x2 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 3] += c_p[1] * x0 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 4] += c_p[1] * x1 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 5] += c_p[1] * x2 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 6] += c_p[2] * x0 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 7] += c_p[2] * x1 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 8] += c_p[2] * x2 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 9] += 1.0 * x0 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 10] += 1.0 * x1 * a_p[2];
						hessian.data()[(a_offset + 6) * size + c_offset + 11] += 1.0 * x2 * a_p[2];

						hessian.data()[(a_offset + 7) * size + c_offset + 0] += c_p[0] * x3 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 1] += c_p[0] * x4 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 2] += c_p[0] * x5 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 3] += c_p[1] * x3 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 4] += c_p[1] * x4 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 5] += c_p[1] * x5 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 6] += c_p[2] * x3 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 7] += c_p[2] * x4 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 8] += c_p[2] * x5 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 9] += 1.0 * x3 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 10] += 1.0 * x4 * a_p[2];
						hessian.data()[(a_offset + 7) * size + c_offset + 11] += 1.0 * x5 * a_p[2];

						hessian.data()[(a_offset + 8) * size + c_offset + 0] += c_p[0] * x6 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 1] += c_p[0] * x7 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 2] += c_p[0] * x8 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 3] += c_p[1] * x6 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 4] += c_p[1] * x7 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 5] += c_p[1] * x8 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 6] += c_p[2] * x6 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 7] += c_p[2] * x7 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 8] += c_p[2] * x8 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 9] += 1.0 * x6 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 10] += 1.0 * x7 * a_p[2];
						hessian.data()[(a_offset + 8) * size + c_offset + 11] += 1.0 * x8 * a_p[2];
						//
						hessian.data()[(a_offset + 9) * size + c_offset + 0] += c_p[0] * x0 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 1] += c_p[0] * x1 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 2] += c_p[0] * x2 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 3] += c_p[1] * x0 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 4] += c_p[1] * x1 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 5] += c_p[1] * x2 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 6] += c_p[2] * x0 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 7] += c_p[2] * x1 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 8] += c_p[2] * x2 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 9] += 1.0 * x0 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 10] += 1.0 * x1 * 1.0;
						hessian.data()[(a_offset + 9) * size + c_offset + 11] += 1.0 * x2 * 1.0;

						hessian.data()[(a_offset + 10) * size + c_offset + 0] += c_p[0] * x3 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 1] += c_p[0] * x4 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 2] += c_p[0] * x5 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 3] += c_p[1] * x3 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 4] += c_p[1] * x4 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 5] += c_p[1] * x5 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 6] += c_p[2] * x3 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 7] += c_p[2] * x4 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 8] += c_p[2] * x5 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 9] += 1.0 * x3 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 10] += 1.0 * x4 * 1.0;
						hessian.data()[(a_offset + 10) * size + c_offset + 11] += 1.0 * x5 * 1.0;

						hessian.data()[(a_offset + 11) * size + c_offset + 0] += c_p[0] * x6 * 1.0;
						hessian.data()[(a_offset + 11) * size + c_offset + 1] += c_p[0] * x7 * 1.0;
						hessian.data()[(a_offset + 11) * size + c_offset + 2] += c_p[0] * x8 * 1.0;
						hessian.data()[(a_offset + 11) * size + c_offset + 3] += c_p[1] * x6 * 1.0;
						hessian.data()[(a_offset + 11) * size + c_offset + 4] += c_p[1] * x7 * 1.0;
						hessian.data()[(a_offset + 11) * size + c_offset + 5] += c_p[1] * x8 * 1.0;
						hessian.data()[(a_offset + 11) * size + c_offset + 6] += c_p[2] * x6 * 1.0;
						hessian.data()[(a_offset + 11) * size + c_offset + 7] += c_p[2] * x7 * 1.0;
						hessian.data()[(a_offset + 11) * size + c_offset + 8] += c_p[2] * x8 * 1.0;
						hessian.data()[(a_offset + 11) * size + c_offset + 9] += 1.0 * x6 * 1.0;
						hessian.data()[(a_offset + 11) * size + (c_offset)+10] += 1.0 * x7 * 1.0;
						hessian.data()[(a_offset + 11) * size + c_offset + 11] += 1.0 * x8 * 1.0;
					}
				}
			}
		}
	}

	void MipcConstraint::fillOverallHessian(qeal S, MatrixX & dbdx2, std::vector<TripletX>& triplet)
	{
		for (int c = 0; c < 4; c++)
		{
			FrameType c_frameType = spheres[c]->center->getFrameType();
			if (c_frameType == FrameType::STATIC)
				continue;
			int c_frameId = spheres[c]->center->getFrameId();
			int c_offset = spheres[c]->center->getOffset();

			qeal c_p[3];
			spheres[c]->center->getOriginalP(c_p);
			for (int a = 0; a < 4; a++)
			{
				FrameType a_frameType = spheres[a]->center->getFrameType();
				if (a_frameType == FrameType::STATIC)
					continue;
				int a_frameId = spheres[a]->center->getFrameId();
				int a_offset = spheres[a]->center->getOffset();
				qeal a_p[3];
				spheres[a]->center->getOriginalP(a_p);

				qeal x0 = dbdx2.data()[12 * 3 * a + 3 * c + 0];
				qeal x1 = dbdx2.data()[12 * 3 * a + 3 * c + 1];
				qeal x2 = dbdx2.data()[12 * 3 * a + 3 * c + 2];

				qeal x3 = dbdx2.data()[12 * (3 * a + 1) + 3 * c + 0];
				qeal x4 = dbdx2.data()[12 * (3 * a + 1) + 3 * c + 1];
				qeal x5 = dbdx2.data()[12 * (3 * a + 1) + 3 * c + 2];

				qeal x6 = dbdx2.data()[12 * (3 * a + 2) + 3 * c + 0];
				qeal x7 = dbdx2.data()[12 * (3 * a + 2) + 3 * c + 1];
				qeal x8 = dbdx2.data()[12 * (3 * a + 2) + 3 * c + 2];

				// hard code
				if (c_frameType == FrameType::LINEAR)
				{
					if (a_frameType == FrameType::LINEAR)
					{
						// 24 x 24
						triplet.push_back(TripletX(c_offset + 0, a_offset + 0, c_p[0] * x0 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 0, c_p[0] * x1 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 0, c_p[0] * x2 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 0, c_p[1] * x0 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 0, c_p[1] * x1 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 0, c_p[1] * x2 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 0, c_p[2] * x0 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 0, c_p[2] * x1 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 0, c_p[2] * x2 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 0, 1.0 * x0 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 0, 1.0 * x1 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 0, 1.0 * x2 * a_p[0]));
						//
						triplet.push_back(TripletX(c_offset + 0, a_offset + 1, c_p[0] * x3 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 1, c_p[0] * x4 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 1, c_p[0] * x5 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 1, c_p[1] * x3 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 1, c_p[1] * x4 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 1, c_p[1] * x5 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 1, c_p[2] * x3 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 1, c_p[2] * x4 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 1, c_p[2] * x5 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 1, 1.0 * x3 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 1, 1.0 * x4 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 1, 1.0 * x5 * a_p[0]));
						//
						triplet.push_back(TripletX(c_offset + 0, a_offset + 2, c_p[0] * x6 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 2, c_p[0] * x7 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 2, c_p[0] * x8 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 2, c_p[1] * x6 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 2, c_p[1] * x7 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 2, c_p[1] * x8 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 2, c_p[2] * x6 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 2, c_p[2] * x7 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 2, c_p[2] * x8 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 2, 1.0 * x6 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 2, 1.0 * x7 * a_p[0]));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 2, 1.0 * x8 * a_p[0]));

						//
						triplet.push_back(TripletX(c_offset + 0, a_offset + 3, c_p[0] * x0 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 3, c_p[0] * x1 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 3, c_p[0] * x2 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 3, c_p[1] * x0 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 3, c_p[1] * x1 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 3, c_p[1] * x2 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 3, c_p[2] * x0 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 3, c_p[2] * x1 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 3, c_p[2] * x2 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 3, 1.0 * x0 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 3, 1.0 * x1 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 3, 1.0 * x2 * a_p[1]));
						//
						triplet.push_back(TripletX(c_offset + 0, a_offset + 4, c_p[0] * x3 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 4, c_p[0] * x4 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 4, c_p[0] * x5 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 4, c_p[1] * x3 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 4, c_p[1] * x4 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 4, c_p[1] * x5 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 4, c_p[2] * x3 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 4, c_p[2] * x4 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 4, c_p[2] * x5 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 4, 1.0 * x3 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 4, 1.0 * x4 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 4, 1.0 * x5 * a_p[1]));

						//
						triplet.push_back(TripletX(c_offset + 0, a_offset + 5, c_p[0] * x6 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 5, c_p[0] * x7 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 5, c_p[0] * x8 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 5, c_p[1] * x6 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 5, c_p[1] * x7 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 5, c_p[1] * x8 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 5, c_p[2] * x6 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 5, c_p[2] * x7 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 5, c_p[2] * x8 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 5, 1.0 * x6 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 5, 1.0 * x7 * a_p[1]));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 5, 1.0 * x8 * a_p[1]));

						//
						triplet.push_back(TripletX(c_offset + 0, a_offset + 6, c_p[0] * x0 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 6, c_p[0] * x1 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 6, c_p[0] * x2 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 6, c_p[1] * x0 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 6, c_p[1] * x1 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 6, c_p[1] * x2 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 6, c_p[2] * x0 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 6, c_p[2] * x1 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 6, c_p[2] * x2 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 6, 1.0 * x0 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 6, 1.0 * x1 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 6, 1.0 * x2 * a_p[2]));

						//
						triplet.push_back(TripletX(c_offset + 0, a_offset + 7, c_p[0] * x3 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 7, c_p[0] * x4 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 7, c_p[0] * x5 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 7, c_p[1] * x3 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 7, c_p[1] * x4 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 7, c_p[1] * x5 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 7, c_p[2] * x3 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 7, c_p[2] * x4 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 7, c_p[2] * x5 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 7, 1.0 * x3 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 7, 1.0 * x4 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 7, 1.0 * x5 * a_p[2]));

						//
						triplet.push_back(TripletX(c_offset + 0, a_offset + 8, c_p[0] * x6 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 8, c_p[0] * x7 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 8, c_p[0] * x8 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 8, c_p[1] * x6 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 8, c_p[1] * x7 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 8, c_p[1] * x8 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 8, c_p[2] * x6 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 8, c_p[2] * x7 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 8, c_p[2] * x8 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 8, 1.0 * x6 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 8, 1.0 * x7 * a_p[2]));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 8, 1.0 * x8 * a_p[2]));

						//
						triplet.push_back(TripletX(c_offset + 0, a_offset + 9, c_p[0] * x0 * 1.0));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 9, c_p[0] * x1 * 1.0));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 9, c_p[0] * x2 * 1.0));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 9, c_p[1] * x0 * 1.0));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 9, c_p[1] * x1 *1.0));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 9, c_p[1] * x2 * 1.0));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 9, c_p[2] * x0 * 1.0));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 9, c_p[2] * x1 * 1.0));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 9, c_p[2] * x2 * 1.0));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 9, 1.0 * x0 * 1.0));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 9, 1.0 * x1 * 1.0));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 9, 1.0 * x2 * 1.0));

						triplet.push_back(TripletX(c_offset + 0, a_offset + 10, c_p[0] * x3 * 1.0));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 10, c_p[0] * x4 * 1.0));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 10, c_p[0] * x5 * 1.0));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 10, c_p[1] * x3 * 1.0));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 10, c_p[1] * x4 *1.0));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 10, c_p[1] * x5 * 1.0));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 10, c_p[2] * x3 * 1.0));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 10, c_p[2] * x4 * 1.0));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 10, c_p[2] * x5 * 1.0));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 10, 1.0 * x3 * 1.0));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 10, 1.0 * x4 * 1.0));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 10, 1.0 * x5 * 1.0));

						triplet.push_back(TripletX(c_offset + 0, a_offset + 11, c_p[0] * x6 * 1.0));
						triplet.push_back(TripletX(c_offset + 1, a_offset + 11, c_p[0] * x7 * 1.0));
						triplet.push_back(TripletX(c_offset + 2, a_offset + 11, c_p[0] * x8 * 1.0));
						triplet.push_back(TripletX(c_offset + 3, a_offset + 11, c_p[1] * x6 * 1.0));
						triplet.push_back(TripletX(c_offset + 4, a_offset + 11, c_p[1] * x7 *1.0));
						triplet.push_back(TripletX(c_offset + 5, a_offset + 11, c_p[1] * x8 * 1.0));
						triplet.push_back(TripletX(c_offset + 6, a_offset + 11, c_p[2] * x6 * 1.0));
						triplet.push_back(TripletX(c_offset + 7, a_offset + 11, c_p[2] * x7 * 1.0));
						triplet.push_back(TripletX(c_offset + 8, a_offset + 11, c_p[2] * x8 * 1.0));
						triplet.push_back(TripletX(c_offset + 9, a_offset + 11, 1.0 * x6 * 1.0));
						triplet.push_back(TripletX(c_offset + 10, a_offset + 11, 1.0 * x7 * 1.0));
						triplet.push_back(TripletX(c_offset + 11, a_offset + 11, 1.0 * x8 * 1.0));
					}
				}

			}
		}
	}

	qeal MipcConeConeConstraint::frictionEnergy()
	{
		Vector3 c11p, c12p, c21p, c22p;
		spheres[0]->center->projectFullspacePreP(c11p.data());
		spheres[1]->center->projectFullspacePreP(c12p.data());
		spheres[2]->center->projectFullspacePreP(c21p.data());
		spheres[3]->center->projectFullspacePreP(c22p.data());

		Vector3 rel_u;
		for (int i = 0; i < 3; i++)
			rel_u.data()[i] = lagAlpha * (spheres[0]->center->getP()[i] - c11p.data()[i]) + (1.0 - lagAlpha) * (spheres[1]->center->getP()[i] - c12p.data()[i]) - (lagBbeta * (spheres[2]->center->getP()[i] - c21p.data()[i]) + (1.0 - lagBbeta) * (spheres[3]->center->getP()[i] - c22p.data()[i]));

		Vector2 rel_uk = lagBasis.transpose() * rel_u;
		qeal energy;
		f0_SF(rel_uk.squaredNorm(), epsvh, energy);

		energy *= mu * lagLamda;

		return energy;
	}

	qeal MipcConeConeConstraint::computeDistance()
	{
		for (int i = 0; i < 3; i++)
		{
			sC1.data()[i] = spheres[0]->center->getP()[i] - spheres[1]->center->getP()[i];
			sC2.data()[i] = spheres[3]->center->getP()[i] - spheres[2]->center->getP()[i];
			sC3.data()[i] = spheres[1]->center->getP()[i] - spheres[3]->center->getP()[i];
		}

		sR1 = *(spheres[0]->radius) - *(spheres[1]->radius);
		sR2 = *(spheres[2]->radius) - *(spheres[3]->radius);
		sR3 = *(spheres[1]->radius) + *(spheres[3]->radius);

		A = sC1.dot(sC1) - sR1 * sR1;
		B = 2.0 * (sC1.dot(sC2) - sR1 * sR2);
		C = sC2.dot(sC2) - sR2 * sR2;
		D = 2.0 * (sC1.dot(sC3) - sR1 * sR3);
		E = 2.0 * (sC2.dot(sC3) - sR2 * sR3);
		F = sC3.dot(sC3) - sR3 * sR3;

		delta = 4 * A * C - B * B;

		qeal temp_alpha, temp_beta;
		qeal temp_dist;

		distanceMode = TWO_ENDPOINTS;
		alpha = 0.0, beta = 0.0;
		distance = valueOfQuadircSurface2D(alpha, beta, A, B, C, D, E, F);

		temp_alpha = 1.0, temp_beta = 0.0;
		temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
		if (distance > temp_dist)
		{
			distance = temp_dist;
			alpha = temp_alpha; beta = temp_beta;
		}

		temp_alpha = 0.0, temp_beta = 1.0;
		temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
		if (distance > temp_dist)
		{
			distance = temp_dist;
			alpha = temp_alpha; beta = temp_beta;
		}

		temp_alpha = 1.0, temp_beta = 1.0;
		temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
		if (distance > temp_dist)
		{
			distance = temp_dist;
			alpha = temp_alpha; beta = temp_beta;
		}

		temp_alpha = 0.0; temp_beta = -E / (2.0 * C);
		if (temp_beta > 0.0 && temp_beta < 1.0)
		{
			temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
			if (distance > temp_dist)
			{
				distance = temp_dist;
				alpha = temp_alpha; beta = temp_beta;
				distanceMode = ALPHA_ZERO;
			}
		}

		temp_alpha = 1.0; temp_beta = -(B + E) / (2.0 *C);
		if (temp_beta > 0.0 && temp_beta < 1.0)
		{
			temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
			if (distance > temp_dist)
			{
				distance = temp_dist;
				alpha = temp_alpha; beta = temp_beta;
				distanceMode = ALPHA_ONE;
			}
		}

		temp_alpha = -D / (2.0 *A); temp_beta = 0.0;
		if (temp_alpha > 0.0 && temp_alpha < 1.0)
		{
			temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
			if (distance > temp_dist)
			{
				distance = temp_dist;
				alpha = temp_alpha; beta = temp_beta;
				distanceMode = BETA_ZERO;
			}
		}

		temp_alpha = -(B + D) / (2.0 *A); temp_beta = 1.0;
		if (temp_alpha > 0.0 && temp_alpha < 1.0)
		{
			temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
			if (distance > temp_dist)
			{
				distance = temp_dist;
				alpha = temp_alpha; beta = temp_beta;
				distanceMode = BETA_ONE;
			}
		}

		if (delta != 0.0)
		{
			temp_alpha = (B * E - 2.0 * C * D) / delta; temp_beta = (B * D - 2.0 * A * E) / delta;
			if (temp_alpha > 0.0 && temp_alpha < 1.0 && temp_beta> 0.0 && temp_beta < 1.0)
			{
				temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
				if (distance > temp_dist)
				{
					distance = temp_dist;
					alpha = temp_alpha; beta = temp_beta;
					distanceMode = ALPHA_BETA;
				}
			}
		}

		//compute cloeset points;
		Vector3 cp, cq;
		for (int i = 0; i < 3; i++)
		{
			cp.data()[i] = alpha * spheres[0]->center->getP()[i] + (1.0 - alpha) * spheres[1]->center->getP()[i];
			cq.data()[i] = beta * spheres[2]->center->getP()[i] + (1.0 - beta) * spheres[3]->center->getP()[i];
		}
		qeal rp, rq;
		rp = alpha * spheres[0]->radius[0] + +(1.0 - alpha) * spheres[1]->radius[0];
		rq = beta * spheres[2]->radius[0] + (1.0 - beta) * spheres[3]->radius[0];

		Vector3 dir = cq - cp;
		dir.normalize();
		cloestPoints[0] = cp + dir * rp;
		cloestPoints[1] = cq - dir * rq;

		return distance;
	}

	void MipcConeConeConstraint::getTanBasis(Eigen::Matrix<qeal, 3, 2>& lagBasis)
	{
		lagBasis.setZero();
		switch (distanceMode)
		{
		case TWO_ENDPOINTS: //PP
		{
			Vector3 p;
			Vector3 q;

			if (alpha == 0.0 && beta == 0.0)
			{
				p = Vector3(spheres[1]->center->getP()[0], spheres[1]->center->getP()[1], spheres[1]->center->getP()[2]);
				q = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);
			}
			else if (alpha == 1.0 && beta == 0.0)
			{
				p = Vector3(spheres[0]->center->getP()[0], spheres[0]->center->getP()[1], spheres[0]->center->getP()[2]);
				q = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);
			}
			else if (alpha == 0.0 && beta == 1.0)
			{
				p = Vector3(spheres[1]->center->getP()[0], spheres[1]->center->getP()[1], spheres[1]->center->getP()[2]);
				q = Vector3(spheres[2]->center->getP()[0], spheres[2]->center->getP()[1], spheres[2]->center->getP()[2]);
			}
			else
			{
				p = Vector3(spheres[0]->center->getP()[0], spheres[0]->center->getP()[1], spheres[0]->center->getP()[2]);
				q = Vector3(spheres[2]->center->getP()[0], spheres[2]->center->getP()[1], spheres[2]->center->getP()[2]);
			}
			Point_Point_Tangent_Basis(p, q, lagBasis);
			break;
		}
		case ALPHA_ZERO: // PE
		{
			Vector3 p = Vector3(spheres[1]->center->getP()[0], spheres[1]->center->getP()[1], spheres[1]->center->getP()[2]);
			Vector3 e1 = Vector3(spheres[2]->center->getP()[0], spheres[2]->center->getP()[1], spheres[2]->center->getP()[2]);
			Vector3 e2 = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);
			Point_Edge_Tangent_Basis(p, e1, e2, lagBasis);
			break;
		}
		case BETA_ZERO: //PE
		{
			Vector3 p = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);
			Vector3 e1 = Vector3(spheres[0]->center->getP()[0], spheres[0]->center->getP()[1], spheres[0]->center->getP()[2]);
			Vector3 e2 = Vector3(spheres[1]->center->getP()[0], spheres[1]->center->getP()[1], spheres[1]->center->getP()[2]);
			Point_Edge_Tangent_Basis(p, e1, e2, lagBasis);
			break;
		}
		case ALPHA_ONE: //PE
		{
			Vector3 p = Vector3(spheres[0]->center->getP()[0], spheres[0]->center->getP()[1], spheres[0]->center->getP()[2]);
			Vector3 e1 = Vector3(spheres[2]->center->getP()[0], spheres[2]->center->getP()[1], spheres[2]->center->getP()[2]);
			Vector3 e2 = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);
			Point_Edge_Tangent_Basis(p, e1, e2, lagBasis);
			break;
		}	case BETA_ONE: //PE
		{
			Vector3 p = Vector3(spheres[2]->center->getP()[0], spheres[2]->center->getP()[1], spheres[2]->center->getP()[2]);
			Vector3 e1 = Vector3(spheres[0]->center->getP()[0], spheres[0]->center->getP()[1], spheres[0]->center->getP()[2]);
			Vector3 e2 = Vector3(spheres[1]->center->getP()[0], spheres[1]->center->getP()[1], spheres[1]->center->getP()[2]);
			Point_Edge_Tangent_Basis(p, e1, e2, lagBasis);
			break;
		}
		case ALPHA_BETA: //EE
		{
			Vector3 e0 = Vector3(spheres[0]->center->getP()[0], spheres[0]->center->getP()[1], spheres[0]->center->getP()[2]);
			Vector3 e1 = Vector3(spheres[1]->center->getP()[0], spheres[1]->center->getP()[1], spheres[1]->center->getP()[2]);
			Vector3 e2 = Vector3(spheres[2]->center->getP()[0], spheres[2]->center->getP()[1], spheres[2]->center->getP()[2]);
			Vector3 e3 = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);
			Edge_Edge_Tangent_Basis(e0, e1, e2, e3, lagBasis);
			break;
		}
		default:
			break;
		};
	}

	void MipcConeConeConstraint::computeLagTangentBasis(const qeal kappa)
	{
		getTanBasis(lagBasis);
		lagAlpha = alpha; lagBbeta = beta;
		lagLamda = -kappa * getBarrierGradient() * 2.0 * sqrt(distance);
	}

	void MipcConeConeConstraint::getGradientAndHessian(qeal kappa, VectorX& gradient, MatrixX& hessian)
	{
		qeal barrierGrad = getBarrierGradient();
		qeal barrierHessian = getBarrierHessian();

		VectorX distGrad;
		distGrad.resize(12);
		distGrad.setZero();

		MatrixX distHessina;
		distHessina.resize(12, 12);
		distHessina.setZero();

		diff_F_x(distGrad);
		//gradient
		qeal scale = -1.0 * kappa * barrierGrad;
		fillOverallGradient(scale, distGrad, gradient);

		//hessian
		switch (distanceMode)
		{
		case TWO_ENDPOINTS:
			endPointsHessina(distHessina);
			break;
		case ALPHA_ZERO:
			alphaIsZeroHessina(distHessina);
			break;
		case ALPHA_ONE:
			alphaIsOneHessina(distHessina);
			break;
		case BETA_ZERO:
			betaIsZeroHessina(distHessina);
			break;
		case BETA_ONE:
			betaIsOneHessina(distHessina);
			break;
		case ALPHA_BETA:
			alphaBetaHessina(distHessina);
			break;
		default:
			break;
		};

		MatrixX hess = barrierGrad * distHessina + barrierHessian * (distGrad * distGrad.transpose());
		makePD(hess);
		hess *= kappa;

		fillOverallHessian(1.0, hess, hessian);
	}

	void MipcConeConeConstraint::getFrictionGradientAndHessian(qeal kappa, VectorX& gradient, MatrixX& hessian)
	{
		Vector3 c11p, c12p, c21p, c22p;
		spheres[0]->center->projectFullspacePreP(c11p.data());
		spheres[1]->center->projectFullspacePreP(c12p.data());
		spheres[2]->center->projectFullspacePreP(c21p.data());
		spheres[3]->center->projectFullspacePreP(c22p.data());

		Vector3 rel_u;
		for (int i = 0; i < 3; i++)
			rel_u.data()[i] = (lagAlpha * (spheres[0]->center->getP()[i] - c11p.data()[i]) + (1.0 - lagAlpha) * (spheres[1]->center->getP()[i] - c12p.data()[i])) - (lagBbeta * (spheres[2]->center->getP()[i] - c21p.data()[i]) + (1.0 - lagBbeta) * (spheres[3]->center->getP()[i] - c22p.data()[i]));

		Vector2 rel_uk = lagBasis.transpose() * rel_u;
		qeal rel_ukSqNorm = rel_uk.squaredNorm();
		qeal rel_ukNorm = std::sqrt(rel_ukSqNorm);

		qeal f1_div_relDXNorm, f2_term;
		f1_SF_Div_RelDXNorm(rel_ukSqNorm, epsvh, f1_div_relDXNorm);
		f2_SF_Term(rel_ukSqNorm, epsvh, f2_term);


		Vector3 fricForce = -1.0 * f1_div_relDXNorm * mu *lagLamda * lagBasis * rel_uk;

		MatrixX HessianI(12, 12);
		HessianI.setZero();
		MatrixX selectMatrix(12, 3);
		selectMatrix.setZero();
		selectMatrix.data()[0] = lagAlpha;   selectMatrix.data()[13] = lagAlpha;   selectMatrix.data()[26] = lagAlpha;
		selectMatrix.data()[3] = 1 - lagAlpha;  selectMatrix.data()[16] = 1 - lagAlpha;  selectMatrix.data()[29] = 1 - lagAlpha;
		selectMatrix.data()[6] = -lagBbeta;   selectMatrix.data()[19] = -lagBbeta;   selectMatrix.data()[32] = -lagBbeta;
		selectMatrix.data()[9] = -(1.0 - lagBbeta); selectMatrix.data()[22] = -(1.0 - lagBbeta);  selectMatrix.data()[35] = -(1.0 - lagBbeta);

		MatrixX TT = lagBasis.transpose() * selectMatrix.transpose();

		if (rel_ukSqNorm >= (epsvh * epsvh)) {
			// no SPD projection needed
			Vector2 ubar;
			ubar.data()[0] = -rel_uk.data()[1];
			ubar.data()[1] = rel_uk.data()[0];
			HessianI = (TT.transpose() * ((mu * lagLamda * f1_div_relDXNorm / rel_ukSqNorm) * ubar)) * (ubar.transpose() * TT);
		}
		else
		{
			if (rel_ukSqNorm == 0) {
				HessianI = ((mu * lagLamda * f1_div_relDXNorm) * TT.transpose()) * TT;
			}
			else
			{
				MatrixX innerMtr = ((f2_term / rel_ukNorm) * rel_uk) * rel_uk.transpose();
				innerMtr.diagonal().array() += f1_div_relDXNorm;
				makePD(innerMtr);
				innerMtr *= mu * lagLamda;
				HessianI = TT.transpose() * innerMtr * TT;
			}
		}

		VectorX ff(12);
		ff.setZero();

		if (spheres[0]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[0]->center->getFrameId();
			ff.data()[0] = lagAlpha * fricForce.data()[0];
			ff.data()[1] = lagAlpha * fricForce.data()[1];
			ff.data()[2] = lagAlpha * fricForce.data()[2];
		}

		if (spheres[1]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[1]->center->getFrameId();
			ff.data()[3] = (1.0 - lagAlpha) *  fricForce.data()[0];
			ff.data()[4] = (1.0 - lagAlpha) *  fricForce.data()[1];
			ff.data()[5] = (1.0 - lagAlpha) *  fricForce.data()[2];
		}

		if (spheres[2]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[2]->center->getFrameId();
			ff.data()[6] = -lagBbeta * fricForce.data()[0];
			ff.data()[7] = -lagBbeta * fricForce.data()[1];
			ff.data()[8] = -lagBbeta * fricForce.data()[2];
		}

		if (spheres[3]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[3]->center->getFrameId();
			ff.data()[9] = -(1.0 - lagBbeta) * fricForce.data()[0];
			ff.data()[10] = -(1.0 - lagBbeta) * fricForce.data()[1];
			ff.data()[11] = -(1.0 - lagBbeta) * fricForce.data()[2];
		}
		fillOverallGradient(1.0, ff, gradient);
		fillOverallHessian(1.0, HessianI, hessian);
	}

	void MipcConeConeConstraint::getGradientAndHessian(qeal kappa, VectorX& gradient, std::vector<TripletX>& triplet)
	{
		Vector3 c11p, c12p, c21p, c22p;
		spheres[0]->center->projectFullspacePreP(c11p.data());
		spheres[1]->center->projectFullspacePreP(c12p.data());
		spheres[2]->center->projectFullspacePreP(c21p.data());
		spheres[3]->center->projectFullspacePreP(c22p.data());

		Vector3 rel_u;
		for (int i = 0; i < 3; i++)
			rel_u.data()[i] = (lagAlpha * (spheres[0]->center->getP()[i] - c11p.data()[i]) + (1.0 - lagAlpha) * (spheres[1]->center->getP()[i] - c12p.data()[i])) - (lagBbeta * (spheres[2]->center->getP()[i] - c21p.data()[i]) + (1.0 - lagBbeta) * (spheres[3]->center->getP()[i] - c22p.data()[i]));

		Vector2 rel_uk = lagBasis.transpose() * rel_u;
		qeal rel_ukSqNorm = rel_uk.squaredNorm();
		qeal rel_ukNorm = std::sqrt(rel_ukSqNorm);

		qeal f1_div_relDXNorm, f2_term;
		f1_SF_Div_RelDXNorm(rel_ukSqNorm, epsvh, f1_div_relDXNorm);
		f2_SF_Term(rel_ukSqNorm, epsvh, f2_term);


		Vector3 fricForce = -1.0 * f1_div_relDXNorm * mu *lagLamda * lagBasis * rel_uk;

		MatrixX HessianI(12, 12);
		HessianI.setZero();
		MatrixX selectMatrix(12, 3);
		selectMatrix.setZero();
		selectMatrix.data()[0] = lagAlpha;   selectMatrix.data()[13] = lagAlpha;   selectMatrix.data()[26] = lagAlpha;
		selectMatrix.data()[3] = 1 - lagAlpha;  selectMatrix.data()[16] = 1 - lagAlpha;  selectMatrix.data()[29] = 1 - lagAlpha;
		selectMatrix.data()[6] = -lagBbeta;   selectMatrix.data()[19] = -lagBbeta;   selectMatrix.data()[32] = -lagBbeta;
		selectMatrix.data()[9] = -(1.0 - lagBbeta); selectMatrix.data()[22] = -(1.0 - lagBbeta);  selectMatrix.data()[35] = -(1.0 - lagBbeta);

		MatrixX TT = lagBasis.transpose() * selectMatrix.transpose();

		if (rel_ukSqNorm >= (epsvh * epsvh)) {
			// no SPD projection needed
			Vector2 ubar;
			ubar.data()[0] = -rel_uk.data()[1];
			ubar.data()[1] = rel_uk.data()[0];
			HessianI = (TT.transpose() * ((mu * lagLamda * f1_div_relDXNorm / rel_ukSqNorm) * ubar)) * (ubar.transpose() * TT);
		}
		else
		{
			if (rel_ukSqNorm == 0) {
				HessianI = ((mu * lagLamda * f1_div_relDXNorm) * TT.transpose()) * TT;
			}
			else
			{
				MatrixX innerMtr = ((f2_term / rel_ukNorm) * rel_uk) * rel_uk.transpose();
				innerMtr.diagonal().array() += f1_div_relDXNorm;
				makePD(innerMtr);
				innerMtr *= mu * lagLamda;
				HessianI = TT.transpose() * innerMtr * TT;
			}
		}

		VectorX ff(12);
		ff.setZero();

		if (spheres[0]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[0]->center->getFrameId();
			ff.data()[0] = lagAlpha * fricForce.data()[0];
			ff.data()[1] = lagAlpha * fricForce.data()[1];
			ff.data()[2] = lagAlpha * fricForce.data()[2];
		}

		if (spheres[1]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[1]->center->getFrameId();
			ff.data()[3] = (1.0 - lagAlpha) *  fricForce.data()[0];
			ff.data()[4] = (1.0 - lagAlpha) *  fricForce.data()[1];
			ff.data()[5] = (1.0 - lagAlpha) *  fricForce.data()[2];
		}

		if (spheres[2]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[2]->center->getFrameId();
			ff.data()[6] = -lagBbeta * fricForce.data()[0];
			ff.data()[7] = -lagBbeta * fricForce.data()[1];
			ff.data()[8] = -lagBbeta * fricForce.data()[2];
		}

		if (spheres[3]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[3]->center->getFrameId();
			ff.data()[9] = -(1.0 - lagBbeta) * fricForce.data()[0];
			ff.data()[10] = -(1.0 - lagBbeta) * fricForce.data()[1];
			ff.data()[11] = -(1.0 - lagBbeta) * fricForce.data()[2];
		}
		fillOverallGradient(1.0, ff, gradient);
		fillOverallHessian(1.0, HessianI, triplet);
	}

	void MipcConeConeConstraint::getFrictionGradientAndHessian(qeal kappa, VectorX& gradient, std::vector<TripletX>& triplet)
	{
		Vector3 c11p, c12p, c21p, c22p;
		spheres[0]->center->projectFullspacePreP(c11p.data());
		spheres[1]->center->projectFullspacePreP(c12p.data());
		spheres[2]->center->projectFullspacePreP(c21p.data());
		spheres[3]->center->projectFullspacePreP(c22p.data());

		Vector3 rel_u;
		for (int i = 0; i < 3; i++)
			rel_u.data()[i] = (lagAlpha * (spheres[0]->center->getP()[i] - c11p.data()[i]) + (1.0 - lagAlpha) * (spheres[1]->center->getP()[i] - c12p.data()[i])) - (lagBbeta * (spheres[2]->center->getP()[i] - c21p.data()[i]) + (1.0 - lagBbeta) * (spheres[3]->center->getP()[i] - c22p.data()[i]));

		Vector2 rel_uk = lagBasis.transpose() * rel_u;
		qeal rel_ukSqNorm = rel_uk.squaredNorm();
		qeal rel_ukNorm = std::sqrt(rel_ukSqNorm);

		qeal f1_div_relDXNorm, f2_term;
		f1_SF_Div_RelDXNorm(rel_ukSqNorm, epsvh, f1_div_relDXNorm);
		f2_SF_Term(rel_ukSqNorm, epsvh, f2_term);


		Vector3 fricForce = -1.0 * f1_div_relDXNorm * mu *lagLamda * lagBasis * rel_uk;

		MatrixX HessianI(12, 12);
		HessianI.setZero();
		MatrixX selectMatrix(12, 3);
		selectMatrix.setZero();
		selectMatrix.data()[0] = lagAlpha;   selectMatrix.data()[13] = lagAlpha;   selectMatrix.data()[26] = lagAlpha;
		selectMatrix.data()[3] = 1 - lagAlpha;  selectMatrix.data()[16] = 1 - lagAlpha;  selectMatrix.data()[29] = 1 - lagAlpha;
		selectMatrix.data()[6] = -lagBbeta;   selectMatrix.data()[19] = -lagBbeta;   selectMatrix.data()[32] = -lagBbeta;
		selectMatrix.data()[9] = -(1.0 - lagBbeta); selectMatrix.data()[22] = -(1.0 - lagBbeta);  selectMatrix.data()[35] = -(1.0 - lagBbeta);

		MatrixX TT = lagBasis.transpose() * selectMatrix.transpose();

		if (rel_ukSqNorm >= (epsvh * epsvh)) {
			// no SPD projection needed
			Vector2 ubar;
			ubar.data()[0] = -rel_uk.data()[1];
			ubar.data()[1] = rel_uk.data()[0];
			HessianI = (TT.transpose() * ((mu * lagLamda * f1_div_relDXNorm / rel_ukSqNorm) * ubar)) * (ubar.transpose() * TT);
		}
		else
		{
			if (rel_ukSqNorm == 0) {
				HessianI = ((mu * lagLamda * f1_div_relDXNorm) * TT.transpose()) * TT;
			}
			else
			{
				MatrixX innerMtr = ((f2_term / rel_ukNorm) * rel_uk) * rel_uk.transpose();
				innerMtr.diagonal().array() += f1_div_relDXNorm;
				makePD(innerMtr);
				innerMtr *= mu * lagLamda;
				HessianI = TT.transpose() * innerMtr * TT;
			}
		}

		VectorX ff(12);
		ff.setZero();

		if (spheres[0]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[0]->center->getFrameId();
			ff.data()[0] = lagAlpha * fricForce.data()[0];
			ff.data()[1] = lagAlpha * fricForce.data()[1];
			ff.data()[2] = lagAlpha * fricForce.data()[2];
		}

		if (spheres[1]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[1]->center->getFrameId();
			ff.data()[3] = (1.0 - lagAlpha) *  fricForce.data()[0];
			ff.data()[4] = (1.0 - lagAlpha) *  fricForce.data()[1];
			ff.data()[5] = (1.0 - lagAlpha) *  fricForce.data()[2];
		}

		if (spheres[2]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[2]->center->getFrameId();
			ff.data()[6] = -lagBbeta * fricForce.data()[0];
			ff.data()[7] = -lagBbeta * fricForce.data()[1];
			ff.data()[8] = -lagBbeta * fricForce.data()[2];
		}

		if (spheres[3]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[3]->center->getFrameId();
			ff.data()[9] = -(1.0 - lagBbeta) * fricForce.data()[0];
			ff.data()[10] = -(1.0 - lagBbeta) * fricForce.data()[1];
			ff.data()[11] = -(1.0 - lagBbeta) * fricForce.data()[2];
		}
		fillOverallGradient(1.0, ff, gradient);
		fillOverallHessian(1.0, HessianI, triplet);
	}

	void MipcConeConeConstraint::diff_F_x(VectorX& diff)
	{
		Vector3 v = 2.0 * (alpha * sC1 + beta * sC2 + sC3);
		diff[0] = alpha * v.data()[0];
		diff[1] = alpha * v.data()[1];
		diff[2] = alpha * v.data()[2];

		diff[3] = (1.0 - alpha) * v.data()[0];
		diff[4] = (1.0 - alpha) * v.data()[1];
		diff[5] = (1.0 - alpha) * v.data()[2];

		diff[6] = -beta * v.data()[0];
		diff[7] = -beta * v.data()[1];
		diff[8] = -beta * v.data()[2];

		diff[9] = -(1.0 - beta) * v.data()[0];
		diff[10] = -(1.0 - beta) * v.data()[1];
		diff[11] = -(1.0 - beta) * v.data()[2];
	}

	void MipcConeConeConstraint::endPointsHessina(MatrixX & hessina)
	{
		// alpha and beta is constant
		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (alpha * C1x + beta * C2x + C3x);
		qeal vy = 2.0 * (alpha * C1y + beta * C2y + C3y);
		qeal vz = 2.0 * (alpha * C1z + beta * C2z + C3z);

		qeal dvxdc11_x = 2.0 * alpha;
		qeal dvxdc11_y = 0;
		qeal dvxdc11_z = 0;

		qeal dvxdc12_x = 2.0 * (1 - alpha);
		qeal dvxdc12_y = 0;
		qeal dvxdc12_z = 0;

		qeal dvxdc21_x = 2.0 * -beta;
		qeal dvxdc21_y = 0;
		qeal dvxdc21_z = 0;

		qeal dvxdc22_x = 2.0 * (beta - 1.0);
		qeal dvxdc22_y = 0;
		qeal dvxdc22_z = 0;
		//
		qeal dvydc11_x = 0;
		qeal dvydc11_y = 2.0 * alpha;
		qeal dvydc11_z = 0;

		qeal dvydc12_x = 0;
		qeal dvydc12_y = 2.0 * (1 - alpha);
		qeal dvydc12_z = 0;

		qeal dvydc21_x = 0;
		qeal dvydc21_y = 2.0 * -beta;
		qeal dvydc21_z = 0;

		qeal dvydc22_x = 0;
		qeal dvydc22_y = 2.0 * (beta - 1.0);
		qeal dvydc22_z = 0;
		//
		qeal dvzdc11_x = 0;
		qeal dvzdc11_y = 0;
		qeal dvzdc11_z = 2.0 * alpha;

		qeal dvzdc12_x = 0;
		qeal dvzdc12_y = 0;
		qeal dvzdc12_z = 2.0 * (1 - alpha);

		qeal dvzdc21_x = 0;
		qeal dvzdc21_y = 0;
		qeal dvzdc21_z = 2.0 * -beta;

		qeal dvzdc22_x = 0;
		qeal dvzdc22_y = 0;
		qeal dvzdc22_z = 2.0 * (beta - 1.0);

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = alpha * dvxdc11_x;
		hessina.data()[1] = alpha * dvydc11_x;
		hessina.data()[2] = alpha * dvzdc11_x;

		hessina.data()[3] = (1 - alpha) * dvxdc11_x;
		hessina.data()[4] = (1 - alpha) * dvydc11_x;
		hessina.data()[5] = (1 - alpha) * dvzdc11_x;

		hessina.data()[6] = -beta * dvxdc11_x;
		hessina.data()[7] = -beta * dvydc11_x;
		hessina.data()[8] = -beta * dvzdc11_x;

		hessina.data()[9] = (beta - 1.0) * dvxdc11_x;
		hessina.data()[10] = (beta - 1.0) * dvydc11_x;
		hessina.data()[11] = (beta - 1.0) * dvzdc11_x;
		// f / c11y
		hessina.data()[12] = alpha * dvxdc11_y;
		hessina.data()[13] = alpha * dvydc11_y;
		hessina.data()[14] = alpha * dvzdc11_y;

		hessina.data()[15] = (1 - alpha) * dvxdc11_y;
		hessina.data()[16] = (1 - alpha) * dvydc11_y;
		hessina.data()[17] = (1 - alpha) * dvzdc11_y;

		hessina.data()[18] = -beta * dvxdc11_y;
		hessina.data()[19] = -beta * dvydc11_y;
		hessina.data()[20] = -beta * dvzdc11_y;

		hessina.data()[21] = (beta - 1.0) * dvxdc11_y;
		hessina.data()[22] = (beta - 1.0) * dvydc11_y;
		hessina.data()[23] = (beta - 1.0) * dvzdc11_y;
		// f / c11z
		hessina.data()[24] = alpha * dvxdc11_z;
		hessina.data()[25] = alpha * dvydc11_z;
		hessina.data()[26] = alpha * dvzdc11_z;

		hessina.data()[27] = (1 - alpha) * dvxdc11_z;
		hessina.data()[28] = (1 - alpha) * dvydc11_z;
		hessina.data()[29] = (1 - alpha) * dvzdc11_z;

		hessina.data()[30] = -beta * dvxdc11_z;
		hessina.data()[31] = -beta * dvydc11_z;
		hessina.data()[32] = -beta * dvzdc11_z;

		hessina.data()[33] = (beta - 1.0) * dvxdc11_z;
		hessina.data()[34] = (beta - 1.0) * dvydc11_z;
		hessina.data()[35] = (beta - 1.0) * dvzdc11_z;

		// f / c12x
		hessina.data()[36] = alpha * dvxdc12_x;
		hessina.data()[37] = alpha * dvydc12_x;
		hessina.data()[38] = alpha * dvzdc12_x;

		hessina.data()[39] = (1 - alpha) * dvxdc12_x;
		hessina.data()[40] = (1 - alpha) * dvydc12_x;
		hessina.data()[41] = (1 - alpha) * dvzdc12_x;

		hessina.data()[42] = -beta * dvxdc12_x;
		hessina.data()[43] = -beta * dvydc12_x;
		hessina.data()[44] = -beta * dvzdc12_x;

		hessina.data()[45] = (beta - 1.0) * dvxdc12_x;
		hessina.data()[46] = (beta - 1.0) * dvydc12_x;
		hessina.data()[47] = (beta - 1.0) * dvzdc12_x;
		// f / c12y
		hessina.data()[48] = alpha * dvxdc12_y;
		hessina.data()[49] = alpha * dvydc12_y;
		hessina.data()[50] = alpha * dvzdc12_y;

		hessina.data()[51] = (1 - alpha) * dvxdc12_y;
		hessina.data()[52] = (1 - alpha) * dvydc12_y;
		hessina.data()[53] = (1 - alpha) * dvzdc12_y;

		hessina.data()[54] = -beta * dvxdc12_y;
		hessina.data()[55] = -beta * dvydc12_y;
		hessina.data()[56] = -beta * dvzdc12_y;

		hessina.data()[57] = (beta - 1.0) * dvxdc12_y;
		hessina.data()[58] = (beta - 1.0) * dvydc12_y;
		hessina.data()[59] = (beta - 1.0) * dvzdc12_y;
		// f / c12z
		hessina.data()[60] = alpha * dvxdc12_z;
		hessina.data()[61] = alpha * dvydc12_z;
		hessina.data()[62] = alpha * dvzdc12_z;

		hessina.data()[63] = (1 - alpha) * dvxdc12_z;
		hessina.data()[64] = (1 - alpha) * dvydc12_z;
		hessina.data()[65] = (1 - alpha) * dvzdc12_z;

		hessina.data()[66] = -beta * dvxdc12_z;
		hessina.data()[67] = -beta * dvydc12_z;
		hessina.data()[68] = -beta * dvzdc12_z;

		hessina.data()[69] = (beta - 1.0) * dvxdc12_z;
		hessina.data()[70] = (beta - 1.0) * dvydc12_z;
		hessina.data()[71] = (beta - 1.0) * dvzdc12_z;

		// f / c21x
		hessina.data()[72] = alpha * dvxdc21_x;
		hessina.data()[73] = alpha * dvydc21_x;
		hessina.data()[74] = alpha * dvzdc21_x;

		hessina.data()[75] = (1 - alpha) * dvxdc21_x;
		hessina.data()[76] = (1 - alpha) * dvydc21_x;
		hessina.data()[77] = (1 - alpha) * dvzdc21_x;

		hessina.data()[78] = -beta * dvxdc21_x;
		hessina.data()[79] = -beta * dvydc21_x;
		hessina.data()[80] = -beta * dvzdc21_x;

		hessina.data()[81] = (beta - 1.0) * dvxdc21_x;
		hessina.data()[82] = (beta - 1.0) * dvydc21_x;
		hessina.data()[83] = (beta - 1.0) * dvzdc21_x;
		// f / c21y
		hessina.data()[84] = alpha * dvxdc21_y;
		hessina.data()[85] = alpha * dvydc21_y;
		hessina.data()[86] = alpha * dvzdc21_y;

		hessina.data()[87] = (1 - alpha) * dvxdc21_y;
		hessina.data()[88] = (1 - alpha) * dvydc21_y;
		hessina.data()[89] = (1 - alpha) * dvzdc21_y;

		hessina.data()[90] = -beta * dvxdc21_y;
		hessina.data()[91] = -beta * dvydc21_y;
		hessina.data()[92] = -beta * dvzdc21_y;

		hessina.data()[93] = (beta - 1.0) * dvxdc21_y;
		hessina.data()[94] = (beta - 1.0) * dvydc21_y;
		hessina.data()[95] = (beta - 1.0) * dvzdc21_y;
		// f / c21z
		hessina.data()[96] = alpha * dvxdc21_z;
		hessina.data()[97] = alpha * dvydc21_z;
		hessina.data()[98] = alpha * dvzdc21_z;

		hessina.data()[99] = (1 - alpha) * dvxdc21_z;
		hessina.data()[100] = (1 - alpha) * dvydc21_z;
		hessina.data()[101] = (1 - alpha) * dvzdc21_z;

		hessina.data()[102] = -beta * dvxdc21_z;
		hessina.data()[103] = -beta * dvydc21_z;
		hessina.data()[104] = -beta * dvzdc21_z;

		hessina.data()[105] = (beta - 1.0) * dvxdc21_z;
		hessina.data()[106] = (beta - 1.0) * dvydc21_z;
		hessina.data()[107] = (beta - 1.0) * dvzdc21_z;

		// f / c22x
		hessina.data()[108] = alpha * dvxdc22_x;
		hessina.data()[109] = alpha * dvydc22_x;
		hessina.data()[110] = alpha * dvzdc22_x;

		hessina.data()[111] = (1 - alpha) * dvxdc22_x;
		hessina.data()[112] = (1 - alpha) * dvydc22_x;
		hessina.data()[113] = (1 - alpha) * dvzdc22_x;

		hessina.data()[114] = -beta * dvxdc22_x;
		hessina.data()[115] = -beta * dvydc22_x;
		hessina.data()[116] = -beta * dvzdc22_x;

		hessina.data()[117] = (beta - 1.0) * dvxdc22_x;
		hessina.data()[118] = (beta - 1.0) * dvydc22_x;
		hessina.data()[119] = (beta - 1.0) * dvzdc22_x;
		// f / c22y
		hessina.data()[120] = alpha * dvxdc22_y;
		hessina.data()[121] = alpha * dvydc22_y;
		hessina.data()[122] = alpha * dvzdc22_y;

		hessina.data()[123] = (1 - alpha) * dvxdc22_y;
		hessina.data()[124] = (1 - alpha) * dvydc22_y;
		hessina.data()[125] = (1 - alpha) * dvzdc22_y;

		hessina.data()[126] = -beta * dvxdc22_y;
		hessina.data()[127] = -beta * dvydc22_y;
		hessina.data()[128] = -beta * dvzdc22_y;

		hessina.data()[129] = (beta - 1.0) * dvxdc22_y;
		hessina.data()[130] = (beta - 1.0) * dvydc22_y;
		hessina.data()[131] = (beta - 1.0) * dvzdc22_y;
		// f / c22z
		hessina.data()[132] = alpha * dvxdc22_z;
		hessina.data()[133] = alpha * dvydc22_z;
		hessina.data()[134] = alpha * dvzdc22_z;

		hessina.data()[135] = (1 - alpha) * dvxdc22_z;
		hessina.data()[136] = (1 - alpha) * dvydc22_z;
		hessina.data()[137] = (1 - alpha) * dvzdc22_z;

		hessina.data()[138] = -beta * dvxdc22_z;
		hessina.data()[139] = -beta * dvydc22_z;
		hessina.data()[140] = -beta * dvzdc22_z;

		hessina.data()[141] = (beta - 1.0) * dvxdc22_z;
		hessina.data()[142] = (beta - 1.0) * dvydc22_z;
		hessina.data()[143] = (beta - 1.0) * dvzdc22_z;
	}

	void MipcConeConeConstraint::alphaIsZeroHessina(MatrixX & hessina)
	{
		qeal dAdc11_x = 2.0 * sC1.data()[0]; qeal dAdc11_y = 2.0 * sC1.data()[1]; qeal dAdc11_z = 2.0 * sC1.data()[2];
		qeal dAdc12_x = -2.0 * sC1.data()[0]; qeal dAdc12_y = -2.0 * sC1.data()[1]; qeal dAdc12_z = -2.0 * sC1.data()[2];

		qeal dBdc11_x = 2.0 * sC2.data()[0]; qeal dBdc11_y = 2.0 * sC2.data()[1]; qeal dBdc11_z = 2.0 * sC2.data()[2];
		qeal dBdc12_x = -2.0 * sC2.data()[0]; qeal dBdc12_y = -2.0 * sC2.data()[1]; qeal dBdc12_z = -2.0 * sC2.data()[2];
		qeal dBdc21_x = -2.0 * sC1.data()[0]; qeal dBdc21_y = -2.0 * sC1.data()[1]; qeal dBdc21_z = -2.0 * sC1.data()[2];
		qeal dBdc22_x = 2.0 * sC1.data()[0]; qeal dBdc22_y = 2.0 * sC1.data()[1]; qeal dBdc22_z = 2.0 * sC1.data()[2];

		qeal dCdc21_x = -2.0 * sC2.data()[0]; qeal dCdc21_y = -2.0 * sC2.data()[1]; qeal dCdc21_z = -2.0 * sC2.data()[2];
		qeal dCdc22_x = 2.0 * sC2.data()[0]; qeal dCdc22_y = 2.0 * sC2.data()[1]; qeal dCdc22_z = 2.0 * sC2.data()[2];

		qeal dDdc11_x = 2.0 * sC3.data()[0]; qeal dDdc11_y = 2.0 * sC3.data()[1]; qeal dDdc11_z = 2.0 * sC3.data()[2];
		qeal dDdc12_x = 2.0 * (sC1.data()[0] - sC3.data()[0]); qeal dDdc12_y = 2.0 * (sC1.data()[1] - sC3.data()[1]); qeal dDdc12_z = 2.0 * (sC1.data()[2] - sC3.data()[2]);
		qeal dDdc22_x = -2.0 * sC1.data()[0]; qeal dDdc22_y = -2.0 * sC1.data()[1]; qeal dDdc22_z = -2.0 * sC1.data()[2];

		qeal dEdc12_x = 2.0 * sC2.data()[0]; qeal dEdc12_y = 2.0 * sC2.data()[1]; qeal dEdc12_z = 2.0 * sC2.data()[2];
		qeal dEdc21_x = -2.0 * sC3.data()[0]; qeal dEdc21_y = -2.0 * sC3.data()[1]; qeal dEdc21_z = -2.0 * sC3.data()[2];
		qeal dEdc22_x = 2.0 * (sC3.data()[0] - sC2.data()[0]); qeal dEdc22_y = 2.0 * (sC3.data()[1] - sC2.data()[1]); qeal dEdc22_z = 2.0 * (sC3.data()[2] - sC2.data()[2]);

		qeal dFdc12_x = 2.0 * sC3.data()[0]; qeal dFdc12_y = 2.0 * sC3.data()[1]; qeal dFdc12_z = 2.0 * sC3.data()[2];
		qeal dFdc22 = -2.0 * sC3.data()[0]; qeal dFdc22_y = -2.0 * sC3.data()[1]; qeal dFd22_z = -2.0 * sC3.data()[2];

		//
		qeal q1 = 0.5 * E / (C * C);  qeal q2 = -0.5 / C;

		qeal dBTdc12_x = q2 * dEdc12_x;
		qeal dBTdc12_y = q2 * dEdc12_y;
		qeal dBTdc12_z = q2 * dEdc12_z;

		qeal dBTdc21_x = q1 * dCdc21_x + q2 * dEdc21_x;
		qeal dBTdc21_y = q1 * dCdc21_y + q2 * dEdc21_y;
		qeal dBTdc21_z = q1 * dCdc21_z + q2 * dEdc21_z;

		qeal dBTdc22_x = q1 * dCdc22_x + q2 * dEdc22_x;
		qeal dBTdc22_y = q1 * dCdc22_y + q2 * dEdc22_y;
		qeal dBTdc22_z = q1 * dCdc22_z + q2 * dEdc22_z;

		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (beta * C2x + C3x);
		qeal vy = 2.0 * (beta * C2y + C3y);
		qeal vz = 2.0 * (beta * C2z + C3z);

		qeal dvxdc12_x = 2.0 * (C2x * dBTdc12_x + 1);
		qeal dvxdc12_y = 2.0 * (C2x * dBTdc12_y);
		qeal dvxdc12_z = 2.0 * (C2x * dBTdc12_z);

		qeal dvxdc21_x = 2.0 * (C2x * dBTdc21_x - beta);
		qeal dvxdc21_y = 2.0 * (C2x * dBTdc21_y);
		qeal dvxdc21_z = 2.0 * (C2x * dBTdc21_z);

		qeal dvxdc22_x = 2.0 * (C2x * dBTdc22_x + beta - 1);
		qeal dvxdc22_y = 2.0 * (C2x * dBTdc22_y);
		qeal dvxdc22_z = 2.0 * (C2x * dBTdc22_z);
		//
		qeal dvydc12_x = 2.0 * (C2y * dBTdc12_x);
		qeal dvydc12_y = 2.0 * (C2y * dBTdc12_y + 1);
		qeal dvydc12_z = 2.0 * (C2y * dBTdc12_z);

		qeal dvydc21_x = 2.0 * (C2y * dBTdc21_x);
		qeal dvydc21_y = 2.0 * (C2y * dBTdc21_y - beta);
		qeal dvydc21_z = 2.0 * (C2y * dBTdc21_z);

		qeal dvydc22_x = 2.0 * (C2y * dBTdc22_x);
		qeal dvydc22_y = 2.0 * (C2y * dBTdc22_y + beta - 1);
		qeal dvydc22_z = 2.0 * (C2y * dBTdc22_z);
		//

		qeal dvzdc12_x = 2.0 * (C2z * dBTdc12_x);
		qeal dvzdc12_y = 2.0 * (C2z * dBTdc12_y);
		qeal dvzdc12_z = 2.0 * (C2z * dBTdc12_z + 1);

		qeal dvzdc21_x = 2.0 * (C2z * dBTdc21_x);
		qeal dvzdc21_y = 2.0 * (C2z * dBTdc21_y);
		qeal dvzdc21_z = 2.0 * (C2z * dBTdc21_z - beta);

		qeal dvzdc22_x = 2.0 * (C2z * dBTdc22_x);
		qeal dvzdc22_y = 2.0 * (C2z * dBTdc22_y);
		qeal dvzdc22_z = 2.0 * (C2z * dBTdc22_z + beta - 1);

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = 0;
		hessina.data()[1] = 0;
		hessina.data()[2] = 0;

		hessina.data()[3] = 0;
		hessina.data()[4] = 0;
		hessina.data()[5] = 0;

		hessina.data()[6] = 0;
		hessina.data()[7] = 0;
		hessina.data()[8] = 0;

		hessina.data()[9] = 0;
		hessina.data()[10] = 0;
		hessina.data()[11] = 0;
		// f / c11y
		hessina.data()[12] = 0;
		hessina.data()[13] = 0;
		hessina.data()[14] = 0;

		hessina.data()[15] = 0;
		hessina.data()[16] = 0;
		hessina.data()[17] = 0;

		hessina.data()[18] = 0;
		hessina.data()[19] = 0;
		hessina.data()[20] = 0;

		hessina.data()[21] = 0;
		hessina.data()[22] = 0;
		hessina.data()[23] = 0;
		// f / c11z
		hessina.data()[24] = 0;
		hessina.data()[25] = 0;
		hessina.data()[26] = 0;

		hessina.data()[27] = 0;
		hessina.data()[28] = 0;
		hessina.data()[29] = 0;

		hessina.data()[30] = 0;
		hessina.data()[31] = 0;
		hessina.data()[32] = 0;

		hessina.data()[33] = 0;
		hessina.data()[34] = 0;
		hessina.data()[35] = 0;

		// f / c12x
		hessina.data()[36] = 0;
		hessina.data()[37] = 0;
		hessina.data()[38] = 0;

		hessina.data()[39] = dvxdc12_x;
		hessina.data()[40] = dvydc12_x;
		hessina.data()[41] = dvzdc12_x;

		hessina.data()[42] = vx * -dBTdc12_x - beta * dvxdc12_x;
		hessina.data()[43] = vy * -dBTdc12_x - beta * dvydc12_x;
		hessina.data()[44] = vz * -dBTdc12_x - beta * dvzdc12_x;

		hessina.data()[45] = vx * dBTdc12_x + (beta - 1.0) * dvxdc12_x;
		hessina.data()[46] = vy * dBTdc12_x + (beta - 1.0) * dvydc12_x;
		hessina.data()[47] = vz * dBTdc12_x + (beta - 1.0) * dvzdc12_x;
		// f / c12y
		hessina.data()[48] = 0;
		hessina.data()[49] = 0;
		hessina.data()[50] = 0;

		hessina.data()[51] = dvxdc12_y;
		hessina.data()[52] = dvydc12_y;
		hessina.data()[53] = dvzdc12_y;

		hessina.data()[54] = vx * -dBTdc12_y - beta * dvxdc12_y;
		hessina.data()[55] = vy * -dBTdc12_y - beta * dvydc12_y;
		hessina.data()[56] = vz * -dBTdc12_y - beta * dvzdc12_y;

		hessina.data()[57] = vx * dBTdc12_y + (beta - 1.0) * dvxdc12_y;
		hessina.data()[58] = vy * dBTdc12_y + (beta - 1.0) * dvydc12_y;
		hessina.data()[59] = vz * dBTdc12_y + (beta - 1.0) * dvzdc12_y;
		// f / c12z
		hessina.data()[60] = 0;
		hessina.data()[61] = 0;
		hessina.data()[62] = 0;

		hessina.data()[63] = dvxdc12_z;
		hessina.data()[64] = dvydc12_z;
		hessina.data()[65] = dvzdc12_z;

		hessina.data()[66] = vx * -dBTdc12_z - beta * dvxdc12_z;
		hessina.data()[67] = vy * -dBTdc12_z - beta * dvydc12_z;
		hessina.data()[68] = vz * -dBTdc12_z - beta * dvzdc12_z;

		hessina.data()[69] = vx * dBTdc12_z + (beta - 1.0) * dvxdc12_z;
		hessina.data()[70] = vy * dBTdc12_z + (beta - 1.0) * dvydc12_z;
		hessina.data()[71] = vz * dBTdc12_z + (beta - 1.0) * dvzdc12_z;

		// f / c21x
		hessina.data()[72] = 0;
		hessina.data()[73] = 0;
		hessina.data()[74] = 0;

		hessina.data()[75] = dvxdc21_x;
		hessina.data()[76] = dvydc21_x;
		hessina.data()[77] = dvzdc21_x;

		hessina.data()[78] = vx * -dBTdc21_x - beta * dvxdc21_x;
		hessina.data()[79] = vy * -dBTdc21_x - beta * dvydc21_x;
		hessina.data()[80] = vz * -dBTdc21_x - beta * dvzdc21_x;

		hessina.data()[81] = vx * dBTdc21_x + (beta - 1.0) * dvxdc21_x;
		hessina.data()[82] = vy * dBTdc21_x + (beta - 1.0) * dvydc21_x;
		hessina.data()[83] = vz * dBTdc21_x + (beta - 1.0) * dvzdc21_x;
		// f / c21y
		hessina.data()[84] = 0;
		hessina.data()[85] = 0;
		hessina.data()[86] = 0;

		hessina.data()[87] = dvxdc21_y;
		hessina.data()[88] = dvydc21_y;
		hessina.data()[89] = dvzdc21_y;

		hessina.data()[90] = vx * -dBTdc21_y - beta * dvxdc21_y;
		hessina.data()[91] = vy * -dBTdc21_y - beta * dvydc21_y;
		hessina.data()[92] = vz * -dBTdc21_y - beta * dvzdc21_y;

		hessina.data()[93] = vx * dBTdc21_y + (beta - 1.0) * dvxdc21_y;
		hessina.data()[94] = vy * dBTdc21_y + (beta - 1.0) * dvydc21_y;
		hessina.data()[95] = vz * dBTdc21_y + (beta - 1.0) * dvzdc21_y;
		// f / c21z
		hessina.data()[96] = 0;
		hessina.data()[97] = 0;
		hessina.data()[98] = 0;

		hessina.data()[99] = dvxdc21_z;
		hessina.data()[100] = dvydc21_z;
		hessina.data()[101] = dvzdc21_z;

		hessina.data()[102] = vx * -dBTdc21_z - beta * dvxdc21_z;
		hessina.data()[103] = vy * -dBTdc21_z - beta * dvydc21_z;
		hessina.data()[104] = vz * -dBTdc21_z - beta * dvzdc21_z;

		hessina.data()[105] = vx * dBTdc21_z + (beta - 1.0) * dvxdc21_z;
		hessina.data()[106] = vy * dBTdc21_z + (beta - 1.0) * dvydc21_z;
		hessina.data()[107] = vz * dBTdc21_z + (beta - 1.0) * dvzdc21_z;

		// f / c22x
		hessina.data()[108] = 0;
		hessina.data()[109] = 0;
		hessina.data()[110] = 0;

		hessina.data()[111] = dvxdc22_x;
		hessina.data()[112] = dvydc22_x;
		hessina.data()[113] = dvzdc22_x;

		hessina.data()[114] = vx * -dBTdc22_x - beta * dvxdc22_x;
		hessina.data()[115] = vy * -dBTdc22_x - beta * dvydc22_x;
		hessina.data()[116] = vz * -dBTdc22_x - beta * dvzdc22_x;

		hessina.data()[117] = vx * dBTdc22_x + (beta - 1.0) * dvxdc22_x;
		hessina.data()[118] = vy * dBTdc22_x + (beta - 1.0) * dvydc22_x;
		hessina.data()[119] = vz * dBTdc22_x + (beta - 1.0) * dvzdc22_x;
		// f / c22y
		hessina.data()[120] = 0;
		hessina.data()[121] = 0;
		hessina.data()[122] = 0;

		hessina.data()[123] = dvxdc22_y;
		hessina.data()[124] = dvydc22_y;
		hessina.data()[125] = dvzdc22_y;

		hessina.data()[126] = vx * -dBTdc22_y - beta * dvxdc22_y;
		hessina.data()[127] = vy * -dBTdc22_y - beta * dvydc22_y;
		hessina.data()[128] = vz * -dBTdc22_y - beta * dvzdc22_y;

		hessina.data()[129] = vx * dBTdc22_y + (beta - 1.0) * dvxdc22_y;
		hessina.data()[130] = vy * dBTdc22_y + (beta - 1.0) * dvydc22_y;
		hessina.data()[131] = vz * dBTdc22_y + (beta - 1.0) * dvzdc22_y;
		// f / c22z
		hessina.data()[132] = 0;
		hessina.data()[133] = 0;
		hessina.data()[134] = 0;

		hessina.data()[135] = dvxdc22_z;
		hessina.data()[136] = dvydc22_z;
		hessina.data()[137] = dvzdc22_z;

		hessina.data()[138] = vx * -dBTdc22_z - beta * dvxdc22_z;
		hessina.data()[139] = vy * -dBTdc22_z - beta * dvydc22_z;
		hessina.data()[140] = vz * -dBTdc22_z - beta * dvzdc22_z;

		hessina.data()[141] = vx * dBTdc22_z + (beta - 1.0) * dvxdc22_z;
		hessina.data()[142] = vy * dBTdc22_z + (beta - 1.0) * dvydc22_z;
		hessina.data()[143] = vz * dBTdc22_z + (beta - 1.0) * dvzdc22_z;

	}

	void MipcConeConeConstraint::alphaIsOneHessina(MatrixX & hessina)
	{
		qeal dAdc11_x = 2.0 * sC1.data()[0]; qeal dAdc11_y = 2.0 * sC1.data()[1]; qeal dAdc11_z = 2.0 * sC1.data()[2];
		qeal dAdc12_x = -2.0 * sC1.data()[0]; qeal dAdc12_y = -2.0 * sC1.data()[1]; qeal dAdc12_z = -2.0 * sC1.data()[2];

		qeal dBdc11_x = 2.0 * sC2.data()[0]; qeal dBdc11_y = 2.0 * sC2.data()[1]; qeal dBdc11_z = 2.0 * sC2.data()[2];
		qeal dBdc12_x = -2.0 * sC2.data()[0]; qeal dBdc12_y = -2.0 * sC2.data()[1]; qeal dBdc12_z = -2.0 * sC2.data()[2];
		qeal dBdc21_x = -2.0 * sC1.data()[0]; qeal dBdc21_y = -2.0 * sC1.data()[1]; qeal dBdc21_z = -2.0 * sC1.data()[2];
		qeal dBdc22_x = 2.0 * sC1.data()[0]; qeal dBdc22_y = 2.0 * sC1.data()[1]; qeal dBdc22_z = 2.0 * sC1.data()[2];

		qeal dCdc21_x = -2.0 * sC2.data()[0]; qeal dCdc21_y = -2.0 * sC2.data()[1]; qeal dCdc21_z = -2.0 * sC2.data()[2];
		qeal dCdc22_x = 2.0 * sC2.data()[0]; qeal dCdc22_y = 2.0 * sC2.data()[1]; qeal dCdc22_z = 2.0 * sC2.data()[2];

		qeal dDdc11_x = 2.0 * sC3.data()[0]; qeal dDdc11_y = 2.0 * sC3.data()[1]; qeal dDdc11_z = 2.0 * sC3.data()[2];
		qeal dDdc12_x = 2.0 * (sC1.data()[0] - sC3.data()[0]); qeal dDdc12_y = 2.0 * (sC1.data()[1] - sC3.data()[1]); qeal dDdc12_z = 2.0 * (sC1.data()[2] - sC3.data()[2]);
		qeal dDdc22_x = -2.0 * sC1.data()[0]; qeal dDdc22_y = -2.0 * sC1.data()[1]; qeal dDdc22_z = -2.0 * sC1.data()[2];

		qeal dEdc12_x = 2.0 * sC2.data()[0]; qeal dEdc12_y = 2.0 * sC2.data()[1]; qeal dEdc12_z = 2.0 * sC2.data()[2];
		qeal dEdc21_x = -2.0 * sC3.data()[0]; qeal dEdc21_y = -2.0 * sC3.data()[1]; qeal dEdc21_z = -2.0 * sC3.data()[2];
		qeal dEdc22_x = 2.0 * (sC3.data()[0] - sC2.data()[0]); qeal dEdc22_y = 2.0 * (sC3.data()[1] - sC2.data()[1]); qeal dEdc22_z = 2.0 * (sC3.data()[2] - sC2.data()[2]);

		qeal dFdc12_x = 2.0 * sC3.data()[0]; qeal dFdc12_y = 2.0 * sC3.data()[1]; qeal dFdc12_z = 2.0 * sC3.data()[2];
		qeal dFdc22 = -2.0 * sC3.data()[0]; qeal dFdc22_y = -2.0 * sC3.data()[1]; qeal dFd22_z = -2.0 * sC3.data()[2];

		//
		qeal q1 = 0.5 * (B + E) / (C * C);  qeal q2 = -0.5 / C;

		qeal dALdc11_x = 0;
		qeal dALdc11_y = 0;
		qeal dALdc11_z = 0;

		qeal dALdc12_x = 0;
		qeal dALdc12_y = 0;
		qeal dALdc12_z = 0;

		qeal dALdc21_x = 0;
		qeal dALdc21_y = 0;
		qeal dALdc21_z = 0;

		qeal dALdc22_x = 0;
		qeal dALdc22_y = 0;
		qeal dALdc22_z = 0;

		qeal dBTdc11_x = q2 * dBdc11_x;
		qeal dBTdc11_y = q2 * dBdc11_y;
		qeal dBTdc11_z = q2 * dBdc11_z;

		qeal dBTdc12_x = q2 * (dBdc12_x + dEdc12_x);
		qeal dBTdc12_y = q2 * (dBdc12_y + dEdc12_y);
		qeal dBTdc12_z = q2 * (dBdc12_z + dEdc12_z);

		qeal dBTdc21_x = q1 * dCdc21_x + q2 * (dBdc21_x + dEdc21_x);
		qeal dBTdc21_y = q1 * dCdc21_y + q2 * (dBdc21_y + dEdc21_y);
		qeal dBTdc21_z = q1 * dCdc21_z + q2 * (dBdc21_z + dEdc21_z);

		qeal dBTdc22_x = q1 * dCdc22_x + q2 * (dBdc22_x + dEdc22_x);
		qeal dBTdc22_y = q1 * dCdc22_y + q2 * (dBdc22_y + dEdc22_y);
		qeal dBTdc22_z = q1 * dCdc22_z + q2 * (dBdc22_z + dEdc22_z);

		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (alpha * C1x + beta * C2x + C3x);
		qeal vy = 2.0 * (alpha * C1y + beta * C2y + C3y);
		qeal vz = 2.0 * (alpha * C1z + beta * C2z + C3z);

		qeal dvxdc11_x = 2.0 * (C1x * dALdc11_x + C2x * dBTdc11_x + alpha);
		qeal dvxdc11_y = 2.0 * (C1x * dALdc11_y + C2x * dBTdc11_y);
		qeal dvxdc11_z = 2.0 * (C1x * dALdc11_z + C2x * dBTdc11_z);

		qeal dvxdc12_x = 2.0 * (C1x * dALdc12_x + C2x * dBTdc12_x - alpha + 1);
		qeal dvxdc12_y = 2.0 * (C1x * dALdc12_y + C2x * dBTdc12_y);
		qeal dvxdc12_z = 2.0 * (C1x * dALdc12_z + C2x * dBTdc12_z);

		qeal dvxdc21_x = 2.0 * (C1x * dALdc21_x + C2x * dBTdc21_x - beta);
		qeal dvxdc21_y = 2.0 * (C1x * dALdc21_y + C2x * dBTdc21_y);
		qeal dvxdc21_z = 2.0 * (C1x * dALdc21_z + C2x * dBTdc21_z);

		qeal dvxdc22_x = 2.0 * (C1x * dALdc22_x + C2x * dBTdc22_x + beta - 1);
		qeal dvxdc22_y = 2.0 * (C1x * dALdc22_y + C2x * dBTdc22_y);
		qeal dvxdc22_z = 2.0 * (C1x * dALdc22_z + C2x * dBTdc22_z);
		//
		qeal dvydc11_x = 2.0 * (C1y * dALdc11_x + C2y * dBTdc11_x);
		qeal dvydc11_y = 2.0 * (C1y * dALdc11_y + C2y * dBTdc11_y + alpha);
		qeal dvydc11_z = 2.0 * (C1y * dALdc11_z + C2y * dBTdc11_z);

		qeal dvydc12_x = 2.0 * (C1y * dALdc12_x + C2y * dBTdc12_x);
		qeal dvydc12_y = 2.0 * (C1y * dALdc12_y + C2y * dBTdc12_y - alpha + 1);
		qeal dvydc12_z = 2.0 * (C1y * dALdc12_z + C2y * dBTdc12_z);

		qeal dvydc21_x = 2.0 * (C1y * dALdc21_x + C2y * dBTdc21_x);
		qeal dvydc21_y = 2.0 * (C1y * dALdc21_y + C2y * dBTdc21_y - beta);
		qeal dvydc21_z = 2.0 * (C1y * dALdc21_z + C2y * dBTdc21_z);

		qeal dvydc22_x = 2.0 * (C1y * dALdc22_x + C2y * dBTdc22_x);
		qeal dvydc22_y = 2.0 * (C1y * dALdc22_y + C2y * dBTdc22_y + beta - 1);
		qeal dvydc22_z = 2.0 * (C1y * dALdc22_z + C2y * dBTdc22_z);
		//
		qeal dvzdc11_x = 2.0 * (C1z * dALdc11_x + C2z * dBTdc11_x);
		qeal dvzdc11_y = 2.0 * (C1z * dALdc11_y + C2z * dBTdc11_y);
		qeal dvzdc11_z = 2.0 * (C1z * dALdc11_z + C2z * dBTdc11_z + alpha);

		qeal dvzdc12_x = 2.0 * (C1z * dALdc12_x + C2z * dBTdc12_x);
		qeal dvzdc12_y = 2.0 * (C1z * dALdc12_y + C2z * dBTdc12_y);
		qeal dvzdc12_z = 2.0 * (C1z * dALdc12_z + C2z * dBTdc12_z - alpha + 1);

		qeal dvzdc21_x = 2.0 * (C1z * dALdc21_x + C2z * dBTdc21_x);
		qeal dvzdc21_y = 2.0 * (C1z * dALdc21_y + C2z * dBTdc21_y);
		qeal dvzdc21_z = 2.0 * (C1z * dALdc21_z + C2z * dBTdc21_z - beta);

		qeal dvzdc22_x = 2.0 * (C1z * dALdc22_x + C2z * dBTdc22_x);
		qeal dvzdc22_y = 2.0 * (C1z * dALdc22_y + C2z * dBTdc22_y);
		qeal dvzdc22_z = 2.0 * (C1z * dALdc22_z + C2z * dBTdc22_z + beta - 1);

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = vx * dALdc11_x + alpha * dvxdc11_x;
		hessina.data()[1] = vy * dALdc11_x + alpha * dvydc11_x;
		hessina.data()[2] = vz * dALdc11_x + alpha * dvzdc11_x;

		hessina.data()[3] = vx * -dALdc11_x + (1 - alpha) * dvxdc11_x;
		hessina.data()[4] = vy * -dALdc11_x + (1 - alpha) * dvydc11_x;
		hessina.data()[5] = vz * -dALdc11_x + (1 - alpha) * dvzdc11_x;

		hessina.data()[6] = vx * -dBTdc11_x - beta * dvxdc11_x;
		hessina.data()[7] = vy * -dBTdc11_x - beta * dvydc11_x;
		hessina.data()[8] = vz * -dBTdc11_x - beta * dvzdc11_x;

		hessina.data()[9] = vx * dBTdc11_x + (beta - 1.0) * dvxdc11_x;
		hessina.data()[10] = vy * dBTdc11_x + (beta - 1.0) * dvydc11_x;
		hessina.data()[11] = vz * dBTdc11_x + (beta - 1.0) * dvzdc11_x;
		// f / c11y
		hessina.data()[12] = vx * dALdc11_y + alpha * dvxdc11_y;
		hessina.data()[13] = vy * dALdc11_y + alpha * dvydc11_y;
		hessina.data()[14] = vz * dALdc11_y + alpha * dvzdc11_y;

		hessina.data()[15] = vx * -dALdc11_y + (1 - alpha) * dvxdc11_y;
		hessina.data()[16] = vy * -dALdc11_y + (1 - alpha) * dvydc11_y;
		hessina.data()[17] = vz * -dALdc11_y + (1 - alpha) * dvzdc11_y;

		hessina.data()[18] = vx * -dBTdc11_y - beta * dvxdc11_y;
		hessina.data()[19] = vy * -dBTdc11_y - beta * dvydc11_y;
		hessina.data()[20] = vz * -dBTdc11_y - beta * dvzdc11_y;

		hessina.data()[21] = vx * dBTdc11_y + (beta - 1.0) * dvxdc11_y;
		hessina.data()[22] = vy * dBTdc11_y + (beta - 1.0) * dvydc11_y;
		hessina.data()[23] = vz * dBTdc11_y + (beta - 1.0) * dvzdc11_y;
		// f / c11z
		hessina.data()[24] = vx * dALdc11_z + alpha * dvxdc11_z;
		hessina.data()[25] = vy * dALdc11_z + alpha * dvydc11_z;
		hessina.data()[26] = vz * dALdc11_z + alpha * dvzdc11_z;

		hessina.data()[27] = vx * -dALdc11_z + (1 - alpha) * dvxdc11_z;
		hessina.data()[28] = vy * -dALdc11_z + (1 - alpha) * dvydc11_z;
		hessina.data()[29] = vz * -dALdc11_z + (1 - alpha) * dvzdc11_z;

		hessina.data()[30] = vx * -dBTdc11_z - beta * dvxdc11_z;
		hessina.data()[31] = vy * -dBTdc11_z - beta * dvydc11_z;
		hessina.data()[32] = vz * -dBTdc11_z - beta * dvzdc11_z;

		hessina.data()[33] = vx * dBTdc11_z + (beta - 1.0) * dvxdc11_z;
		hessina.data()[34] = vy * dBTdc11_z + (beta - 1.0) * dvydc11_z;
		hessina.data()[35] = vz * dBTdc11_z + (beta - 1.0) * dvzdc11_z;

		// f / c12x
		hessina.data()[36] = vx * dALdc12_x + alpha * dvxdc12_x;
		hessina.data()[37] = vy * dALdc12_x + alpha * dvydc12_x;
		hessina.data()[38] = vz * dALdc12_x + alpha * dvzdc12_x;

		hessina.data()[39] = vx * -dALdc12_x + (1 - alpha) * dvxdc12_x;
		hessina.data()[40] = vy * -dALdc12_x + (1 - alpha) * dvydc12_x;
		hessina.data()[41] = vz * -dALdc12_x + (1 - alpha) * dvzdc12_x;

		hessina.data()[42] = vx * -dBTdc12_x - beta * dvxdc12_x;
		hessina.data()[43] = vy * -dBTdc12_x - beta * dvydc12_x;
		hessina.data()[44] = vz * -dBTdc12_x - beta * dvzdc12_x;

		hessina.data()[45] = vx * dBTdc12_x + (beta - 1.0) * dvxdc12_x;
		hessina.data()[46] = vy * dBTdc12_x + (beta - 1.0) * dvydc12_x;
		hessina.data()[47] = vz * dBTdc12_x + (beta - 1.0) * dvzdc12_x;
		// f / c12y
		hessina.data()[48] = vx * dALdc12_y + alpha * dvxdc12_y;
		hessina.data()[49] = vy * dALdc12_y + alpha * dvydc12_y;
		hessina.data()[50] = vz * dALdc12_y + alpha * dvzdc12_y;

		hessina.data()[51] = vx * -dALdc12_y + (1 - alpha) * dvxdc12_y;
		hessina.data()[52] = vy * -dALdc12_y + (1 - alpha) * dvydc12_y;
		hessina.data()[53] = vz * -dALdc12_y + (1 - alpha) * dvzdc12_y;

		hessina.data()[54] = vx * -dBTdc12_y - beta * dvxdc12_y;
		hessina.data()[55] = vy * -dBTdc12_y - beta * dvydc12_y;
		hessina.data()[56] = vz * -dBTdc12_y - beta * dvzdc12_y;

		hessina.data()[57] = vx * dBTdc12_y + (beta - 1.0) * dvxdc12_y;
		hessina.data()[58] = vy * dBTdc12_y + (beta - 1.0) * dvydc12_y;
		hessina.data()[59] = vz * dBTdc12_y + (beta - 1.0) * dvzdc12_y;
		// f / c12z
		hessina.data()[60] = vx * dALdc12_z + alpha * dvxdc12_z;
		hessina.data()[61] = vy * dALdc12_z + alpha * dvydc12_z;
		hessina.data()[62] = vz * dALdc12_z + alpha * dvzdc12_z;

		hessina.data()[63] = vx * -dALdc12_z + (1 - alpha) * dvxdc12_z;
		hessina.data()[64] = vy * -dALdc12_z + (1 - alpha) * dvydc12_z;
		hessina.data()[65] = vz * -dALdc12_z + (1 - alpha) * dvzdc12_z;

		hessina.data()[66] = vx * -dBTdc12_z - beta * dvxdc12_z;
		hessina.data()[67] = vy * -dBTdc12_z - beta * dvydc12_z;
		hessina.data()[68] = vz * -dBTdc12_z - beta * dvzdc12_z;

		hessina.data()[69] = vx * dBTdc12_z + (beta - 1.0) * dvxdc12_z;
		hessina.data()[70] = vy * dBTdc12_z + (beta - 1.0) * dvydc12_z;
		hessina.data()[71] = vz * dBTdc12_z + (beta - 1.0) * dvzdc12_z;

		// f / c21x
		hessina.data()[72] = vx * dALdc21_x + alpha * dvxdc21_x;
		hessina.data()[73] = vy * dALdc21_x + alpha * dvydc21_x;
		hessina.data()[74] = vz * dALdc21_x + alpha * dvzdc21_x;

		hessina.data()[75] = vx * -dALdc21_x + (1 - alpha) * dvxdc21_x;
		hessina.data()[76] = vy * -dALdc21_x + (1 - alpha) * dvydc21_x;
		hessina.data()[77] = vz * -dALdc21_x + (1 - alpha) * dvzdc21_x;

		hessina.data()[78] = vx * -dBTdc21_x - beta * dvxdc21_x;
		hessina.data()[79] = vy * -dBTdc21_x - beta * dvydc21_x;
		hessina.data()[80] = vz * -dBTdc21_x - beta * dvzdc21_x;

		hessina.data()[81] = vx * dBTdc21_x + (beta - 1.0) * dvxdc21_x;
		hessina.data()[82] = vy * dBTdc21_x + (beta - 1.0) * dvydc21_x;
		hessina.data()[83] = vz * dBTdc21_x + (beta - 1.0) * dvzdc21_x;
		// f / c21y
		hessina.data()[84] = vx * dALdc21_y + alpha * dvxdc21_y;
		hessina.data()[85] = vy * dALdc21_y + alpha * dvydc21_y;
		hessina.data()[86] = vz * dALdc21_y + alpha * dvzdc21_y;

		hessina.data()[87] = vx * -dALdc21_y + (1 - alpha) * dvxdc21_y;
		hessina.data()[88] = vy * -dALdc21_y + (1 - alpha) * dvydc21_y;
		hessina.data()[89] = vz * -dALdc21_y + (1 - alpha) * dvzdc21_y;

		hessina.data()[90] = vx * -dBTdc21_y - beta * dvxdc21_y;
		hessina.data()[91] = vy * -dBTdc21_y - beta * dvydc21_y;
		hessina.data()[92] = vz * -dBTdc21_y - beta * dvzdc21_y;

		hessina.data()[93] = vx * dBTdc21_y + (beta - 1.0) * dvxdc21_y;
		hessina.data()[94] = vy * dBTdc21_y + (beta - 1.0) * dvydc21_y;
		hessina.data()[95] = vz * dBTdc21_y + (beta - 1.0) * dvzdc21_y;
		// f / c21z
		hessina.data()[96] = vx * dALdc21_z + alpha * dvxdc21_z;
		hessina.data()[97] = vy * dALdc21_z + alpha * dvydc21_z;
		hessina.data()[98] = vz * dALdc21_z + alpha * dvzdc21_z;

		hessina.data()[99] = vx * -dALdc21_z + (1 - alpha) * dvxdc21_z;
		hessina.data()[100] = vy * -dALdc21_z + (1 - alpha) * dvydc21_z;
		hessina.data()[101] = vz * -dALdc21_z + (1 - alpha) * dvzdc21_z;

		hessina.data()[102] = vx * -dBTdc21_z - beta * dvxdc21_z;
		hessina.data()[103] = vy * -dBTdc21_z - beta * dvydc21_z;
		hessina.data()[104] = vz * -dBTdc21_z - beta * dvzdc21_z;

		hessina.data()[105] = vx * dBTdc21_z + (beta - 1.0) * dvxdc21_z;
		hessina.data()[106] = vy * dBTdc21_z + (beta - 1.0) * dvydc21_z;
		hessina.data()[107] = vz * dBTdc21_z + (beta - 1.0) * dvzdc21_z;

		// f / c22x
		hessina.data()[108] = vx * dALdc22_x + alpha * dvxdc22_x;
		hessina.data()[109] = vy * dALdc22_x + alpha * dvydc22_x;
		hessina.data()[110] = vz * dALdc22_x + alpha * dvzdc22_x;

		hessina.data()[111] = vx * -dALdc22_x + (1 - alpha) * dvxdc22_x;
		hessina.data()[112] = vy * -dALdc22_x + (1 - alpha) * dvydc22_x;
		hessina.data()[113] = vz * -dALdc22_x + (1 - alpha) * dvzdc22_x;

		hessina.data()[114] = vx * -dBTdc22_x - beta * dvxdc22_x;
		hessina.data()[115] = vy * -dBTdc22_x - beta * dvydc22_x;
		hessina.data()[116] = vz * -dBTdc22_x - beta * dvzdc22_x;

		hessina.data()[117] = vx * dBTdc22_x + (beta - 1.0) * dvxdc22_x;
		hessina.data()[118] = vy * dBTdc22_x + (beta - 1.0) * dvydc22_x;
		hessina.data()[119] = vz * dBTdc22_x + (beta - 1.0) * dvzdc22_x;
		// f / c22y
		hessina.data()[120] = vx * dALdc22_y + alpha * dvxdc22_y;
		hessina.data()[121] = vy * dALdc22_y + alpha * dvydc22_y;
		hessina.data()[122] = vz * dALdc22_y + alpha * dvzdc22_y;

		hessina.data()[123] = vx * -dALdc22_y + (1 - alpha) * dvxdc22_y;
		hessina.data()[124] = vy * -dALdc22_y + (1 - alpha) * dvydc22_y;
		hessina.data()[125] = vz * -dALdc22_y + (1 - alpha) * dvzdc22_y;

		hessina.data()[126] = vx * -dBTdc22_y - beta * dvxdc22_y;
		hessina.data()[127] = vy * -dBTdc22_y - beta * dvydc22_y;
		hessina.data()[128] = vz * -dBTdc22_y - beta * dvzdc22_y;

		hessina.data()[129] = vx * dBTdc22_y + (beta - 1.0) * dvxdc22_y;
		hessina.data()[130] = vy * dBTdc22_y + (beta - 1.0) * dvydc22_y;
		hessina.data()[131] = vz * dBTdc22_y + (beta - 1.0) * dvzdc22_y;
		// f / c22z
		hessina.data()[132] = vx * dALdc22_z + alpha * dvxdc22_z;
		hessina.data()[133] = vy * dALdc22_z + alpha * dvydc22_z;
		hessina.data()[134] = vz * dALdc22_z + alpha * dvzdc22_z;

		hessina.data()[135] = vx * -dALdc22_z + (1 - alpha) * dvxdc22_z;
		hessina.data()[136] = vy * -dALdc22_z + (1 - alpha) * dvydc22_z;
		hessina.data()[137] = vz * -dALdc22_z + (1 - alpha) * dvzdc22_z;

		hessina.data()[138] = vx * -dBTdc22_z - beta * dvxdc22_z;
		hessina.data()[139] = vy * -dBTdc22_z - beta * dvydc22_z;
		hessina.data()[140] = vz * -dBTdc22_z - beta * dvzdc22_z;

		hessina.data()[141] = vx * dBTdc22_z + (beta - 1.0) * dvxdc22_z;
		hessina.data()[142] = vy * dBTdc22_z + (beta - 1.0) * dvydc22_z;
		hessina.data()[143] = vz * dBTdc22_z + (beta - 1.0) * dvzdc22_z;
	}

	void MipcConeConeConstraint::betaIsZeroHessina(MatrixX & hessina)
	{
		qeal dAdc11_x = 2.0 * sC1.data()[0]; qeal dAdc11_y = 2.0 * sC1.data()[1]; qeal dAdc11_z = 2.0 * sC1.data()[2];
		qeal dAdc12_x = -2.0 * sC1.data()[0]; qeal dAdc12_y = -2.0 * sC1.data()[1]; qeal dAdc12_z = -2.0 * sC1.data()[2];

		qeal dBdc11_x = 2.0 * sC2.data()[0]; qeal dBdc11_y = 2.0 * sC2.data()[1]; qeal dBdc11_z = 2.0 * sC2.data()[2];
		qeal dBdc12_x = -2.0 * sC2.data()[0]; qeal dBdc12_y = -2.0 * sC2.data()[1]; qeal dBdc12_z = -2.0 * sC2.data()[2];
		qeal dBdc21_x = -2.0 * sC1.data()[0]; qeal dBdc21_y = -2.0 * sC1.data()[1]; qeal dBdc21_z = -2.0 * sC1.data()[2];
		qeal dBdc22_x = 2.0 * sC1.data()[0]; qeal dBdc22_y = 2.0 * sC1.data()[1]; qeal dBdc22_z = 2.0 * sC1.data()[2];

		qeal dCdc21_x = -2.0 * sC2.data()[0]; qeal dCdc21_y = -2.0 * sC2.data()[1]; qeal dCdc21_z = -2.0 * sC2.data()[2];
		qeal dCdc22_x = 2.0 * sC2.data()[0]; qeal dCdc22_y = 2.0 * sC2.data()[1]; qeal dCdc22_z = 2.0 * sC2.data()[2];

		qeal dDdc11_x = 2.0 * sC3.data()[0]; qeal dDdc11_y = 2.0 * sC3.data()[1]; qeal dDdc11_z = 2.0 * sC3.data()[2];
		qeal dDdc12_x = 2.0 * (sC1.data()[0] - sC3.data()[0]); qeal dDdc12_y = 2.0 * (sC1.data()[1] - sC3.data()[1]); qeal dDdc12_z = 2.0 * (sC1.data()[2] - sC3.data()[2]);
		qeal dDdc22_x = -2.0 * sC1.data()[0]; qeal dDdc22_y = -2.0 * sC1.data()[1]; qeal dDdc22_z = -2.0 * sC1.data()[2];

		qeal dEdc12_x = 2.0 * sC2.data()[0]; qeal dEdc12_y = 2.0 * sC2.data()[1]; qeal dEdc12_z = 2.0 * sC2.data()[2];
		qeal dEdc21_x = -2.0 * sC3.data()[0]; qeal dEdc21_y = -2.0 * sC3.data()[1]; qeal dEdc21_z = -2.0 * sC3.data()[2];
		qeal dEdc22_x = 2.0 * (sC3.data()[0] - sC2.data()[0]); qeal dEdc22_y = 2.0 * (sC3.data()[1] - sC2.data()[1]); qeal dEdc22_z = 2.0 * (sC3.data()[2] - sC2.data()[2]);

		qeal dFdc12_x = 2.0 * sC3.data()[0]; qeal dFdc12_y = 2.0 * sC3.data()[1]; qeal dFdc12_z = 2.0 * sC3.data()[2];
		qeal dFdc22 = -2.0 * sC3.data()[0]; qeal dFdc22_y = -2.0 * sC3.data()[1]; qeal dFd22_z = -2.0 * sC3.data()[2];

		//
		qeal q1 = 0.5 * D / (A * A);  qeal q2 = -0.5 / A;

		qeal dALdc11_x = q1 * dAdc11_x + q2 * dDdc11_x;
		qeal dALdc11_y = q1 * dAdc11_y + q2 * dDdc11_y;
		qeal dALdc11_z = q1 * dAdc11_z + q2 * dDdc11_z;

		qeal dALdc12_x = q1 * dAdc12_x + q2 * dDdc12_x;
		qeal dALdc12_y = q1 * dAdc12_y + q2 * dDdc12_y;
		qeal dALdc12_z = q1 * dAdc12_z + q2 * dDdc12_z;

		qeal dALdc21_x = 0;
		qeal dALdc21_y = 0;
		qeal dALdc21_z = 0;

		qeal dALdc22_x = q2 * dDdc22_x;
		qeal dALdc22_y = q2 * dDdc22_y;
		qeal dALdc22_z = q2 * dDdc22_z;

		qeal dBTdc11_x = 0;
		qeal dBTdc11_y = 0;
		qeal dBTdc11_z = 0;

		qeal dBTdc12_x = 0;
		qeal dBTdc12_y = 0;
		qeal dBTdc12_z = 0;

		qeal dBTdc21_x = 0;
		qeal dBTdc21_y = 0;
		qeal dBTdc21_z = 0;

		qeal dBTdc22_x = 0;
		qeal dBTdc22_y = 0;
		qeal dBTdc22_z = 0;

		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (alpha * C1x + C3x);
		qeal vy = 2.0 * (alpha * C1y + C3y);
		qeal vz = 2.0 * (alpha * C1z + C3z);

		qeal dvxdc11_x = 2.0 * (C1x * dALdc11_x + C2x * dBTdc11_x + alpha);
		qeal dvxdc11_y = 2.0 * (C1x * dALdc11_y + C2x * dBTdc11_y);
		qeal dvxdc11_z = 2.0 * (C1x * dALdc11_z + C2x * dBTdc11_z);

		qeal dvxdc12_x = 2.0 * (C1x * dALdc12_x + C2x * dBTdc12_x - alpha + 1);
		qeal dvxdc12_y = 2.0 * (C1x * dALdc12_y + C2x * dBTdc12_y);
		qeal dvxdc12_z = 2.0 * (C1x * dALdc12_z + C2x * dBTdc12_z);

		qeal dvxdc21_x = 0;
		qeal dvxdc21_y = 0;
		qeal dvxdc21_z = 0;

		qeal dvxdc22_x = 2.0 * (C1x * dALdc22_x + C2x * dBTdc22_x + beta - 1);
		qeal dvxdc22_y = 2.0 * (C1x * dALdc22_y + C2x * dBTdc22_y);
		qeal dvxdc22_z = 2.0 * (C1x * dALdc22_z + C2x * dBTdc22_z);
		//
		qeal dvydc11_x = 2.0 * (C1y * dALdc11_x + C2y * dBTdc11_x);
		qeal dvydc11_y = 2.0 * (C1y * dALdc11_y + C2y * dBTdc11_y + alpha);
		qeal dvydc11_z = 2.0 * (C1y * dALdc11_z + C2y * dBTdc11_z);

		qeal dvydc12_x = 2.0 * (C1y * dALdc12_x + C2y * dBTdc12_x);
		qeal dvydc12_y = 2.0 * (C1y * dALdc12_y + C2y * dBTdc12_y - alpha + 1);
		qeal dvydc12_z = 2.0 * (C1y * dALdc12_z + C2y * dBTdc12_z);

		qeal dvydc21_x = 0;
		qeal dvydc21_y = 0;
		qeal dvydc21_z = 0;

		qeal dvydc22_x = 2.0 * (C1y * dALdc22_x + C2y * dBTdc22_x);
		qeal dvydc22_y = 2.0 * (C1y * dALdc22_y + C2y * dBTdc22_y + beta - 1);
		qeal dvydc22_z = 2.0 * (C1y * dALdc22_z + C2y * dBTdc22_z);
		//
		qeal dvzdc11_x = 2.0 * (C1z * dALdc11_x + C2z * dBTdc11_x);
		qeal dvzdc11_y = 2.0 * (C1z * dALdc11_y + C2z * dBTdc11_y);
		qeal dvzdc11_z = 2.0 * (C1z * dALdc11_z + C2z * dBTdc11_z + alpha);

		qeal dvzdc12_x = 2.0 * (C1z * dALdc12_x + C2z * dBTdc12_x);
		qeal dvzdc12_y = 2.0 * (C1z * dALdc12_y + C2z * dBTdc12_y);
		qeal dvzdc12_z = 2.0 * (C1z * dALdc12_z + C2z * dBTdc12_z - alpha + 1);

		qeal dvzdc21_x = 0;
		qeal dvzdc21_y = 0;
		qeal dvzdc21_z = 0;

		qeal dvzdc22_x = 2.0 * (C1z * dALdc22_x + C2z * dBTdc22_x);
		qeal dvzdc22_y = 2.0 * (C1z * dALdc22_y + C2z * dBTdc22_y);
		qeal dvzdc22_z = 2.0 * (C1z * dALdc22_z + C2z * dBTdc22_z + beta - 1);

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = vx * dALdc11_x + alpha * dvxdc11_x;
		hessina.data()[1] = vy * dALdc11_x + alpha * dvydc11_x;
		hessina.data()[2] = vz * dALdc11_x + alpha * dvzdc11_x;

		hessina.data()[3] = vx * -dALdc11_x + (1 - alpha) * dvxdc11_x;
		hessina.data()[4] = vy * -dALdc11_x + (1 - alpha) * dvydc11_x;
		hessina.data()[5] = vz * -dALdc11_x + (1 - alpha) * dvzdc11_x;

		hessina.data()[6] = vx * -dBTdc11_x - beta * dvxdc11_x;
		hessina.data()[7] = vy * -dBTdc11_x - beta * dvydc11_x;
		hessina.data()[8] = vz * -dBTdc11_x - beta * dvzdc11_x;

		hessina.data()[9] = vx * dBTdc11_x + (beta - 1.0) * dvxdc11_x;
		hessina.data()[10] = vy * dBTdc11_x + (beta - 1.0) * dvydc11_x;
		hessina.data()[11] = vz * dBTdc11_x + (beta - 1.0) * dvzdc11_x;
		// f / c11y
		hessina.data()[12] = vx * dALdc11_y + alpha * dvxdc11_y;
		hessina.data()[13] = vy * dALdc11_y + alpha * dvydc11_y;
		hessina.data()[14] = vz * dALdc11_y + alpha * dvzdc11_y;

		hessina.data()[15] = vx * -dALdc11_y + (1 - alpha) * dvxdc11_y;
		hessina.data()[16] = vy * -dALdc11_y + (1 - alpha) * dvydc11_y;
		hessina.data()[17] = vz * -dALdc11_y + (1 - alpha) * dvzdc11_y;

		hessina.data()[18] = vx * -dBTdc11_y - beta * dvxdc11_y;
		hessina.data()[19] = vy * -dBTdc11_y - beta * dvydc11_y;
		hessina.data()[20] = vz * -dBTdc11_y - beta * dvzdc11_y;

		hessina.data()[21] = vx * dBTdc11_y + (beta - 1.0) * dvxdc11_y;
		hessina.data()[22] = vy * dBTdc11_y + (beta - 1.0) * dvydc11_y;
		hessina.data()[23] = vz * dBTdc11_y + (beta - 1.0) * dvzdc11_y;
		// f / c11z
		hessina.data()[24] = vx * dALdc11_z + alpha * dvxdc11_z;
		hessina.data()[25] = vy * dALdc11_z + alpha * dvydc11_z;
		hessina.data()[26] = vz * dALdc11_z + alpha * dvzdc11_z;

		hessina.data()[27] = vx * -dALdc11_z + (1 - alpha) * dvxdc11_z;
		hessina.data()[28] = vy * -dALdc11_z + (1 - alpha) * dvydc11_z;
		hessina.data()[29] = vz * -dALdc11_z + (1 - alpha) * dvzdc11_z;

		hessina.data()[30] = vx * -dBTdc11_z - beta * dvxdc11_z;
		hessina.data()[31] = vy * -dBTdc11_z - beta * dvydc11_z;
		hessina.data()[32] = vz * -dBTdc11_z - beta * dvzdc11_z;

		hessina.data()[33] = vx * dBTdc11_z + (beta - 1.0) * dvxdc11_z;
		hessina.data()[34] = vy * dBTdc11_z + (beta - 1.0) * dvydc11_z;
		hessina.data()[35] = vz * dBTdc11_z + (beta - 1.0) * dvzdc11_z;

		// f / c12x
		hessina.data()[36] = vx * dALdc12_x + alpha * dvxdc12_x;
		hessina.data()[37] = vy * dALdc12_x + alpha * dvydc12_x;
		hessina.data()[38] = vz * dALdc12_x + alpha * dvzdc12_x;

		hessina.data()[39] = vx * -dALdc12_x + (1 - alpha) * dvxdc12_x;
		hessina.data()[40] = vy * -dALdc12_x + (1 - alpha) * dvydc12_x;
		hessina.data()[41] = vz * -dALdc12_x + (1 - alpha) * dvzdc12_x;

		hessina.data()[42] = vx * -dBTdc12_x - beta * dvxdc12_x;
		hessina.data()[43] = vy * -dBTdc12_x - beta * dvydc12_x;
		hessina.data()[44] = vz * -dBTdc12_x - beta * dvzdc12_x;

		hessina.data()[45] = vx * dBTdc12_x + (beta - 1.0) * dvxdc12_x;
		hessina.data()[46] = vy * dBTdc12_x + (beta - 1.0) * dvydc12_x;
		hessina.data()[47] = vz * dBTdc12_x + (beta - 1.0) * dvzdc12_x;
		// f / c12y
		hessina.data()[48] = vx * dALdc12_y + alpha * dvxdc12_y;
		hessina.data()[49] = vy * dALdc12_y + alpha * dvydc12_y;
		hessina.data()[50] = vz * dALdc12_y + alpha * dvzdc12_y;

		hessina.data()[51] = vx * -dALdc12_y + (1 - alpha) * dvxdc12_y;
		hessina.data()[52] = vy * -dALdc12_y + (1 - alpha) * dvydc12_y;
		hessina.data()[53] = vz * -dALdc12_y + (1 - alpha) * dvzdc12_y;

		hessina.data()[54] = vx * -dBTdc12_y - beta * dvxdc12_y;
		hessina.data()[55] = vy * -dBTdc12_y - beta * dvydc12_y;
		hessina.data()[56] = vz * -dBTdc12_y - beta * dvzdc12_y;

		hessina.data()[57] = vx * dBTdc12_y + (beta - 1.0) * dvxdc12_y;
		hessina.data()[58] = vy * dBTdc12_y + (beta - 1.0) * dvydc12_y;
		hessina.data()[59] = vz * dBTdc12_y + (beta - 1.0) * dvzdc12_y;
		// f / c12z
		hessina.data()[60] = vx * dALdc12_z + alpha * dvxdc12_z;
		hessina.data()[61] = vy * dALdc12_z + alpha * dvydc12_z;
		hessina.data()[62] = vz * dALdc12_z + alpha * dvzdc12_z;

		hessina.data()[63] = vx * -dALdc12_z + (1 - alpha) * dvxdc12_z;
		hessina.data()[64] = vy * -dALdc12_z + (1 - alpha) * dvydc12_z;
		hessina.data()[65] = vz * -dALdc12_z + (1 - alpha) * dvzdc12_z;

		hessina.data()[66] = vx * -dBTdc12_z - beta * dvxdc12_z;
		hessina.data()[67] = vy * -dBTdc12_z - beta * dvydc12_z;
		hessina.data()[68] = vz * -dBTdc12_z - beta * dvzdc12_z;

		hessina.data()[69] = vx * dBTdc12_z + (beta - 1.0) * dvxdc12_z;
		hessina.data()[70] = vy * dBTdc12_z + (beta - 1.0) * dvydc12_z;
		hessina.data()[71] = vz * dBTdc12_z + (beta - 1.0) * dvzdc12_z;

		// f / c21x
		hessina.data()[72] = vx * dALdc21_x + alpha * dvxdc21_x;
		hessina.data()[73] = vy * dALdc21_x + alpha * dvydc21_x;
		hessina.data()[74] = vz * dALdc21_x + alpha * dvzdc21_x;

		hessina.data()[75] = vx * -dALdc21_x + (1 - alpha) * dvxdc21_x;
		hessina.data()[76] = vy * -dALdc21_x + (1 - alpha) * dvydc21_x;
		hessina.data()[77] = vz * -dALdc21_x + (1 - alpha) * dvzdc21_x;

		hessina.data()[78] = vx * -dBTdc21_x - beta * dvxdc21_x;
		hessina.data()[79] = vy * -dBTdc21_x - beta * dvydc21_x;
		hessina.data()[80] = vz * -dBTdc21_x - beta * dvzdc21_x;

		hessina.data()[81] = vx * dBTdc21_x + (beta - 1.0) * dvxdc21_x;
		hessina.data()[82] = vy * dBTdc21_x + (beta - 1.0) * dvydc21_x;
		hessina.data()[83] = vz * dBTdc21_x + (beta - 1.0) * dvzdc21_x;
		// f / c21y
		hessina.data()[84] = vx * dALdc21_y + alpha * dvxdc21_y;
		hessina.data()[85] = vy * dALdc21_y + alpha * dvydc21_y;
		hessina.data()[86] = vz * dALdc21_y + alpha * dvzdc21_y;

		hessina.data()[87] = vx * -dALdc21_y + (1 - alpha) * dvxdc21_y;
		hessina.data()[88] = vy * -dALdc21_y + (1 - alpha) * dvydc21_y;
		hessina.data()[89] = vz * -dALdc21_y + (1 - alpha) * dvzdc21_y;

		hessina.data()[90] = vx * -dBTdc21_y - beta * dvxdc21_y;
		hessina.data()[91] = vy * -dBTdc21_y - beta * dvydc21_y;
		hessina.data()[92] = vz * -dBTdc21_y - beta * dvzdc21_y;

		hessina.data()[93] = vx * dBTdc21_y + (beta - 1.0) * dvxdc21_y;
		hessina.data()[94] = vy * dBTdc21_y + (beta - 1.0) * dvydc21_y;
		hessina.data()[95] = vz * dBTdc21_y + (beta - 1.0) * dvzdc21_y;
		// f / c21z
		hessina.data()[96] = vx * dALdc21_z + alpha * dvxdc21_z;
		hessina.data()[97] = vy * dALdc21_z + alpha * dvydc21_z;
		hessina.data()[98] = vz * dALdc21_z + alpha * dvzdc21_z;

		hessina.data()[99] = vx * -dALdc21_z + (1 - alpha) * dvxdc21_z;
		hessina.data()[100] = vy * -dALdc21_z + (1 - alpha) * dvydc21_z;
		hessina.data()[101] = vz * -dALdc21_z + (1 - alpha) * dvzdc21_z;

		hessina.data()[102] = vx * -dBTdc21_z - beta * dvxdc21_z;
		hessina.data()[103] = vy * -dBTdc21_z - beta * dvydc21_z;
		hessina.data()[104] = vz * -dBTdc21_z - beta * dvzdc21_z;

		hessina.data()[105] = vx * dBTdc21_z + (beta - 1.0) * dvxdc21_z;
		hessina.data()[106] = vy * dBTdc21_z + (beta - 1.0) * dvydc21_z;
		hessina.data()[107] = vz * dBTdc21_z + (beta - 1.0) * dvzdc21_z;

		// f / c22x
		hessina.data()[108] = vx * dALdc22_x + alpha * dvxdc22_x;
		hessina.data()[109] = vy * dALdc22_x + alpha * dvydc22_x;
		hessina.data()[110] = vz * dALdc22_x + alpha * dvzdc22_x;

		hessina.data()[111] = vx * -dALdc22_x + (1 - alpha) * dvxdc22_x;
		hessina.data()[112] = vy * -dALdc22_x + (1 - alpha) * dvydc22_x;
		hessina.data()[113] = vz * -dALdc22_x + (1 - alpha) * dvzdc22_x;

		hessina.data()[114] = vx * -dBTdc22_x - beta * dvxdc22_x;
		hessina.data()[115] = vy * -dBTdc22_x - beta * dvydc22_x;
		hessina.data()[116] = vz * -dBTdc22_x - beta * dvzdc22_x;

		hessina.data()[117] = vx * dBTdc22_x + (beta - 1.0) * dvxdc22_x;
		hessina.data()[118] = vy * dBTdc22_x + (beta - 1.0) * dvydc22_x;
		hessina.data()[119] = vz * dBTdc22_x + (beta - 1.0) * dvzdc22_x;
		// f / c22y
		hessina.data()[120] = vx * dALdc22_y + alpha * dvxdc22_y;
		hessina.data()[121] = vy * dALdc22_y + alpha * dvydc22_y;
		hessina.data()[122] = vz * dALdc22_y + alpha * dvzdc22_y;

		hessina.data()[123] = vx * -dALdc22_y + (1 - alpha) * dvxdc22_y;
		hessina.data()[124] = vy * -dALdc22_y + (1 - alpha) * dvydc22_y;
		hessina.data()[125] = vz * -dALdc22_y + (1 - alpha) * dvzdc22_y;

		hessina.data()[126] = vx * -dBTdc22_y - beta * dvxdc22_y;
		hessina.data()[127] = vy * -dBTdc22_y - beta * dvydc22_y;
		hessina.data()[128] = vz * -dBTdc22_y - beta * dvzdc22_y;

		hessina.data()[129] = vx * dBTdc22_y + (beta - 1.0) * dvxdc22_y;
		hessina.data()[130] = vy * dBTdc22_y + (beta - 1.0) * dvydc22_y;
		hessina.data()[131] = vz * dBTdc22_y + (beta - 1.0) * dvzdc22_y;
		// f / c22z
		hessina.data()[132] = vx * dALdc22_z + alpha * dvxdc22_z;
		hessina.data()[133] = vy * dALdc22_z + alpha * dvydc22_z;
		hessina.data()[134] = vz * dALdc22_z + alpha * dvzdc22_z;

		hessina.data()[135] = vx * -dALdc22_z + (1 - alpha) * dvxdc22_z;
		hessina.data()[136] = vy * -dALdc22_z + (1 - alpha) * dvydc22_z;
		hessina.data()[137] = vz * -dALdc22_z + (1 - alpha) * dvzdc22_z;

		hessina.data()[138] = vx * -dBTdc22_z - beta * dvxdc22_z;
		hessina.data()[139] = vy * -dBTdc22_z - beta * dvydc22_z;
		hessina.data()[140] = vz * -dBTdc22_z - beta * dvzdc22_z;

		hessina.data()[141] = vx * dBTdc22_z + (beta - 1.0) * dvxdc22_z;
		hessina.data()[142] = vy * dBTdc22_z + (beta - 1.0) * dvydc22_z;
		hessina.data()[143] = vz * dBTdc22_z + (beta - 1.0) * dvzdc22_z;
	}

	void MipcConeConeConstraint::betaIsOneHessina(MatrixX & hessina)
	{
		qeal dAdc11_x = 2.0 * sC1.data()[0]; qeal dAdc11_y = 2.0 * sC1.data()[1]; qeal dAdc11_z = 2.0 * sC1.data()[2];
		qeal dAdc12_x = -2.0 * sC1.data()[0]; qeal dAdc12_y = -2.0 * sC1.data()[1]; qeal dAdc12_z = -2.0 * sC1.data()[2];

		qeal dBdc11_x = 2.0 * sC2.data()[0]; qeal dBdc11_y = 2.0 * sC2.data()[1]; qeal dBdc11_z = 2.0 * sC2.data()[2];
		qeal dBdc12_x = -2.0 * sC2.data()[0]; qeal dBdc12_y = -2.0 * sC2.data()[1]; qeal dBdc12_z = -2.0 * sC2.data()[2];
		qeal dBdc21_x = -2.0 * sC1.data()[0]; qeal dBdc21_y = -2.0 * sC1.data()[1]; qeal dBdc21_z = -2.0 * sC1.data()[2];
		qeal dBdc22_x = 2.0 * sC1.data()[0]; qeal dBdc22_y = 2.0 * sC1.data()[1]; qeal dBdc22_z = 2.0 * sC1.data()[2];

		qeal dCdc21_x = -2.0 * sC2.data()[0]; qeal dCdc21_y = -2.0 * sC2.data()[1]; qeal dCdc21_z = -2.0 * sC2.data()[2];
		qeal dCdc22_x = 2.0 * sC2.data()[0]; qeal dCdc22_y = 2.0 * sC2.data()[1]; qeal dCdc22_z = 2.0 * sC2.data()[2];

		qeal dDdc11_x = 2.0 * sC3.data()[0]; qeal dDdc11_y = 2.0 * sC3.data()[1]; qeal dDdc11_z = 2.0 * sC3.data()[2];
		qeal dDdc12_x = 2.0 * (sC1.data()[0] - sC3.data()[0]); qeal dDdc12_y = 2.0 * (sC1.data()[1] - sC3.data()[1]); qeal dDdc12_z = 2.0 * (sC1.data()[2] - sC3.data()[2]);
		qeal dDdc22_x = -2.0 * sC1.data()[0]; qeal dDdc22_y = -2.0 * sC1.data()[1]; qeal dDdc22_z = -2.0 * sC1.data()[2];

		qeal dEdc12_x = 2.0 * sC2.data()[0]; qeal dEdc12_y = 2.0 * sC2.data()[1]; qeal dEdc12_z = 2.0 * sC2.data()[2];
		qeal dEdc21_x = -2.0 * sC3.data()[0]; qeal dEdc21_y = -2.0 * sC3.data()[1]; qeal dEdc21_z = -2.0 * sC3.data()[2];
		qeal dEdc22_x = 2.0 * (sC3.data()[0] - sC2.data()[0]); qeal dEdc22_y = 2.0 * (sC3.data()[1] - sC2.data()[1]); qeal dEdc22_z = 2.0 * (sC3.data()[2] - sC2.data()[2]);

		qeal dFdc12_x = 2.0 * sC3.data()[0]; qeal dFdc12_y = 2.0 * sC3.data()[1]; qeal dFdc12_z = 2.0 * sC3.data()[2];
		qeal dFdc22 = -2.0 * sC3.data()[0]; qeal dFdc22_y = -2.0 * sC3.data()[1]; qeal dFd22_z = -2.0 * sC3.data()[2];

		//
		qeal q1 = 0.5 * (B + D) / (A * A);  qeal q2 = -0.5 / A;
		qeal dALdc11_x = q1 * dAdc11_x + q2 * (dBdc11_x + dDdc11_x);
		qeal dALdc11_y = q1 * dAdc11_y + q2 * (dBdc11_y + dDdc11_y);
		qeal dALdc11_z = q1 * dAdc11_z + q2 * (dBdc11_z + dDdc11_z);

		qeal dALdc12_x = q1 * dAdc12_x + q2 * (dBdc12_x + dDdc12_x);
		qeal dALdc12_y = q1 * dAdc12_y + q2 * (dBdc12_y + dDdc12_y);
		qeal dALdc12_z = q1 * dAdc12_z + q2 * (dBdc12_z + dDdc12_z);

		qeal dALdc21_x = q2 * dBdc21_x;
		qeal dALdc21_y = q2 * dBdc21_y;
		qeal dALdc21_z = q2 * dBdc21_z;

		qeal dALdc22_x = q2 * (dBdc22_x + dDdc22_x);
		qeal dALdc22_y = q2 * (dBdc22_x + dDdc22_y);
		qeal dALdc22_z = q2 * (dBdc22_x + dDdc22_z);

		qeal dBTdc11_x = 0;
		qeal dBTdc11_y = 0;
		qeal dBTdc11_z = 0;

		qeal dBTdc12_x = 0;
		qeal dBTdc12_y = 0;
		qeal dBTdc12_z = 0;

		qeal dBTdc21_x = 0;
		qeal dBTdc21_y = 0;
		qeal dBTdc21_z = 0;

		qeal dBTdc22_x = 0;
		qeal dBTdc22_y = 0;
		qeal dBTdc22_z = 0;

		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (alpha * C1x + beta * C2x + C3x);
		qeal vy = 2.0 * (alpha * C1y + beta * C2y + C3y);
		qeal vz = 2.0 * (alpha * C1z + beta * C2z + C3z);

		qeal dvxdc11_x = 2.0 * (C1x * dALdc11_x + C2x * dBTdc11_x + alpha);
		qeal dvxdc11_y = 2.0 * (C1x * dALdc11_y + C2x * dBTdc11_y);
		qeal dvxdc11_z = 2.0 * (C1x * dALdc11_z + C2x * dBTdc11_z);

		qeal dvxdc12_x = 2.0 * (C1x * dALdc12_x + C2x * dBTdc12_x - alpha + 1);
		qeal dvxdc12_y = 2.0 * (C1x * dALdc12_y + C2x * dBTdc12_y);
		qeal dvxdc12_z = 2.0 * (C1x * dALdc12_z + C2x * dBTdc12_z);

		qeal dvxdc21_x = 2.0 * (C1x * dALdc21_x + C2x * dBTdc21_x - beta);
		qeal dvxdc21_y = 2.0 * (C1x * dALdc21_y + C2x * dBTdc21_y);
		qeal dvxdc21_z = 2.0 * (C1x * dALdc21_z + C2x * dBTdc21_z);

		qeal dvxdc22_x = 2.0 * (C1x * dALdc22_x + C2x * dBTdc22_x + beta - 1);
		qeal dvxdc22_y = 2.0 * (C1x * dALdc22_y + C2x * dBTdc22_y);
		qeal dvxdc22_z = 2.0 * (C1x * dALdc22_z + C2x * dBTdc22_z);
		//
		qeal dvydc11_x = 2.0 * (C1y * dALdc11_x + C2y * dBTdc11_x);
		qeal dvydc11_y = 2.0 * (C1y * dALdc11_y + C2y * dBTdc11_y + alpha);
		qeal dvydc11_z = 2.0 * (C1y * dALdc11_z + C2y * dBTdc11_z);

		qeal dvydc12_x = 2.0 * (C1y * dALdc12_x + C2y * dBTdc12_x);
		qeal dvydc12_y = 2.0 * (C1y * dALdc12_y + C2y * dBTdc12_y - alpha + 1);
		qeal dvydc12_z = 2.0 * (C1y * dALdc12_z + C2y * dBTdc12_z);

		qeal dvydc21_x = 2.0 * (C1y * dALdc21_x + C2y * dBTdc21_x);
		qeal dvydc21_y = 2.0 * (C1y * dALdc21_y + C2y * dBTdc21_y - beta);
		qeal dvydc21_z = 2.0 * (C1y * dALdc21_z + C2y * dBTdc21_z);

		qeal dvydc22_x = 2.0 * (C1y * dALdc22_x + C2y * dBTdc22_x);
		qeal dvydc22_y = 2.0 * (C1y * dALdc22_y + C2y * dBTdc22_y + beta - 1);
		qeal dvydc22_z = 2.0 * (C1y * dALdc22_z + C2y * dBTdc22_z);
		//
		qeal dvzdc11_x = 2.0 * (C1z * dALdc11_x + C2z * dBTdc11_x);
		qeal dvzdc11_y = 2.0 * (C1z * dALdc11_y + C2z * dBTdc11_y);
		qeal dvzdc11_z = 2.0 * (C1z * dALdc11_z + C2z * dBTdc11_z + alpha);

		qeal dvzdc12_x = 2.0 * (C1z * dALdc12_x + C2z * dBTdc12_x);
		qeal dvzdc12_y = 2.0 * (C1z * dALdc12_y + C2z * dBTdc12_y);
		qeal dvzdc12_z = 2.0 * (C1z * dALdc12_z + C2z * dBTdc12_z - alpha + 1);

		qeal dvzdc21_x = 2.0 * (C1z * dALdc21_x + C2z * dBTdc21_x);
		qeal dvzdc21_y = 2.0 * (C1z * dALdc21_y + C2z * dBTdc21_y);
		qeal dvzdc21_z = 2.0 * (C1z * dALdc21_z + C2z * dBTdc21_z - beta);

		qeal dvzdc22_x = 2.0 * (C1z * dALdc22_x + C2z * dBTdc22_x);
		qeal dvzdc22_y = 2.0 * (C1z * dALdc22_y + C2z * dBTdc22_y);
		qeal dvzdc22_z = 2.0 * (C1z * dALdc22_z + C2z * dBTdc22_z + beta - 1);

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = vx * dALdc11_x + alpha * dvxdc11_x;
		hessina.data()[1] = vy * dALdc11_x + alpha * dvydc11_x;
		hessina.data()[2] = vz * dALdc11_x + alpha * dvzdc11_x;

		hessina.data()[3] = vx * -dALdc11_x + (1 - alpha) * dvxdc11_x;
		hessina.data()[4] = vy * -dALdc11_x + (1 - alpha) * dvydc11_x;
		hessina.data()[5] = vz * -dALdc11_x + (1 - alpha) * dvzdc11_x;

		hessina.data()[6] = vx * -dBTdc11_x - beta * dvxdc11_x;
		hessina.data()[7] = vy * -dBTdc11_x - beta * dvydc11_x;
		hessina.data()[8] = vz * -dBTdc11_x - beta * dvzdc11_x;

		hessina.data()[9] = vx * dBTdc11_x + (beta - 1.0) * dvxdc11_x;
		hessina.data()[10] = vy * dBTdc11_x + (beta - 1.0) * dvydc11_x;
		hessina.data()[11] = vz * dBTdc11_x + (beta - 1.0) * dvzdc11_x;
		// f / c11y
		hessina.data()[12] = vx * dALdc11_y + alpha * dvxdc11_y;
		hessina.data()[13] = vy * dALdc11_y + alpha * dvydc11_y;
		hessina.data()[14] = vz * dALdc11_y + alpha * dvzdc11_y;

		hessina.data()[15] = vx * -dALdc11_y + (1 - alpha) * dvxdc11_y;
		hessina.data()[16] = vy * -dALdc11_y + (1 - alpha) * dvydc11_y;
		hessina.data()[17] = vz * -dALdc11_y + (1 - alpha) * dvzdc11_y;

		hessina.data()[18] = vx * -dBTdc11_y - beta * dvxdc11_y;
		hessina.data()[19] = vy * -dBTdc11_y - beta * dvydc11_y;
		hessina.data()[20] = vz * -dBTdc11_y - beta * dvzdc11_y;

		hessina.data()[21] = vx * dBTdc11_y + (beta - 1.0) * dvxdc11_y;
		hessina.data()[22] = vy * dBTdc11_y + (beta - 1.0) * dvydc11_y;
		hessina.data()[23] = vz * dBTdc11_y + (beta - 1.0) * dvzdc11_y;
		// f / c11z
		hessina.data()[24] = vx * dALdc11_z + alpha * dvxdc11_z;
		hessina.data()[25] = vy * dALdc11_z + alpha * dvydc11_z;
		hessina.data()[26] = vz * dALdc11_z + alpha * dvzdc11_z;

		hessina.data()[27] = vx * -dALdc11_z + (1 - alpha) * dvxdc11_z;
		hessina.data()[28] = vy * -dALdc11_z + (1 - alpha) * dvydc11_z;
		hessina.data()[29] = vz * -dALdc11_z + (1 - alpha) * dvzdc11_z;

		hessina.data()[30] = vx * -dBTdc11_z - beta * dvxdc11_z;
		hessina.data()[31] = vy * -dBTdc11_z - beta * dvydc11_z;
		hessina.data()[32] = vz * -dBTdc11_z - beta * dvzdc11_z;

		hessina.data()[33] = vx * dBTdc11_z + (beta - 1.0) * dvxdc11_z;
		hessina.data()[34] = vy * dBTdc11_z + (beta - 1.0) * dvydc11_z;
		hessina.data()[35] = vz * dBTdc11_z + (beta - 1.0) * dvzdc11_z;

		// f / c12x
		hessina.data()[36] = vx * dALdc12_x + alpha * dvxdc12_x;
		hessina.data()[37] = vy * dALdc12_x + alpha * dvydc12_x;
		hessina.data()[38] = vz * dALdc12_x + alpha * dvzdc12_x;

		hessina.data()[39] = vx * -dALdc12_x + (1 - alpha) * dvxdc12_x;
		hessina.data()[40] = vy * -dALdc12_x + (1 - alpha) * dvydc12_x;
		hessina.data()[41] = vz * -dALdc12_x + (1 - alpha) * dvzdc12_x;

		hessina.data()[42] = vx * -dBTdc12_x - beta * dvxdc12_x;
		hessina.data()[43] = vy * -dBTdc12_x - beta * dvydc12_x;
		hessina.data()[44] = vz * -dBTdc12_x - beta * dvzdc12_x;

		hessina.data()[45] = vx * dBTdc12_x + (beta - 1.0) * dvxdc12_x;
		hessina.data()[46] = vy * dBTdc12_x + (beta - 1.0) * dvydc12_x;
		hessina.data()[47] = vz * dBTdc12_x + (beta - 1.0) * dvzdc12_x;
		// f / c12y
		hessina.data()[48] = vx * dALdc12_y + alpha * dvxdc12_y;
		hessina.data()[49] = vy * dALdc12_y + alpha * dvydc12_y;
		hessina.data()[50] = vz * dALdc12_y + alpha * dvzdc12_y;

		hessina.data()[51] = vx * -dALdc12_y + (1 - alpha) * dvxdc12_y;
		hessina.data()[52] = vy * -dALdc12_y + (1 - alpha) * dvydc12_y;
		hessina.data()[53] = vz * -dALdc12_y + (1 - alpha) * dvzdc12_y;

		hessina.data()[54] = vx * -dBTdc12_y - beta * dvxdc12_y;
		hessina.data()[55] = vy * -dBTdc12_y - beta * dvydc12_y;
		hessina.data()[56] = vz * -dBTdc12_y - beta * dvzdc12_y;

		hessina.data()[57] = vx * dBTdc12_y + (beta - 1.0) * dvxdc12_y;
		hessina.data()[58] = vy * dBTdc12_y + (beta - 1.0) * dvydc12_y;
		hessina.data()[59] = vz * dBTdc12_y + (beta - 1.0) * dvzdc12_y;
		// f / c12z
		hessina.data()[60] = vx * dALdc12_z + alpha * dvxdc12_z;
		hessina.data()[61] = vy * dALdc12_z + alpha * dvydc12_z;
		hessina.data()[62] = vz * dALdc12_z + alpha * dvzdc12_z;

		hessina.data()[63] = vx * -dALdc12_z + (1 - alpha) * dvxdc12_z;
		hessina.data()[64] = vy * -dALdc12_z + (1 - alpha) * dvydc12_z;
		hessina.data()[65] = vz * -dALdc12_z + (1 - alpha) * dvzdc12_z;

		hessina.data()[66] = vx * -dBTdc12_z - beta * dvxdc12_z;
		hessina.data()[67] = vy * -dBTdc12_z - beta * dvydc12_z;
		hessina.data()[68] = vz * -dBTdc12_z - beta * dvzdc12_z;

		hessina.data()[69] = vx * dBTdc12_z + (beta - 1.0) * dvxdc12_z;
		hessina.data()[70] = vy * dBTdc12_z + (beta - 1.0) * dvydc12_z;
		hessina.data()[71] = vz * dBTdc12_z + (beta - 1.0) * dvzdc12_z;

		// f / c21x
		hessina.data()[72] = vx * dALdc21_x + alpha * dvxdc21_x;
		hessina.data()[73] = vy * dALdc21_x + alpha * dvydc21_x;
		hessina.data()[74] = vz * dALdc21_x + alpha * dvzdc21_x;

		hessina.data()[75] = vx * -dALdc21_x + (1 - alpha) * dvxdc21_x;
		hessina.data()[76] = vy * -dALdc21_x + (1 - alpha) * dvydc21_x;
		hessina.data()[77] = vz * -dALdc21_x + (1 - alpha) * dvzdc21_x;

		hessina.data()[78] = vx * -dBTdc21_x - beta * dvxdc21_x;
		hessina.data()[79] = vy * -dBTdc21_x - beta * dvydc21_x;
		hessina.data()[80] = vz * -dBTdc21_x - beta * dvzdc21_x;

		hessina.data()[81] = vx * dBTdc21_x + (beta - 1.0) * dvxdc21_x;
		hessina.data()[82] = vy * dBTdc21_x + (beta - 1.0) * dvydc21_x;
		hessina.data()[83] = vz * dBTdc21_x + (beta - 1.0) * dvzdc21_x;
		// f / c21y
		hessina.data()[84] = vx * dALdc21_y + alpha * dvxdc21_y;
		hessina.data()[85] = vy * dALdc21_y + alpha * dvydc21_y;
		hessina.data()[86] = vz * dALdc21_y + alpha * dvzdc21_y;

		hessina.data()[87] = vx * -dALdc21_y + (1 - alpha) * dvxdc21_y;
		hessina.data()[88] = vy * -dALdc21_y + (1 - alpha) * dvydc21_y;
		hessina.data()[89] = vz * -dALdc21_y + (1 - alpha) * dvzdc21_y;

		hessina.data()[90] = vx * -dBTdc21_y - beta * dvxdc21_y;
		hessina.data()[91] = vy * -dBTdc21_y - beta * dvydc21_y;
		hessina.data()[92] = vz * -dBTdc21_y - beta * dvzdc21_y;

		hessina.data()[93] = vx * dBTdc21_y + (beta - 1.0) * dvxdc21_y;
		hessina.data()[94] = vy * dBTdc21_y + (beta - 1.0) * dvydc21_y;
		hessina.data()[95] = vz * dBTdc21_y + (beta - 1.0) * dvzdc21_y;
		// f / c21z
		hessina.data()[96] = vx * dALdc21_z + alpha * dvxdc21_z;
		hessina.data()[97] = vy * dALdc21_z + alpha * dvydc21_z;
		hessina.data()[98] = vz * dALdc21_z + alpha * dvzdc21_z;

		hessina.data()[99] = vx * -dALdc21_z + (1 - alpha) * dvxdc21_z;
		hessina.data()[100] = vy * -dALdc21_z + (1 - alpha) * dvydc21_z;
		hessina.data()[101] = vz * -dALdc21_z + (1 - alpha) * dvzdc21_z;

		hessina.data()[102] = vx * -dBTdc21_z - beta * dvxdc21_z;
		hessina.data()[103] = vy * -dBTdc21_z - beta * dvydc21_z;
		hessina.data()[104] = vz * -dBTdc21_z - beta * dvzdc21_z;

		hessina.data()[105] = vx * dBTdc21_z + (beta - 1.0) * dvxdc21_z;
		hessina.data()[106] = vy * dBTdc21_z + (beta - 1.0) * dvydc21_z;
		hessina.data()[107] = vz * dBTdc21_z + (beta - 1.0) * dvzdc21_z;

		// f / c22x
		hessina.data()[108] = vx * dALdc22_x + alpha * dvxdc22_x;
		hessina.data()[109] = vy * dALdc22_x + alpha * dvydc22_x;
		hessina.data()[110] = vz * dALdc22_x + alpha * dvzdc22_x;

		hessina.data()[111] = vx * -dALdc22_x + (1 - alpha) * dvxdc22_x;
		hessina.data()[112] = vy * -dALdc22_x + (1 - alpha) * dvydc22_x;
		hessina.data()[113] = vz * -dALdc22_x + (1 - alpha) * dvzdc22_x;

		hessina.data()[114] = vx * -dBTdc22_x - beta * dvxdc22_x;
		hessina.data()[115] = vy * -dBTdc22_x - beta * dvydc22_x;
		hessina.data()[116] = vz * -dBTdc22_x - beta * dvzdc22_x;

		hessina.data()[117] = vx * dBTdc22_x + (beta - 1.0) * dvxdc22_x;
		hessina.data()[118] = vy * dBTdc22_x + (beta - 1.0) * dvydc22_x;
		hessina.data()[119] = vz * dBTdc22_x + (beta - 1.0) * dvzdc22_x;
		// f / c22y
		hessina.data()[120] = vx * dALdc22_y + alpha * dvxdc22_y;
		hessina.data()[121] = vy * dALdc22_y + alpha * dvydc22_y;
		hessina.data()[122] = vz * dALdc22_y + alpha * dvzdc22_y;

		hessina.data()[123] = vx * -dALdc22_y + (1 - alpha) * dvxdc22_y;
		hessina.data()[124] = vy * -dALdc22_y + (1 - alpha) * dvydc22_y;
		hessina.data()[125] = vz * -dALdc22_y + (1 - alpha) * dvzdc22_y;

		hessina.data()[126] = vx * -dBTdc22_y - beta * dvxdc22_y;
		hessina.data()[127] = vy * -dBTdc22_y - beta * dvydc22_y;
		hessina.data()[128] = vz * -dBTdc22_y - beta * dvzdc22_y;

		hessina.data()[129] = vx * dBTdc22_y + (beta - 1.0) * dvxdc22_y;
		hessina.data()[130] = vy * dBTdc22_y + (beta - 1.0) * dvydc22_y;
		hessina.data()[131] = vz * dBTdc22_y + (beta - 1.0) * dvzdc22_y;
		// f / c22z
		hessina.data()[132] = vx * dALdc22_z + alpha * dvxdc22_z;
		hessina.data()[133] = vy * dALdc22_z + alpha * dvydc22_z;
		hessina.data()[134] = vz * dALdc22_z + alpha * dvzdc22_z;

		hessina.data()[135] = vx * -dALdc22_z + (1 - alpha) * dvxdc22_z;
		hessina.data()[136] = vy * -dALdc22_z + (1 - alpha) * dvydc22_z;
		hessina.data()[137] = vz * -dALdc22_z + (1 - alpha) * dvzdc22_z;

		hessina.data()[138] = vx * -dBTdc22_z - beta * dvxdc22_z;
		hessina.data()[139] = vy * -dBTdc22_z - beta * dvydc22_z;
		hessina.data()[140] = vz * -dBTdc22_z - beta * dvzdc22_z;

		hessina.data()[141] = vx * dBTdc22_z + (beta - 1.0) * dvxdc22_z;
		hessina.data()[142] = vy * dBTdc22_z + (beta - 1.0) * dvydc22_z;
		hessina.data()[143] = vz * dBTdc22_z + (beta - 1.0) * dvzdc22_z;
	}

	void MipcConeConeConstraint::alphaBetaHessina(MatrixX & hessina)
	{
		qeal dAdc11_x = 2.0 * sC1.data()[0]; qeal dAdc11_y = 2.0 * sC1.data()[1]; qeal dAdc11_z = 2.0 * sC1.data()[2];
		qeal dAdc12_x = -2.0 * sC1.data()[0]; qeal dAdc12_y = -2.0 * sC1.data()[1]; qeal dAdc12_z = -2.0 * sC1.data()[2];

		qeal dBdc11_x = 2.0 * sC2.data()[0]; qeal dBdc11_y = 2.0 * sC2.data()[1]; qeal dBdc11_z = 2.0 * sC2.data()[2];
		qeal dBdc12_x = -2.0 * sC2.data()[0]; qeal dBdc12_y = -2.0 * sC2.data()[1]; qeal dBdc12_z = -2.0 * sC2.data()[2];
		qeal dBdc21_x = -2.0 * sC1.data()[0]; qeal dBdc21_y = -2.0 * sC1.data()[1]; qeal dBdc21_z = -2.0 * sC1.data()[2];
		qeal dBdc22_x = 2.0 * sC1.data()[0]; qeal dBdc22_y = 2.0 * sC1.data()[1]; qeal dBdc22_z = 2.0 * sC1.data()[2];

		qeal dCdc21_x = -2.0 * sC2.data()[0]; qeal dCdc21_y = -2.0 * sC2.data()[1]; qeal dCdc21_z = -2.0 * sC2.data()[2];
		qeal dCdc22_x = 2.0 * sC2.data()[0]; qeal dCdc22_y = 2.0 * sC2.data()[1]; qeal dCdc22_z = 2.0 * sC2.data()[2];

		qeal dDdc11_x = 2.0 * sC3.data()[0]; qeal dDdc11_y = 2.0 * sC3.data()[1]; qeal dDdc11_z = 2.0 * sC3.data()[2];
		qeal dDdc12_x = 2.0 * (sC1.data()[0] - sC3.data()[0]); qeal dDdc12_y = 2.0 * (sC1.data()[1] - sC3.data()[1]); qeal dDdc12_z = 2.0 * (sC1.data()[2] - sC3.data()[2]);
		qeal dDdc22_x = -2.0 * sC1.data()[0]; qeal dDdc22_y = -2.0 * sC1.data()[1]; qeal dDdc22_z = -2.0 * sC1.data()[2];

		qeal dEdc12_x = 2.0 * sC2.data()[0]; qeal dEdc12_y = 2.0 * sC2.data()[1]; qeal dEdc12_z = 2.0 * sC2.data()[2];
		qeal dEdc21_x = -2.0 * sC3.data()[0]; qeal dEdc21_y = -2.0 * sC3.data()[1]; qeal dEdc21_z = -2.0 * sC3.data()[2];
		qeal dEdc22_x = 2.0 * (sC3.data()[0] - sC2.data()[0]); qeal dEdc22_y = 2.0 * (sC3.data()[1] - sC2.data()[1]); qeal dEdc22_z = 2.0 * (sC3.data()[2] - sC2.data()[2]);

		qeal dFdc12_x = 2.0 * sC3.data()[0]; qeal dFdc12_y = 2.0 * sC3.data()[1]; qeal dFdc12_z = 2.0 * sC3.data()[2];
		qeal dFdc22 = -2.0 * sC3.data()[0]; qeal dFdc22_y = -2.0 * sC3.data()[1]; qeal dFd22_z = -2.0 * sC3.data()[2];

		//
		qeal q1 = 1.0 / delta; qeal q2 = alpha / delta; qeal q3 = beta / delta;

		qeal dALdc11_x = q1 * (E * dBdc11_x - 2.0 * C * dDdc11_x) - q2 * (4.0 * C * dAdc11_x - 2.0 * B * dBdc11_x);
		qeal dALdc11_y = q1 * (E * dBdc11_y - 2.0 * C * dDdc11_y) - q2 * (4.0 * C * dAdc11_y - 2.0 * B * dBdc11_y);
		qeal dALdc11_z = q1 * (E * dBdc11_z - 2.0 * C * dDdc11_z) - q2 * (4.0 * C * dAdc11_z - 2.0 * B * dBdc11_z);

		qeal dALdc12_x = q1 * (E * dBdc12_x + B * dEdc12_x - 2.0 * C * dDdc12_x) - q2 * (4.0 * C * dAdc12_x - 2.0 * B * dBdc12_x);
		qeal dALdc12_y = q1 * (E * dBdc12_y + B * dEdc12_y - 2.0 * C * dDdc12_y) - q2 * (4.0 * C * dAdc12_y - 2.0 * B * dBdc12_y);
		qeal dALdc12_z = q1 * (E * dBdc12_z + B * dEdc12_z - 2.0 * C * dDdc12_z) - q2 * (4.0 * C * dAdc12_z - 2.0 * B * dBdc12_z);

		qeal dALdc21_x = q1 * (E * dBdc21_x + B * dEdc21_x - 2.0 * D * dCdc21_x) - q2 * (4.0 * A * dCdc21_x - 2.0 * B * dBdc21_x);
		qeal dALdc21_y = q1 * (E * dBdc21_y + B * dEdc21_y - 2.0 * D * dCdc21_y) - q2 * (4.0 * A * dCdc21_y - 2.0 * B * dBdc21_y);
		qeal dALdc21_z = q1 * (E * dBdc21_z + B * dEdc21_z - 2.0 * D * dCdc21_z) - q2 * (4.0 * A * dCdc21_z - 2.0 * B * dBdc21_z);

		qeal dALdc22_x = q1 * (E * dBdc22_x + B * dEdc22_x - 2.0 * (D * dCdc22_x + C * dDdc22_x)) - q2 * (4.0 * A * dCdc22_x - 2.0 * B * dBdc22_x);
		qeal dALdc22_y = q1 * (E * dBdc22_y + B * dEdc22_y - 2.0 * (D * dCdc22_y + C * dDdc22_y)) - q2 * (4.0 * A * dCdc22_y - 2.0 * B * dBdc22_y);
		qeal dALdc22_z = q1 * (E * dBdc22_z + B * dEdc22_z - 2.0 * (D * dCdc22_z + C * dDdc22_z)) - q2 * (4.0 * A * dCdc22_z - 2.0 * B * dBdc22_z);

		qeal dBTdc11_x = q1 * (D * dBdc11_x + B * dDdc11_x - 2.0 * E * dAdc11_x) - q3 * (4.0 * C * dAdc11_x - 2.0 * B * dBdc11_x);
		qeal dBTdc11_y = q1 * (D * dBdc11_y + B * dDdc11_y - 2.0 * E * dAdc11_y) - q3 * (4.0 * C * dAdc11_y - 2.0 * B * dBdc11_y);
		qeal dBTdc11_z = q1 * (D * dBdc11_z + B * dDdc11_z - 2.0 * E * dAdc11_z) - q3 * (4.0 * C * dAdc11_z - 2.0 * B * dBdc11_z);

		qeal dBTdc12_x = q1 * (D * dBdc12_x + B * dDdc12_x - 2.0 * (E * dAdc12_x + A * dEdc12_x)) - q3 * (4.0 * C * dAdc12_x - 2.0 * B * dBdc12_x);
		qeal dBTdc12_y = q1 * (D * dBdc12_y + B * dDdc12_y - 2.0 * (E * dAdc12_y + A * dEdc12_y)) - q3 * (4.0 * C * dAdc12_y - 2.0 * B * dBdc12_y);
		qeal dBTdc12_z = q1 * (D * dBdc12_z + B * dDdc12_z - 2.0 * (E * dAdc12_z + A * dEdc12_z)) - q3 * (4.0 * C * dAdc12_z - 2.0 * B * dBdc12_z);

		qeal dBTdc21_x = q1 * (D * dBdc21_x - 2.0 * A * dEdc21_x) - q3 * (4.0 * A * dCdc21_x - 2.0 * B * dBdc21_x);
		qeal dBTdc21_y = q1 * (D * dBdc21_y - 2.0 * A * dEdc21_y) - q3 * (4.0 * A * dCdc21_y - 2.0 * B * dBdc21_y);
		qeal dBTdc21_z = q1 * (D * dBdc21_z - 2.0 * A * dEdc21_z) - q3 * (4.0 * A * dCdc21_z - 2.0 * B * dBdc21_z);

		qeal dBTdc22_x = q1 * (D * dBdc22_x + B * dDdc22_x - 2.0 * A * dEdc22_x) - q3 * (4.0 * A * dCdc22_x - 2.0 * B * dBdc22_x);
		qeal dBTdc22_y = q1 * (D * dBdc22_y + B * dDdc22_y - 2.0 * A * dEdc22_y) - q3 * (4.0 * A * dCdc22_y - 2.0 * B * dBdc22_y);
		qeal dBTdc22_z = q1 * (D * dBdc22_z + B * dDdc22_z - 2.0 * A * dEdc22_z) - q3 * (4.0 * A * dCdc22_z - 2.0 * B * dBdc22_z);

		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (alpha * C1x + beta * C2x + C3x);
		qeal vy = 2.0 * (alpha * C1y + beta * C2y + C3y);
		qeal vz = 2.0 * (alpha * C1z + beta * C2z + C3z);

		qeal dvxdc11_x = 2.0 * (C1x * dALdc11_x + C2x * dBTdc11_x + alpha);
		qeal dvxdc11_y = 2.0 * (C1x * dALdc11_y + C2x * dBTdc11_y);
		qeal dvxdc11_z = 2.0 * (C1x * dALdc11_z + C2x * dBTdc11_z);

		qeal dvxdc12_x = 2.0 * (C1x * dALdc12_x + C2x * dBTdc12_x - alpha + 1);
		qeal dvxdc12_y = 2.0 * (C1x * dALdc12_y + C2x * dBTdc12_y);
		qeal dvxdc12_z = 2.0 * (C1x * dALdc12_z + C2x * dBTdc12_z);

		qeal dvxdc21_x = 2.0 * (C1x * dALdc21_x + C2x * dBTdc21_x - beta);
		qeal dvxdc21_y = 2.0 * (C1x * dALdc21_y + C2x * dBTdc21_y);
		qeal dvxdc21_z = 2.0 * (C1x * dALdc21_z + C2x * dBTdc21_z);

		qeal dvxdc22_x = 2.0 * (C1x * dALdc22_x + C2x * dBTdc22_x + beta - 1);
		qeal dvxdc22_y = 2.0 * (C1x * dALdc22_y + C2x * dBTdc22_y);
		qeal dvxdc22_z = 2.0 * (C1x * dALdc22_z + C2x * dBTdc22_z);
		//
		qeal dvydc11_x = 2.0 * (C1y * dALdc11_x + C2y * dBTdc11_x);
		qeal dvydc11_y = 2.0 * (C1y * dALdc11_y + C2y * dBTdc11_y + alpha);
		qeal dvydc11_z = 2.0 * (C1y * dALdc11_z + C2y * dBTdc11_z);

		qeal dvydc12_x = 2.0 * (C1y * dALdc12_x + C2y * dBTdc12_x);
		qeal dvydc12_y = 2.0 * (C1y * dALdc12_y + C2y * dBTdc12_y - alpha + 1);
		qeal dvydc12_z = 2.0 * (C1y * dALdc12_z + C2y * dBTdc12_z);

		qeal dvydc21_x = 2.0 * (C1y * dALdc21_x + C2y * dBTdc21_x);
		qeal dvydc21_y = 2.0 * (C1y * dALdc21_y + C2y * dBTdc21_y - beta);
		qeal dvydc21_z = 2.0 * (C1y * dALdc21_z + C2y * dBTdc21_z);

		qeal dvydc22_x = 2.0 * (C1y * dALdc22_x + C2y * dBTdc22_x);
		qeal dvydc22_y = 2.0 * (C1y * dALdc22_y + C2y * dBTdc22_y + beta - 1);
		qeal dvydc22_z = 2.0 * (C1y * dALdc22_z + C2y * dBTdc22_z);
		//
		qeal dvzdc11_x = 2.0 * (C1z * dALdc11_x + C2z * dBTdc11_x);
		qeal dvzdc11_y = 2.0 * (C1z * dALdc11_y + C2z * dBTdc11_y);
		qeal dvzdc11_z = 2.0 * (C1z * dALdc11_z + C2z * dBTdc11_z + alpha);

		qeal dvzdc12_x = 2.0 * (C1z * dALdc12_x + C2z * dBTdc12_x);
		qeal dvzdc12_y = 2.0 * (C1z * dALdc12_y + C2z * dBTdc12_y);
		qeal dvzdc12_z = 2.0 * (C1z * dALdc12_z + C2z * dBTdc12_z - alpha + 1);

		qeal dvzdc21_x = 2.0 * (C1z * dALdc21_x + C2z * dBTdc21_x);
		qeal dvzdc21_y = 2.0 * (C1z * dALdc21_y + C2z * dBTdc21_y);
		qeal dvzdc21_z = 2.0 * (C1z * dALdc21_z + C2z * dBTdc21_z - beta);

		qeal dvzdc22_x = 2.0 * (C1z * dALdc22_x + C2z * dBTdc22_x);
		qeal dvzdc22_y = 2.0 * (C1z * dALdc22_y + C2z * dBTdc22_y);
		qeal dvzdc22_z = 2.0 * (C1z * dALdc22_z + C2z * dBTdc22_z + beta - 1);

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = vx * dALdc11_x + alpha * dvxdc11_x;
		hessina.data()[1] = vy * dALdc11_x + alpha * dvydc11_x;
		hessina.data()[2] = vz * dALdc11_x + alpha * dvzdc11_x;

		hessina.data()[3] = vx * -dALdc11_x + (1 - alpha) * dvxdc11_x;
		hessina.data()[4] = vy * -dALdc11_x + (1 - alpha) * dvydc11_x;
		hessina.data()[5] = vz * -dALdc11_x + (1 - alpha) * dvzdc11_x;

		hessina.data()[6] = vx * -dBTdc11_x - beta * dvxdc11_x;
		hessina.data()[7] = vy * -dBTdc11_x - beta * dvydc11_x;
		hessina.data()[8] = vz * -dBTdc11_x - beta * dvzdc11_x;

		hessina.data()[9] = vx * dBTdc11_x + (beta - 1.0) * dvxdc11_x;
		hessina.data()[10] = vy * dBTdc11_x + (beta - 1.0) * dvydc11_x;
		hessina.data()[11] = vz * dBTdc11_x + (beta - 1.0) * dvzdc11_x;
		// f / c11y
		hessina.data()[12] = vx * dALdc11_y + alpha * dvxdc11_y;
		hessina.data()[13] = vy * dALdc11_y + alpha * dvydc11_y;
		hessina.data()[14] = vz * dALdc11_y + alpha * dvzdc11_y;

		hessina.data()[15] = vx * -dALdc11_y + (1 - alpha) * dvxdc11_y;
		hessina.data()[16] = vy * -dALdc11_y + (1 - alpha) * dvydc11_y;
		hessina.data()[17] = vz * -dALdc11_y + (1 - alpha) * dvzdc11_y;

		hessina.data()[18] = vx * -dBTdc11_y - beta * dvxdc11_y;
		hessina.data()[19] = vy * -dBTdc11_y - beta * dvydc11_y;
		hessina.data()[20] = vz * -dBTdc11_y - beta * dvzdc11_y;

		hessina.data()[21] = vx * dBTdc11_y + (beta - 1.0) * dvxdc11_y;
		hessina.data()[22] = vy * dBTdc11_y + (beta - 1.0) * dvydc11_y;
		hessina.data()[23] = vz * dBTdc11_y + (beta - 1.0) * dvzdc11_y;
		// f / c11z
		hessina.data()[24] = vx * dALdc11_z + alpha * dvxdc11_z;
		hessina.data()[25] = vy * dALdc11_z + alpha * dvydc11_z;
		hessina.data()[26] = vz * dALdc11_z + alpha * dvzdc11_z;

		hessina.data()[27] = vx * -dALdc11_z + (1 - alpha) * dvxdc11_z;
		hessina.data()[28] = vy * -dALdc11_z + (1 - alpha) * dvydc11_z;
		hessina.data()[29] = vz * -dALdc11_z + (1 - alpha) * dvzdc11_z;

		hessina.data()[30] = vx * -dBTdc11_z - beta * dvxdc11_z;
		hessina.data()[31] = vy * -dBTdc11_z - beta * dvydc11_z;
		hessina.data()[32] = vz * -dBTdc11_z - beta * dvzdc11_z;

		hessina.data()[33] = vx * dBTdc11_z + (beta - 1.0) * dvxdc11_z;
		hessina.data()[34] = vy * dBTdc11_z + (beta - 1.0) * dvydc11_z;
		hessina.data()[35] = vz * dBTdc11_z + (beta - 1.0) * dvzdc11_z;

		// f / c12x
		hessina.data()[36] = vx * dALdc12_x + alpha * dvxdc12_x;
		hessina.data()[37] = vy * dALdc12_x + alpha * dvydc12_x;
		hessina.data()[38] = vz * dALdc12_x + alpha * dvzdc12_x;

		hessina.data()[39] = vx * -dALdc12_x + (1 - alpha) * dvxdc12_x;
		hessina.data()[40] = vy * -dALdc12_x + (1 - alpha) * dvydc12_x;
		hessina.data()[41] = vz * -dALdc12_x + (1 - alpha) * dvzdc12_x;

		hessina.data()[42] = vx * -dBTdc12_x - beta * dvxdc12_x;
		hessina.data()[43] = vy * -dBTdc12_x - beta * dvydc12_x;
		hessina.data()[44] = vz * -dBTdc12_x - beta * dvzdc12_x;

		hessina.data()[45] = vx * dBTdc12_x + (beta - 1.0) * dvxdc12_x;
		hessina.data()[46] = vy * dBTdc12_x + (beta - 1.0) * dvydc12_x;
		hessina.data()[47] = vz * dBTdc12_x + (beta - 1.0) * dvzdc12_x;
		// f / c12y
		hessina.data()[48] = vx * dALdc12_y + alpha * dvxdc12_y;
		hessina.data()[49] = vy * dALdc12_y + alpha * dvydc12_y;
		hessina.data()[50] = vz * dALdc12_y + alpha * dvzdc12_y;

		hessina.data()[51] = vx * -dALdc12_y + (1 - alpha) * dvxdc12_y;
		hessina.data()[52] = vy * -dALdc12_y + (1 - alpha) * dvydc12_y;
		hessina.data()[53] = vz * -dALdc12_y + (1 - alpha) * dvzdc12_y;

		hessina.data()[54] = vx * -dBTdc12_y - beta * dvxdc12_y;
		hessina.data()[55] = vy * -dBTdc12_y - beta * dvydc12_y;
		hessina.data()[56] = vz * -dBTdc12_y - beta * dvzdc12_y;

		hessina.data()[57] = vx * dBTdc12_y + (beta - 1.0) * dvxdc12_y;
		hessina.data()[58] = vy * dBTdc12_y + (beta - 1.0) * dvydc12_y;
		hessina.data()[59] = vz * dBTdc12_y + (beta - 1.0) * dvzdc12_y;
		// f / c12z
		hessina.data()[60] = vx * dALdc12_z + alpha * dvxdc12_z;
		hessina.data()[61] = vy * dALdc12_z + alpha * dvydc12_z;
		hessina.data()[62] = vz * dALdc12_z + alpha * dvzdc12_z;

		hessina.data()[63] = vx * -dALdc12_z + (1 - alpha) * dvxdc12_z;
		hessina.data()[64] = vy * -dALdc12_z + (1 - alpha) * dvydc12_z;
		hessina.data()[65] = vz * -dALdc12_z + (1 - alpha) * dvzdc12_z;

		hessina.data()[66] = vx * -dBTdc12_z - beta * dvxdc12_z;
		hessina.data()[67] = vy * -dBTdc12_z - beta * dvydc12_z;
		hessina.data()[68] = vz * -dBTdc12_z - beta * dvzdc12_z;

		hessina.data()[69] = vx * dBTdc12_z + (beta - 1.0) * dvxdc12_z;
		hessina.data()[70] = vy * dBTdc12_z + (beta - 1.0) * dvydc12_z;
		hessina.data()[71] = vz * dBTdc12_z + (beta - 1.0) * dvzdc12_z;

		// f / c21x
		hessina.data()[72] = vx * dALdc21_x + alpha * dvxdc21_x;
		hessina.data()[73] = vy * dALdc21_x + alpha * dvydc21_x;
		hessina.data()[74] = vz * dALdc21_x + alpha * dvzdc21_x;

		hessina.data()[75] = vx * -dALdc21_x + (1 - alpha) * dvxdc21_x;
		hessina.data()[76] = vy * -dALdc21_x + (1 - alpha) * dvydc21_x;
		hessina.data()[77] = vz * -dALdc21_x + (1 - alpha) * dvzdc21_x;

		hessina.data()[78] = vx * -dBTdc21_x - beta * dvxdc21_x;
		hessina.data()[79] = vy * -dBTdc21_x - beta * dvydc21_x;
		hessina.data()[80] = vz * -dBTdc21_x - beta * dvzdc21_x;

		hessina.data()[81] = vx * dBTdc21_x + (beta - 1.0) * dvxdc21_x;
		hessina.data()[82] = vy * dBTdc21_x + (beta - 1.0) * dvydc21_x;
		hessina.data()[83] = vz * dBTdc21_x + (beta - 1.0) * dvzdc21_x;
		// f / c21y
		hessina.data()[84] = vx * dALdc21_y + alpha * dvxdc21_y;
		hessina.data()[85] = vy * dALdc21_y + alpha * dvydc21_y;
		hessina.data()[86] = vz * dALdc21_y + alpha * dvzdc21_y;

		hessina.data()[87] = vx * -dALdc21_y + (1 - alpha) * dvxdc21_y;
		hessina.data()[88] = vy * -dALdc21_y + (1 - alpha) * dvydc21_y;
		hessina.data()[89] = vz * -dALdc21_y + (1 - alpha) * dvzdc21_y;

		hessina.data()[90] = vx * -dBTdc21_y - beta * dvxdc21_y;
		hessina.data()[91] = vy * -dBTdc21_y - beta * dvydc21_y;
		hessina.data()[92] = vz * -dBTdc21_y - beta * dvzdc21_y;

		hessina.data()[93] = vx * dBTdc21_y + (beta - 1.0) * dvxdc21_y;
		hessina.data()[94] = vy * dBTdc21_y + (beta - 1.0) * dvydc21_y;
		hessina.data()[95] = vz * dBTdc21_y + (beta - 1.0) * dvzdc21_y;
		// f / c21z
		hessina.data()[96] = vx * dALdc21_z + alpha * dvxdc21_z;
		hessina.data()[97] = vy * dALdc21_z + alpha * dvydc21_z;
		hessina.data()[98] = vz * dALdc21_z + alpha * dvzdc21_z;

		hessina.data()[99] = vx * -dALdc21_z + (1 - alpha) * dvxdc21_z;
		hessina.data()[100] = vy * -dALdc21_z + (1 - alpha) * dvydc21_z;
		hessina.data()[101] = vz * -dALdc21_z + (1 - alpha) * dvzdc21_z;

		hessina.data()[102] = vx * -dBTdc21_z - beta * dvxdc21_z;
		hessina.data()[103] = vy * -dBTdc21_z - beta * dvydc21_z;
		hessina.data()[104] = vz * -dBTdc21_z - beta * dvzdc21_z;

		hessina.data()[105] = vx * dBTdc21_z + (beta - 1.0) * dvxdc21_z;
		hessina.data()[106] = vy * dBTdc21_z + (beta - 1.0) * dvydc21_z;
		hessina.data()[107] = vz * dBTdc21_z + (beta - 1.0) * dvzdc21_z;

		// f / c22x
		hessina.data()[108] = vx * dALdc22_x + alpha * dvxdc22_x;
		hessina.data()[109] = vy * dALdc22_x + alpha * dvydc22_x;
		hessina.data()[110] = vz * dALdc22_x + alpha * dvzdc22_x;

		hessina.data()[111] = vx * -dALdc22_x + (1 - alpha) * dvxdc22_x;
		hessina.data()[112] = vy * -dALdc22_x + (1 - alpha) * dvydc22_x;
		hessina.data()[113] = vz * -dALdc22_x + (1 - alpha) * dvzdc22_x;

		hessina.data()[114] = vx * -dBTdc22_x - beta * dvxdc22_x;
		hessina.data()[115] = vy * -dBTdc22_x - beta * dvydc22_x;
		hessina.data()[116] = vz * -dBTdc22_x - beta * dvzdc22_x;

		hessina.data()[117] = vx * dBTdc22_x + (beta - 1.0) * dvxdc22_x;
		hessina.data()[118] = vy * dBTdc22_x + (beta - 1.0) * dvydc22_x;
		hessina.data()[119] = vz * dBTdc22_x + (beta - 1.0) * dvzdc22_x;
		// f / c22y
		hessina.data()[120] = vx * dALdc22_y + alpha * dvxdc22_y;
		hessina.data()[121] = vy * dALdc22_y + alpha * dvydc22_y;
		hessina.data()[122] = vz * dALdc22_y + alpha * dvzdc22_y;

		hessina.data()[123] = vx * -dALdc22_y + (1 - alpha) * dvxdc22_y;
		hessina.data()[124] = vy * -dALdc22_y + (1 - alpha) * dvydc22_y;
		hessina.data()[125] = vz * -dALdc22_y + (1 - alpha) * dvzdc22_y;

		hessina.data()[126] = vx * -dBTdc22_y - beta * dvxdc22_y;
		hessina.data()[127] = vy * -dBTdc22_y - beta * dvydc22_y;
		hessina.data()[128] = vz * -dBTdc22_y - beta * dvzdc22_y;

		hessina.data()[129] = vx * dBTdc22_y + (beta - 1.0) * dvxdc22_y;
		hessina.data()[130] = vy * dBTdc22_y + (beta - 1.0) * dvydc22_y;
		hessina.data()[131] = vz * dBTdc22_y + (beta - 1.0) * dvzdc22_y;
		// f / c22z
		hessina.data()[132] = vx * dALdc22_z + alpha * dvxdc22_z;
		hessina.data()[133] = vy * dALdc22_z + alpha * dvydc22_z;
		hessina.data()[134] = vz * dALdc22_z + alpha * dvzdc22_z;

		hessina.data()[135] = vx * -dALdc22_z + (1 - alpha) * dvxdc22_z;
		hessina.data()[136] = vy * -dALdc22_z + (1 - alpha) * dvydc22_z;
		hessina.data()[137] = vz * -dALdc22_z + (1 - alpha) * dvzdc22_z;

		hessina.data()[138] = vx * -dBTdc22_z - beta * dvxdc22_z;
		hessina.data()[139] = vy * -dBTdc22_z - beta * dvydc22_z;
		hessina.data()[140] = vz * -dBTdc22_z - beta * dvzdc22_z;

		hessina.data()[141] = vx * dBTdc22_z + (beta - 1.0) * dvxdc22_z;
		hessina.data()[142] = vy * dBTdc22_z + (beta - 1.0) * dvydc22_z;
		hessina.data()[143] = vz * dBTdc22_z + (beta - 1.0) * dvzdc22_z;
	}

	qeal MipcSlabSphereConstraint::frictionEnergy()
	{
		Vector3 c11p, c12p, c13p, csp;
		spheres[0]->center->projectFullspacePreP(c11p.data());
		spheres[1]->center->projectFullspacePreP(c12p.data());
		spheres[2]->center->projectFullspacePreP(c13p.data());
		spheres[3]->center->projectFullspacePreP(csp.data());

		Vector3 rel_u;
		for (int i = 0; i < 3; i++)
			rel_u.data()[i] = lagAlpha * (spheres[0]->center->getP()[i] - c11p.data()[i]) + lagBbeta * (spheres[1]->center->getP()[i] - c12p.data()[i]) + (1.0 - lagAlpha - lagBbeta) * (spheres[2]->center->getP()[i] - c13p.data()[i]) - (spheres[3]->center->getP()[i] - csp.data()[i]);

		Vector2 rel_uk = lagBasis.transpose() * rel_u;
		qeal energy;
		f0_SF(rel_uk.squaredNorm(), epsvh, energy);

		energy *= mu * lagLamda;
		return energy;
	}

	qeal MipcSlabSphereConstraint::computeDistance()
	{
		for (int i = 0; i < 3; i++)
		{
			sC1.data()[i] = spheres[0]->center->getP()[i] - spheres[2]->center->getP()[i];
			sC2.data()[i] = spheres[1]->center->getP()[i] - spheres[2]->center->getP()[i];
			sC3.data()[i] = spheres[2]->center->getP()[i] - spheres[3]->center->getP()[i];
		}

		sR1 = *(spheres[0]->radius) - *(spheres[2]->radius);
		sR2 = *(spheres[1]->radius) - *(spheres[2]->radius);
		sR3 = *(spheres[2]->radius) + *(spheres[3]->radius);

		A = sC1.dot(sC1) - sR1 * sR1;
		B = 2.0 * (sC1.dot(sC2) - sR1 * sR2);
		C = sC2.dot(sC2) - sR2 * sR2;
		D = 2.0 * (sC1.dot(sC3) - sR1 * sR3);
		E = 2.0 * (sC2.dot(sC3) - sR2 * sR3);
		F = sC3.dot(sC3) - sR3 * sR3;

		delta = 4 * A * C - B * B;

		qeal temp_alpha, temp_beta;
		qeal temp_dist;
		distanceMode = TWO_ENDPOINTS;
		alpha = 0.0, beta = 0.0;
		distance = valueOfQuadircSurface2D(alpha, beta, A, B, C, D, E, F);

		temp_alpha = 1.0; temp_beta = 0.0;
		temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
		if (distance > temp_dist)
		{
			distance = temp_dist;
			alpha = temp_alpha; beta = temp_beta;
		}

		temp_alpha = 0.0; temp_beta = 1.0;
		temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
		if (distance > temp_dist)
		{
			distance = temp_dist;
			alpha = temp_alpha; beta = temp_beta;
		}

		temp_alpha = 0.0; temp_beta = -E / (2.0 *C);
		if (temp_beta > 0.0 && temp_beta < 1.0)
		{
			temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
			if (distance > temp_dist)
			{
				distance = temp_dist;
				alpha = temp_alpha; beta = temp_beta;
				distanceMode = ALPHA_ZERO;
			}
		}

		temp_alpha = -D / (2.0 *A); temp_beta = 0.0;
		if (temp_alpha > 0.0 && temp_alpha < 1.0)
		{
			temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
			if (distance > temp_dist)
			{
				distance = temp_dist;
				alpha = temp_alpha; beta = temp_beta;
				distanceMode = BETA_ZERO;
			}
		}

		temp_alpha = 0.5 * (2.0 * C + E - B - D) / (A - B + C); temp_beta = 1.0 - temp_alpha;
		if (temp_alpha > 0.0 && temp_alpha < 1.0)
		{
			temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
			if (distance > temp_dist)
			{
				distance = temp_dist;
				alpha = temp_alpha; beta = temp_beta;
				distanceMode = ALPHA_BETA_ONE;
			}
		}

		// can be ignored
		if (delta != 0.0)
		{
			temp_alpha = (B * E - 2.0 * C * D) / delta; temp_beta = (B * D - 2.0 * A * E) / delta;
			if (temp_alpha > 0.0 && temp_alpha < 1.0 && temp_beta> 0.0 && temp_beta < 1.0 && temp_alpha + temp_beta < 1.0)
			{
				temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
				if (distance > temp_dist)
				{
					distance = temp_dist;
					alpha = temp_alpha; beta = temp_beta;
					distanceMode = ALPHA_BETA_ONE;
				}
			}
		}

		Vector3 cp, cq;
		for (int i = 0; i < 3; i++)
		{
			cp.data()[i] = alpha * spheres[0]->center->getP()[i] + beta * spheres[1]->center->getP()[i] + (1.0 - alpha - beta) * spheres[2]->center->getP()[i];
			cq.data()[i] = spheres[3]->center->getP()[i];
		}
		cloestPoints[0] = cp;
		cloestPoints[1] = cq;

		return distance;
	}

	void MipcSlabSphereConstraint::getTanBasis(Eigen::Matrix<qeal, 3, 2>& lagBasis)
	{
		lagBasis.setZero();

		switch (distanceMode)
		{
		case TWO_ENDPOINTS: //PP
		{
			//	std::cout << "PP " << std::endl;
			Vector3 p;
			Vector3 q = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);

			if (alpha == 1.0)
			{
				p = Vector3(spheres[0]->center->getP()[0], spheres[0]->center->getP()[1], spheres[0]->center->getP()[2]);
			}
			else if (beta == 1.0)
			{
				p = Vector3(spheres[1]->center->getP()[0], spheres[1]->center->getP()[1], spheres[1]->center->getP()[2]);
			}
			else
			{
				p = Vector3(spheres[2]->center->getP()[0], spheres[2]->center->getP()[1], spheres[2]->center->getP()[2]);
			}
			Point_Point_Tangent_Basis(p, q, lagBasis);
			break;
		}
		case ALPHA_ZERO: // PE
		{
			//	std::cout << "PE " << std::endl;
			Vector3 e1 = Vector3(spheres[1]->center->getP()[0], spheres[1]->center->getP()[1], spheres[1]->center->getP()[2]);
			Vector3 e2 = Vector3(spheres[2]->center->getP()[0], spheres[2]->center->getP()[1], spheres[2]->center->getP()[2]);
			Vector3 q = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);
			Point_Edge_Tangent_Basis(q, e1, e2, lagBasis);
			break;
		}
		case BETA_ZERO: //PE
		{
			//	std::cout << "PE " << std::endl;
			Vector3 e0 = Vector3(spheres[0]->center->getP()[0], spheres[0]->center->getP()[1], spheres[0]->center->getP()[2]);
			Vector3 e2 = Vector3(spheres[2]->center->getP()[0], spheres[2]->center->getP()[1], spheres[2]->center->getP()[2]);
			Vector3 q = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);
			Point_Edge_Tangent_Basis(q, e0, e2, lagBasis);
			break;
		}
		case ALPHA_BETA_ONE: //PE
		{
			//	std::cout << "PE " << std::endl;
			Vector3 e0 = Vector3(spheres[0]->center->getP()[0], spheres[0]->center->getP()[1], spheres[0]->center->getP()[2]);
			Vector3 e1 = Vector3(spheres[1]->center->getP()[0], spheres[1]->center->getP()[1], spheres[1]->center->getP()[2]);
			Vector3 q = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);
			Point_Edge_Tangent_Basis(q, e0, e1, lagBasis);
			break;
		}
		case ALPHA_BETA: //PT
		{
			//	std::cout << "PT " << std::endl;
			Vector3 e0 = Vector3(spheres[0]->center->getP()[0], spheres[0]->center->getP()[1], spheres[0]->center->getP()[2]);
			Vector3 e1 = Vector3(spheres[1]->center->getP()[0], spheres[1]->center->getP()[1], spheres[1]->center->getP()[2]);
			Vector3 e2 = Vector3(spheres[2]->center->getP()[0], spheres[2]->center->getP()[1], spheres[2]->center->getP()[2]);
			Vector3 q = Vector3(spheres[3]->center->getP()[0], spheres[3]->center->getP()[1], spheres[3]->center->getP()[2]);
			Point_Triangle_Tangent_Basis(q, e2, e0, e1, lagBasis);
			break;
		}
		default:
			break;
		};
	}

	void MipcSlabSphereConstraint::computeLagTangentBasis(const qeal kappa)
	{
		getTanBasis(lagBasis);
		lagAlpha = alpha; lagBbeta = beta;
		lagLamda = -kappa * getBarrierGradient() * 2.0 * sqrt(distance);
	}

	void MipcSlabSphereConstraint::getGradientAndHessian(qeal kappa, VectorX& gradient, MatrixX& hessian)
	{
		qeal barrierGrad = getBarrierGradient();
		qeal barrierHessian = getBarrierHessian();
		VectorX distGrad;
		distGrad.resize(12);
		distGrad.setZero();
		MatrixX distHessina;
		distHessina.resize(12, 12);
		distHessina.setZero();
		diff_F_x(distGrad);

		fillOverallGradient(-1.0 * kappa * barrierGrad, distGrad, gradient);

		switch (distanceMode)
		{
		case TWO_ENDPOINTS:
			endPointsHessina(distHessina);
			break;
		case ALPHA_ZERO:
			alphaIsZeroHessina(distHessina);
			break;
		case BETA_ZERO:
			betaIsZeroHessina(distHessina);
			break;
		case ALPHA_BETA_ONE:
			alphaBetaPlusOneHessina(distHessina);
			break;
		case ALPHA_BETA:
			alphaBetaHessina(distHessina);
			break;
		default:
			break;
		};

		MatrixX hess = barrierGrad * distHessina + barrierHessian * (distGrad * distGrad.transpose());
		makePD(hess);

		hess *= kappa;
		fillOverallHessian(1.0, hess, hessian);
	}

	void MipcSlabSphereConstraint::getFrictionGradientAndHessian(qeal kappa, VectorX& gradient, MatrixX& hessian)
	{
		Vector3 c11p, c12p, c13p, csp;
		spheres[0]->center->projectFullspacePreP(c11p.data());
		spheres[1]->center->projectFullspacePreP(c12p.data());
		spheres[2]->center->projectFullspacePreP(c13p.data());
		spheres[3]->center->projectFullspacePreP(csp.data());

		Vector3 rel_u;
		for (int i = 0; i < 3; i++)
			rel_u.data()[i] = (lagAlpha * (spheres[0]->center->getP()[i] - c11p.data()[i]) + lagBbeta * (spheres[1]->center->getP()[i] - c12p.data()[i]) + (1.0 - lagAlpha - lagBbeta) * (spheres[2]->center->getP()[i] - c13p.data()[i])) - (spheres[3]->center->getP()[i] - csp.data()[i]);

		Vector2 rel_uk = lagBasis.transpose() * rel_u;
		qeal rel_ukSqNorm = rel_uk.squaredNorm();
		qeal rel_ukNorm = std::sqrt(rel_ukSqNorm);
		qeal f1_div_relDXNorm, f2_term;
		f1_SF_Div_RelDXNorm(rel_ukSqNorm, epsvh, f1_div_relDXNorm);
		f2_SF_Term(rel_ukSqNorm, epsvh, f2_term);

		Vector3 fricForce = -1.0 * f1_div_relDXNorm * mu *lagLamda * lagBasis * rel_uk;

		MatrixX HessianI(12, 12);
		HessianI.setZero();
		MatrixX selectMatrix(12, 3);
		selectMatrix.setZero();
		selectMatrix.data()[0] = lagAlpha;   selectMatrix.data()[13] = lagAlpha;   selectMatrix.data()[26] = lagAlpha;
		selectMatrix.data()[3] = lagBbeta;  selectMatrix.data()[16] = lagBbeta;  selectMatrix.data()[29] = lagBbeta;
		selectMatrix.data()[6] = 1.0 - lagAlpha - lagBbeta;   selectMatrix.data()[19] = 1.0 - lagAlpha - lagBbeta;   selectMatrix.data()[32] = 1.0 - lagAlpha - lagBbeta;
		selectMatrix.data()[9] = -1.0; selectMatrix.data()[22] = -1.0;  selectMatrix.data()[35] = -1.0;

		MatrixX TT = lagBasis.transpose() * selectMatrix.transpose();

		if (rel_ukSqNorm >= (epsvh * epsvh)) {
			// no SPD projection needed
			Vector2 ubar;
			ubar.data()[0] = -rel_uk.data()[1];
			ubar.data()[1] = rel_uk.data()[0];
			HessianI = (TT.transpose() * ((mu * lagLamda * f1_div_relDXNorm / rel_ukSqNorm) * ubar)) * (ubar.transpose() * TT);
		}
		else
		{
			if (rel_ukSqNorm == 0) {
				// no SPD projection needed
				HessianI = ((mu * lagLamda * f1_div_relDXNorm) * TT.transpose()) * TT;
			}
			else
			{
				// only need to project the inner 2x2 matrix to SPD
				MatrixX innerMtr = ((f2_term / rel_ukNorm) * rel_uk) * rel_uk.transpose();
				innerMtr.diagonal().array() += f1_div_relDXNorm;
				makePD(innerMtr);
				innerMtr *= mu * lagLamda;
				// tensor product:
				HessianI = TT.transpose() * innerMtr * TT;
			}
		}

		VectorX ff(12);
		ff.setZero();
		if (spheres[0]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[0]->center->getFrameId();

			ff.data()[0] = lagAlpha * fricForce.data()[0];
			ff.data()[1] = lagAlpha * fricForce.data()[1];
			ff.data()[2] = lagAlpha * fricForce.data()[2];
		}

		if (spheres[1]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[1]->center->getFrameId();

			ff.data()[3] = lagBbeta * fricForce.data()[0];
			ff.data()[4] = lagBbeta * fricForce.data()[1];
			ff.data()[5] = lagBbeta * fricForce.data()[2];
		}

		if (spheres[2]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[2]->center->getFrameId();

			ff.data()[6] = (1.0 - lagAlpha - lagBbeta) * fricForce.data()[0];
			ff.data()[7] = (1.0 - lagAlpha - lagBbeta) * fricForce.data()[1];
			ff.data()[8] = (1.0 - lagAlpha - lagBbeta) * fricForce.data()[2];
		}

		if (spheres[3]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[3]->center->getFrameId();

			ff.data()[9] = -1.0 * fricForce.data()[0];
			ff.data()[10] = -1.0 * fricForce.data()[1];
			ff.data()[11] = -1.0 * fricForce.data()[2];
		}
		fillOverallGradient(1.0, ff, gradient);
		fillOverallHessian(1.0, HessianI, hessian);
	}

	void MipcSlabSphereConstraint::getGradientAndHessian(qeal kappa, VectorX& gradient, std::vector<TripletX>& triplet)
	{
		qeal barrierGrad = getBarrierGradient();
		qeal barrierHessian = getBarrierHessian();
		VectorX distGrad;
		distGrad.resize(12);
		distGrad.setZero();
		MatrixX distHessina;
		distHessina.resize(12, 12);
		distHessina.setZero();
		diff_F_x(distGrad);

		fillOverallGradient(-1.0 * kappa * barrierGrad, distGrad, gradient);

		switch (distanceMode)
		{
		case TWO_ENDPOINTS:
			endPointsHessina(distHessina);
			break;
		case ALPHA_ZERO:
			alphaIsZeroHessina(distHessina);
			break;
		case BETA_ZERO:
			betaIsZeroHessina(distHessina);
			break;
		case ALPHA_BETA_ONE:
			alphaBetaPlusOneHessina(distHessina);
			break;
		case ALPHA_BETA:
			alphaBetaHessina(distHessina);
			break;
		default:
			break;
		};

		MatrixX hess = barrierGrad * distHessina + barrierHessian * (distGrad * distGrad.transpose());
		makePD(hess);

		hess *= kappa;
		fillOverallHessian(1.0, hess, triplet);
	}

	void MipcSlabSphereConstraint::getFrictionGradientAndHessian(qeal kappa, VectorX& gradient, std::vector<TripletX>& triplet)
	{
		Vector3 c11p, c12p, c13p, csp;
		spheres[0]->center->projectFullspacePreP(c11p.data());
		spheres[1]->center->projectFullspacePreP(c12p.data());
		spheres[2]->center->projectFullspacePreP(c13p.data());
		spheres[3]->center->projectFullspacePreP(csp.data());

		Vector3 rel_u;
		for (int i = 0; i < 3; i++)
			rel_u.data()[i] = (lagAlpha * (spheres[0]->center->getP()[i] - c11p.data()[i]) + lagBbeta * (spheres[1]->center->getP()[i] - c12p.data()[i]) + (1.0 - lagAlpha - lagBbeta) * (spheres[2]->center->getP()[i] - c13p.data()[i])) - (spheres[3]->center->getP()[i] - csp.data()[i]);

		Vector2 rel_uk = lagBasis.transpose() * rel_u;
		qeal rel_ukSqNorm = rel_uk.squaredNorm();
		qeal rel_ukNorm = std::sqrt(rel_ukSqNorm);
		qeal f1_div_relDXNorm, f2_term;
		f1_SF_Div_RelDXNorm(rel_ukSqNorm, epsvh, f1_div_relDXNorm);
		f2_SF_Term(rel_ukSqNorm, epsvh, f2_term);

		Vector3 fricForce = -1.0 * f1_div_relDXNorm * mu *lagLamda * lagBasis * rel_uk;

		MatrixX HessianI(12, 12);
		HessianI.setZero();
		MatrixX selectMatrix(12, 3);
		selectMatrix.setZero();
		selectMatrix.data()[0] = lagAlpha;   selectMatrix.data()[13] = lagAlpha;   selectMatrix.data()[26] = lagAlpha;
		selectMatrix.data()[3] = lagBbeta;  selectMatrix.data()[16] = lagBbeta;  selectMatrix.data()[29] = lagBbeta;
		selectMatrix.data()[6] = 1.0 - lagAlpha - lagBbeta;   selectMatrix.data()[19] = 1.0 - lagAlpha - lagBbeta;   selectMatrix.data()[32] = 1.0 - lagAlpha - lagBbeta;
		selectMatrix.data()[9] = -1.0; selectMatrix.data()[22] = -1.0;  selectMatrix.data()[35] = -1.0;

		MatrixX TT = lagBasis.transpose() * selectMatrix.transpose();

		if (rel_ukSqNorm >= (epsvh * epsvh)) {
			// no SPD projection needed
			Vector2 ubar;
			ubar.data()[0] = -rel_uk.data()[1];
			ubar.data()[1] = rel_uk.data()[0];
			HessianI = (TT.transpose() * ((mu * lagLamda * f1_div_relDXNorm / rel_ukSqNorm) * ubar)) * (ubar.transpose() * TT);
		}
		else
		{
			if (rel_ukSqNorm == 0) {
				// no SPD projection needed
				HessianI = ((mu * lagLamda * f1_div_relDXNorm) * TT.transpose()) * TT;
			}
			else
			{
				// only need to project the inner 2x2 matrix to SPD
				MatrixX innerMtr = ((f2_term / rel_ukNorm) * rel_uk) * rel_uk.transpose();
				innerMtr.diagonal().array() += f1_div_relDXNorm;
				makePD(innerMtr);
				innerMtr *= mu * lagLamda;
				// tensor product:
				HessianI = TT.transpose() * innerMtr * TT;
			}
		}

		VectorX ff(12);
		ff.setZero();
		if (spheres[0]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[0]->center->getFrameId();

			ff.data()[0] = lagAlpha * fricForce.data()[0];
			ff.data()[1] = lagAlpha * fricForce.data()[1];
			ff.data()[2] = lagAlpha * fricForce.data()[2];
		}

		if (spheres[1]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[1]->center->getFrameId();

			ff.data()[3] = lagBbeta * fricForce.data()[0];
			ff.data()[4] = lagBbeta * fricForce.data()[1];
			ff.data()[5] = lagBbeta * fricForce.data()[2];
		}

		if (spheres[2]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[2]->center->getFrameId();

			ff.data()[6] = (1.0 - lagAlpha - lagBbeta) * fricForce.data()[0];
			ff.data()[7] = (1.0 - lagAlpha - lagBbeta) * fricForce.data()[1];
			ff.data()[8] = (1.0 - lagAlpha - lagBbeta) * fricForce.data()[2];
		}

		if (spheres[3]->center->getFrameType() != FrameType::STATIC)
		{
			int frameId = spheres[3]->center->getFrameId();

			ff.data()[9] = -1.0 * fricForce.data()[0];
			ff.data()[10] = -1.0 * fricForce.data()[1];
			ff.data()[11] = -1.0 * fricForce.data()[2];
		}
		fillOverallGradient(1.0, ff, gradient);
		fillOverallHessian(1.0, HessianI, triplet);
	}

	void MipcSlabSphereConstraint::diff_F_x(VectorX & diff)
	{
		Vector3 v = 2.0 * (alpha * sC1 + beta * sC2 + sC3);

		diff[0] = (alpha * v.data()[0]);
		diff[1] = (alpha * v.data()[1]);
		diff[2] = (alpha * v.data()[2]);

		diff[3] = (beta * v.data()[0]);
		diff[4] = (beta * v.data()[1]);
		diff[5] = (beta * v.data()[2]);

		diff[6] = ((1.0 - alpha - beta) * v.data()[0]);
		diff[7] = ((1.0 - alpha - beta) * v.data()[1]);
		diff[8] = ((1.0 - alpha - beta) * v.data()[2]);

		diff[9] = (-v.data()[0]);
		diff[10] = (-v.data()[1]);
		diff[11] = (-v.data()[2]);
	}

	void MipcSlabSphereConstraint::endPointsHessina(MatrixX & hessina)
	{
		// alpha and beta is constant
		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (alpha * C1x + beta * C2x + C3x);
		qeal vy = 2.0 * (alpha * C1y + beta * C2y + C3y);
		qeal vz = 2.0 * (alpha * C1z + beta * C2z + C3z);

		qeal dvxdc11_x = 2.0 * alpha;
		qeal dvxdc11_y = 0;
		qeal dvxdc11_z = 0;

		qeal dvxdc12_x = 2.0 * beta;
		qeal dvxdc12_y = 0;
		qeal dvxdc12_z = 0;

		qeal dvxdc13_x = 2.0 * (1.0 - alpha - beta);
		qeal dvxdc13_y = 0;
		qeal dvxdc13_z = 0;

		qeal dvxdcp_x = 2.0 * -1.0;
		qeal dvxdcp_y = 0;
		qeal dvxdcp_z = 0;
		//
		qeal dvydc11_x = 0;
		qeal dvydc11_y = 2.0 * alpha;
		qeal dvydc11_z = 0;

		qeal dvydc12_x = 0;
		qeal dvydc12_y = 2.0 * beta;
		qeal dvydc12_z = 0;

		qeal dvydc13_x = 0;
		qeal dvydc13_y = 2.0 * (1.0 - alpha - beta);
		qeal dvydc13_z = 0;

		qeal dvydcp_x = 0;
		qeal dvydcp_y = 2.0 * -1.0;
		qeal dvydcp_z = 0;
		//
		qeal dvzdc11_x = 0;
		qeal dvzdc11_y = 0;
		qeal dvzdc11_z = 2.0 * alpha;

		qeal dvzdc12_x = 0;
		qeal dvzdc12_y = 0;
		qeal dvzdc12_z = 2.0 * beta;

		qeal dvzdc13_x = 0;
		qeal dvzdc13_y = 0;
		qeal dvzdc13_z = 2.0 * (1.0 - alpha - beta);

		qeal dvzdcp_x = 0;
		qeal dvzdcp_y = 0;
		qeal dvzdcp_z = 2.0 * -1.0;

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = alpha * dvxdc11_x;
		hessina.data()[1] = alpha * dvydc11_x;
		hessina.data()[2] = alpha * dvzdc11_x;

		hessina.data()[3] = beta * dvxdc11_x;
		hessina.data()[4] = beta * dvydc11_x;
		hessina.data()[5] = beta * dvzdc11_x;

		hessina.data()[6] = (1.0 - alpha - beta) * dvxdc11_x;
		hessina.data()[7] = (1.0 - alpha - beta) * dvydc11_x;
		hessina.data()[8] = (1.0 - alpha - beta) * dvzdc11_x;

		hessina.data()[9] = -dvxdc11_x;
		hessina.data()[10] = -dvydc11_x;
		hessina.data()[11] = -dvzdc11_x;
		// f / c11y
		hessina.data()[12] = alpha * dvxdc11_y;
		hessina.data()[13] = alpha * dvydc11_y;
		hessina.data()[14] = alpha * dvzdc11_y;

		hessina.data()[15] = beta * dvxdc11_y;
		hessina.data()[16] = beta * dvydc11_y;
		hessina.data()[17] = beta * dvzdc11_y;

		hessina.data()[18] = (1.0 - alpha - beta) * dvxdc11_y;
		hessina.data()[19] = (1.0 - alpha - beta) * dvydc11_y;
		hessina.data()[20] = (1.0 - alpha - beta) * dvzdc11_y;

		hessina.data()[21] = -dvxdc11_y;
		hessina.data()[22] = -dvydc11_y;
		hessina.data()[23] = -dvzdc11_y;
		// f / c11z
		hessina.data()[24] = alpha * dvxdc11_z;
		hessina.data()[25] = alpha * dvydc11_z;
		hessina.data()[26] = alpha * dvzdc11_z;

		hessina.data()[27] = beta * dvxdc11_z;
		hessina.data()[28] = beta * dvydc11_z;
		hessina.data()[29] = beta * dvzdc11_z;

		hessina.data()[30] = (1.0 - alpha - beta) * dvxdc11_z;
		hessina.data()[31] = (1.0 - alpha - beta) * dvydc11_z;
		hessina.data()[32] = (1.0 - alpha - beta) * dvzdc11_z;

		hessina.data()[33] = -dvxdc11_z;
		hessina.data()[34] = -dvydc11_z;
		hessina.data()[35] = -dvzdc11_z;
		// f / c12x
		hessina.data()[36] = alpha * dvxdc12_x;
		hessina.data()[37] = alpha * dvydc12_x;
		hessina.data()[38] = alpha * dvzdc12_x;

		hessina.data()[39] = beta * dvxdc12_x;
		hessina.data()[40] = beta * dvydc12_x;
		hessina.data()[41] = beta * dvzdc12_x;

		hessina.data()[42] = (1.0 - alpha - beta) * dvxdc12_x;
		hessina.data()[43] = (1.0 - alpha - beta) * dvydc12_x;
		hessina.data()[44] = (1.0 - alpha - beta) * dvzdc12_x;

		hessina.data()[45] = -dvxdc12_x;
		hessina.data()[46] = -dvydc12_x;
		hessina.data()[47] = -dvzdc12_x;
		// f / c12y
		hessina.data()[48] = alpha * dvxdc12_y;
		hessina.data()[49] = alpha * dvydc12_y;
		hessina.data()[50] = alpha * dvzdc12_y;

		hessina.data()[51] = beta * dvxdc12_y;
		hessina.data()[52] = beta * dvydc12_y;
		hessina.data()[53] = beta * dvzdc12_y;

		hessina.data()[54] = (1.0 - alpha - beta) * dvxdc12_y;
		hessina.data()[55] = (1.0 - alpha - beta) * dvydc12_y;
		hessina.data()[56] = (1.0 - alpha - beta) * dvzdc12_y;

		hessina.data()[57] = -dvxdc12_y;
		hessina.data()[58] = -dvydc12_y;
		hessina.data()[59] = -dvzdc12_y;
		// f / c12z
		hessina.data()[60] = alpha * dvxdc12_z;
		hessina.data()[61] = alpha * dvydc12_z;
		hessina.data()[62] = alpha * dvzdc12_z;

		hessina.data()[63] = beta * dvxdc12_z;
		hessina.data()[64] = beta * dvydc12_z;
		hessina.data()[65] = beta * dvzdc12_z;

		hessina.data()[66] = (1.0 - alpha - beta) * dvxdc12_z;
		hessina.data()[67] = (1.0 - alpha - beta) * dvydc12_z;
		hessina.data()[68] = (1.0 - alpha - beta) * dvzdc12_z;

		hessina.data()[69] = -dvxdc12_z;
		hessina.data()[70] = -dvydc12_z;
		hessina.data()[71] = -dvzdc12_z;
		//	f / c13x
		hessina.data()[72] = alpha * dvxdc13_x;
		hessina.data()[73] = alpha * dvydc13_x;
		hessina.data()[74] = alpha * dvzdc13_x;

		hessina.data()[75] = beta * dvxdc13_x;
		hessina.data()[76] = beta * dvydc13_x;
		hessina.data()[77] = beta * dvzdc13_x;

		hessina.data()[78] = (1.0 - alpha - beta) * dvxdc13_x;
		hessina.data()[79] = (1.0 - alpha - beta) * dvydc13_x;
		hessina.data()[80] = (1.0 - alpha - beta) * dvzdc13_x;

		hessina.data()[81] = -dvxdc13_x;
		hessina.data()[82] = -dvydc13_x;
		hessina.data()[83] = -dvzdc13_x;
		//	f / c13y
		hessina.data()[84] = alpha * dvxdc13_y;
		hessina.data()[85] = alpha * dvydc13_y;
		hessina.data()[86] = alpha * dvzdc13_y;

		hessina.data()[87] = beta * dvxdc13_y;
		hessina.data()[88] = beta * dvydc13_y;
		hessina.data()[89] = beta * dvzdc13_y;

		hessina.data()[90] = (1.0 - alpha - beta) * dvxdc13_y;
		hessina.data()[91] = (1.0 - alpha - beta) * dvydc13_y;
		hessina.data()[92] = (1.0 - alpha - beta) * dvzdc13_y;

		hessina.data()[93] = -dvxdc13_y;
		hessina.data()[94] = -dvydc13_y;
		hessina.data()[95] = -dvzdc13_y;
		//	f / c13z
		hessina.data()[96] = alpha * dvxdc13_z;
		hessina.data()[97] = alpha * dvydc13_z;
		hessina.data()[98] = alpha * dvzdc13_z;

		hessina.data()[99] = beta * dvxdc13_z;
		hessina.data()[100] = beta * dvydc13_z;
		hessina.data()[101] = beta * dvzdc13_z;

		hessina.data()[102] = (1.0 - alpha - beta) * dvxdc13_z;
		hessina.data()[103] = (1.0 - alpha - beta) * dvydc13_z;
		hessina.data()[104] = (1.0 - alpha - beta) * dvzdc13_z;

		hessina.data()[105] = -dvxdc13_z;
		hessina.data()[106] = -dvydc13_z;
		hessina.data()[107] = -dvzdc13_z;
		//	f / cpx
		hessina.data()[108] = alpha * dvxdcp_x;
		hessina.data()[109] = alpha * dvydcp_x;
		hessina.data()[110] = alpha * dvzdcp_x;

		hessina.data()[111] = beta * dvxdcp_x;
		hessina.data()[112] = beta * dvydcp_x;
		hessina.data()[113] = beta * dvzdcp_x;

		hessina.data()[114] = (1.0 - alpha - beta) * dvxdcp_x;
		hessina.data()[115] = (1.0 - alpha - beta) * dvydcp_x;
		hessina.data()[116] = (1.0 - alpha - beta) * dvzdcp_x;

		hessina.data()[117] = -dvxdcp_x;
		hessina.data()[118] = -dvydcp_x;
		hessina.data()[119] = -dvzdcp_x;
		//	f / cpy
		hessina.data()[120] = alpha * dvxdcp_y;
		hessina.data()[121] = alpha * dvydcp_y;
		hessina.data()[122] = alpha * dvzdcp_y;

		hessina.data()[123] = beta * dvxdcp_y;
		hessina.data()[124] = beta * dvydcp_y;
		hessina.data()[125] = beta * dvzdcp_y;

		hessina.data()[126] = (1.0 - alpha - beta) * dvxdcp_y;
		hessina.data()[127] = (1.0 - alpha - beta) * dvydcp_y;
		hessina.data()[128] = (1.0 - alpha - beta) * dvzdcp_y;

		hessina.data()[129] = -dvxdcp_y;
		hessina.data()[130] = -dvydcp_y;
		hessina.data()[131] = -dvzdcp_y;
		//	f / cpz
		hessina.data()[132] = alpha * dvxdcp_z;
		hessina.data()[133] = alpha * dvydcp_z;
		hessina.data()[134] = alpha * dvzdcp_z;

		hessina.data()[135] = beta * dvxdcp_z;
		hessina.data()[136] = beta * dvydcp_z;
		hessina.data()[137] = beta * dvzdcp_z;

		hessina.data()[138] = (1.0 - alpha - beta) * dvxdcp_z;
		hessina.data()[139] = (1.0 - alpha - beta) * dvydcp_z;
		hessina.data()[140] = (1.0 - alpha - beta) * dvzdcp_z;

		hessina.data()[141] = -dvxdcp_z;
		hessina.data()[142] = -dvydcp_z;
		hessina.data()[143] = -dvzdcp_z;
	}

	void MipcSlabSphereConstraint::alphaIsZeroHessina(MatrixX & hessina)
	{
		qeal dAdc11_x = 2.0 * sC1.data()[0]; qeal dAdc11_y = 2.0 * sC1.data()[1]; qeal dAdc11_z = 2.0 * sC1.data()[2];
		qeal dAdc13_x = -2.0 * sC1.data()[0]; qeal dAdc13_y = -2.0 * sC1.data()[1]; qeal dAdc13_z = -2.0 * sC1.data()[2];

		qeal dBdc11_x = 2.0 * sC2.data()[0]; qeal dBdc11_y = 2.0 * sC2.data()[1]; qeal dBdc11_z = 2.0 * sC2.data()[2];
		qeal dBdc12_x = 2.0 * sC1.data()[0]; qeal dBdc12_y = 2.0 * sC1.data()[1]; qeal dBdc12_z = 2.0 * sC1.data()[2];
		qeal dBdc13_x = -2.0 * (sC1.data()[0] + sC2.data()[0]); qeal dBdc13_y = -2.0 * (sC1.data()[1] + sC2.data()[1]); qeal dBdc13_z = -2.0 * (sC1.data()[2] + sC2.data()[2]);

		qeal dCdc12_x = 2.0 * sC2.data()[0]; qeal dCdc12_y = 2.0 * sC2.data()[1]; qeal dCdc12_z = 2.0 * sC2.data()[2];
		qeal dCdc13_x = -2.0 * sC2.data()[0]; qeal dCdc13_y = -2.0 * sC2.data()[1]; qeal dCdc13_z = -2.0 * sC2.data()[2];

		qeal dDdc11_x = 2.0 * sC3.data()[0]; qeal dDdc11_y = 2.0 * sC3.data()[1]; qeal dDdc11_z = 2.0 * sC3.data()[2];
		qeal dDdc13_x = 2.0 * (sC1.data()[0] - sC3.data()[0]); qeal dDdc13_y = 2.0 * (sC1.data()[1] - sC3.data()[1]); qeal dDdc13_z = 2.0 * (sC1.data()[2] - sC3.data()[2]);
		qeal dDdcp_x = -2.0 * sC1.data()[0]; qeal dDdcp_y = -2.0 * sC1.data()[1]; qeal dDdcp_z = -2.0 * sC1.data()[2];

		qeal dEdc12_x = 2.0 * sC3.data()[0]; qeal dEdc12_y = 2.0 * sC3.data()[1]; qeal dEdc12_z = 2.0 * sC3.data()[2];
		qeal dEdc13_x = 2.0 * (sC2.data()[0] - sC3.data()[0]); qeal dEdc13_y = 2.0 * (sC2.data()[1] - sC3.data()[1]); qeal dEdc13_z = 2.0 * (sC2.data()[2] - sC3.data()[2]);
		qeal dEdcp_x = -2.0 * sC2.data()[0]; qeal dEdcp_y = -2.0 * sC2.data()[1]; qeal dEdcp_z = -2.0 * sC2.data()[2];

		qeal dFdc13_x = 2.0 * sC3.data()[0]; qeal dFdc13_y = 2.0 * sC3.data()[1]; qeal dFdc13_z = 2.0 * sC3.data()[2];
		qeal dFdcp_x = -2.0 * sC3.data()[0]; qeal dFdcp_y = -2.0 * sC3.data()[1]; qeal dFdcp_z = -2.0 * sC3.data()[2];

		//
		qeal q1 = 0.5 * E / (C * C);  qeal q2 = -0.5 / C;

		qeal dBTdc12_x = q1 * dCdc12_x + q2 * dEdc12_x;
		qeal dBTdc12_y = q1 * dCdc12_y + q2 * dEdc12_y;
		qeal dBTdc12_z = q1 * dCdc12_z + q2 * dEdc12_z;

		qeal dBTdc13_x = q1 * dCdc13_x + q2 * dEdc13_x;
		qeal dBTdc13_y = q1 * dCdc13_y + q2 * dEdc13_y;
		qeal dBTdc13_z = q1 * dCdc13_z + q2 * dEdc13_z;

		qeal dBTdcp_x = q2 * dEdcp_x;
		qeal dBTdcp_y = q2 * dEdcp_y;
		qeal dBTdcp_z = q2 * dEdcp_z;

		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (beta * C2x + C3x);
		qeal vy = 2.0 * (beta * C2y + C3y);
		qeal vz = 2.0 * (beta * C2z + C3z);

		qeal dvxdc12_x = 2.0 * (C2x * dBTdc12_x + beta);
		qeal dvxdc12_y = 2.0 * (C2x * dBTdc12_y);
		qeal dvxdc12_z = 2.0 * (C2x * dBTdc12_z);

		qeal dvxdc13_x = 2.0 * (C2x * dBTdc13_x + 1 - beta);
		qeal dvxdc13_y = 2.0 * (C2x * dBTdc13_y);
		qeal dvxdc13_z = 2.0 * (C2x * dBTdc13_z);

		qeal dvxdcp_x = 2.0 * (C2x * dBTdcp_x - 1);
		qeal dvxdcp_y = 2.0 * (C2x * dBTdcp_y);
		qeal dvxdcp_z = 2.0 * (C2x * dBTdcp_z);
		//

		qeal dvydc12_x = 2.0 * (C2y * dBTdc12_x);
		qeal dvydc12_y = 2.0 * (C2y * dBTdc12_y + beta);
		qeal dvydc12_z = 2.0 * (C2y * dBTdc12_z);

		qeal dvydc13_x = 2.0 * (C2y * dBTdc13_x);
		qeal dvydc13_y = 2.0 * (C2y * dBTdc13_y + 1 - beta);
		qeal dvydc13_z = 2.0 * (C2y * dBTdc13_z);

		qeal dvydcp_x = 2.0 * (C2y * dBTdcp_x);
		qeal dvydcp_y = 2.0 * (C2y * dBTdcp_y - 1);
		qeal dvydcp_z = 2.0 * (C2y * dBTdcp_z);
		//

		qeal dvzdc12_x = 2.0 * (C2z * dBTdc12_x);
		qeal dvzdc12_y = 2.0 * (C2z * dBTdc12_y);
		qeal dvzdc12_z = 2.0 * (C2z * dBTdc12_z + beta);

		qeal dvzdc13_x = 2.0 * (C2z * dBTdc13_x);
		qeal dvzdc13_y = 2.0 * (C2z * dBTdc13_y);
		qeal dvzdc13_z = 2.0 * (C2z * dBTdc13_z + 1 - beta);

		qeal dvzdcp_x = 2.0 * (C2z * dBTdcp_x);
		qeal dvzdcp_y = 2.0 * (C2z * dBTdcp_y);
		qeal dvzdcp_z = 2.0 * (C2z * dBTdcp_z - 1);

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = 0;
		hessina.data()[1] = 0;
		hessina.data()[2] = 0;

		hessina.data()[3] = 0;
		hessina.data()[4] = 0;
		hessina.data()[5] = 0;

		hessina.data()[6] = 0;
		hessina.data()[7] = 0;
		hessina.data()[8] = 0;

		hessina.data()[9] = 0;
		hessina.data()[10] = 0;
		hessina.data()[11] = 0;
		// f / c11y
		hessina.data()[12] = 0;
		hessina.data()[13] = 0;
		hessina.data()[14] = 0;

		hessina.data()[15] = 0;
		hessina.data()[16] = 0;
		hessina.data()[17] = 0;

		hessina.data()[18] = 0;
		hessina.data()[19] = 0;
		hessina.data()[20] = 0;

		hessina.data()[21] = 0;
		hessina.data()[22] = 0;
		hessina.data()[23] = 0;
		// f / c11z
		hessina.data()[24] = 0;
		hessina.data()[25] = 0;
		hessina.data()[26] = 0;

		hessina.data()[27] = 0;
		hessina.data()[28] = 0;
		hessina.data()[29] = 0;

		hessina.data()[30] = 0;
		hessina.data()[31] = 0;
		hessina.data()[32] = 0;

		hessina.data()[33] = 0;
		hessina.data()[34] = 0;
		hessina.data()[35] = 0;
		// f / c12x
		hessina.data()[36] = 0;
		hessina.data()[37] = 0;
		hessina.data()[38] = 0;

		hessina.data()[39] = vx * dBTdc12_x + beta * dvxdc12_x;
		hessina.data()[40] = vy * dBTdc12_x + beta * dvydc12_x;  //
		hessina.data()[41] = vz * dBTdc12_x + beta * dvzdc12_x;

		hessina.data()[42] = vx * (-dBTdc12_x) + (1.0 - beta) * dvxdc12_x;
		hessina.data()[43] = vy * (-dBTdc12_x) + (1.0 - beta) * dvydc12_x;
		hessina.data()[44] = vz * (-dBTdc12_x) + (1.0 - beta) * dvzdc12_x;

		hessina.data()[45] = -dvxdc12_x;
		hessina.data()[46] = -dvydc12_x;
		hessina.data()[47] = -dvzdc12_x;
		// f / c12y
		hessina.data()[48] = 0;
		hessina.data()[49] = 0;
		hessina.data()[50] = 0;

		hessina.data()[51] = vx * dBTdc12_y + beta * dvxdc12_y;  //
		hessina.data()[52] = vy * dBTdc12_y + beta * dvydc12_y;
		hessina.data()[53] = vz * dBTdc12_y + beta * dvzdc12_y;

		hessina.data()[54] = vx * (-dBTdc12_y) + (1.0 - beta) * dvxdc12_y;
		hessina.data()[55] = vy * (-dBTdc12_y) + (1.0 - beta) * dvydc12_y;
		hessina.data()[56] = vz * (-dBTdc12_y) + (1.0 - beta) * dvzdc12_y;

		hessina.data()[57] = -dvxdc12_y;
		hessina.data()[58] = -dvydc12_y;
		hessina.data()[59] = -dvzdc12_y;
		// f / c12z
		hessina.data()[60] = 0;
		hessina.data()[61] = 0;
		hessina.data()[62] = 0;

		hessina.data()[63] = vx * dBTdc12_z + beta * dvxdc12_z;
		hessina.data()[64] = vy * dBTdc12_z + beta * dvydc12_z;
		hessina.data()[65] = vz * dBTdc12_z + beta * dvzdc12_z;

		hessina.data()[66] = vx * (-dBTdc12_z) + (1.0 - beta) * dvxdc12_z;
		hessina.data()[67] = vy * (-dBTdc12_z) + (1.0 - beta) * dvydc12_z;
		hessina.data()[68] = vz * (-dBTdc12_z) + (1.0 - beta) * dvzdc12_z;

		hessina.data()[69] = -dvxdc12_z;
		hessina.data()[70] = -dvydc12_z;
		hessina.data()[71] = -dvzdc12_z;
		//	f / c13x
		hessina.data()[72] = 0;
		hessina.data()[73] = 0;
		hessina.data()[74] = 0;

		hessina.data()[75] = vx * dBTdc13_x + beta * dvxdc13_x;
		hessina.data()[76] = vy * dBTdc13_x + beta * dvydc13_x;
		hessina.data()[77] = vz * dBTdc13_x + beta * dvzdc13_x;

		hessina.data()[78] = vx * (-dBTdc13_x) + (1.0 - beta) * dvxdc13_x;
		hessina.data()[79] = vy * (-dBTdc13_x) + (1.0 - beta) * dvydc13_x;
		hessina.data()[80] = vz * (-dBTdc13_x) + (1.0 - beta) * dvzdc13_x;

		hessina.data()[81] = -dvxdc13_x;
		hessina.data()[82] = -dvydc13_x;
		hessina.data()[83] = -dvzdc13_x;
		//	f / c13y
		hessina.data()[84] = 0;
		hessina.data()[85] = 0;
		hessina.data()[86] = 0;

		hessina.data()[87] = vx * dBTdc13_y + beta * dvxdc13_y;
		hessina.data()[88] = vy * dBTdc13_y + beta * dvydc13_y;
		hessina.data()[89] = vz * dBTdc13_y + beta * dvzdc13_y;

		hessina.data()[90] = vx * (-dBTdc13_y) + (1.0 - beta) * dvxdc13_y;
		hessina.data()[91] = vy * (-dBTdc13_y) + (1.0 - beta) * dvydc13_y;
		hessina.data()[92] = vz * (-dBTdc13_y) + (1.0 - beta) * dvzdc13_y;

		hessina.data()[93] = -dvxdc13_y;
		hessina.data()[94] = -dvydc13_y;
		hessina.data()[95] = -dvzdc13_y;
		//	f / c13z
		hessina.data()[96] = 0;
		hessina.data()[97] = 0;
		hessina.data()[98] = 0;

		hessina.data()[99] = vx * dBTdc13_z + beta * dvxdc13_z;
		hessina.data()[100] = vy * dBTdc13_z + beta * dvydc13_z;
		hessina.data()[101] = vz * dBTdc13_z + beta * dvzdc13_z;

		hessina.data()[102] = vx * (-dBTdc13_z) + (1.0 - beta) * dvxdc13_z;
		hessina.data()[103] = vy * (-dBTdc13_z) + (1.0 - beta) * dvydc13_z;
		hessina.data()[104] = vz * (-dBTdc13_z) + (1.0 - beta) * dvzdc13_z;

		hessina.data()[105] = -dvxdc13_z;
		hessina.data()[106] = -dvydc13_z;
		hessina.data()[107] = -dvzdc13_z;
		//	f / cpx
		hessina.data()[108] = 0;
		hessina.data()[109] = 0;
		hessina.data()[110] = 0;

		hessina.data()[111] = vx * dBTdcp_x + beta * dvxdcp_x;
		hessina.data()[112] = vy * dBTdcp_x + beta * dvydcp_x;
		hessina.data()[113] = vz * dBTdcp_x + beta * dvzdcp_x;

		hessina.data()[114] = vx * (-dBTdcp_x) + (1.0 - beta) * dvxdcp_x;
		hessina.data()[115] = vy * (-dBTdcp_x) + (1.0 - beta) * dvydcp_x;
		hessina.data()[116] = vz * (-dBTdcp_x) + (1.0 - beta) * dvzdcp_x;

		hessina.data()[117] = -dvxdcp_x;
		hessina.data()[118] = -dvydcp_x;
		hessina.data()[119] = -dvzdcp_x;
		//	f / cpy
		hessina.data()[120] = 0;
		hessina.data()[121] = 0;
		hessina.data()[122] = 0;

		hessina.data()[123] = vx * dBTdcp_y + beta * dvxdcp_y;
		hessina.data()[124] = vy * dBTdcp_y + beta * dvydcp_y;
		hessina.data()[125] = vz * dBTdcp_y + beta * dvzdcp_y;

		hessina.data()[126] = vx * (-dBTdcp_y) + (1.0 - beta) * dvxdcp_y;
		hessina.data()[127] = vy * (-dBTdcp_y) + (1.0 - beta) * dvydcp_y;
		hessina.data()[128] = vz * (-dBTdcp_y) + (1.0 - beta) * dvzdcp_y;

		hessina.data()[129] = -dvxdcp_y;
		hessina.data()[130] = -dvydcp_y;
		hessina.data()[131] = -dvzdcp_y;
		//	f / cpz
		hessina.data()[132] = 0;
		hessina.data()[133] = 0;
		hessina.data()[134] = 0;

		hessina.data()[135] = vx * dBTdcp_z + beta * dvxdcp_z;
		hessina.data()[136] = vy * dBTdcp_z + beta * dvydcp_z;
		hessina.data()[137] = vz * dBTdcp_z + beta * dvzdcp_z;

		hessina.data()[138] = vx * (-dBTdcp_z) + (1.0 - beta) * dvxdcp_z;
		hessina.data()[139] = vy * (-dBTdcp_z) + (1.0 - beta) * dvydcp_z;
		hessina.data()[140] = vz * (-dBTdcp_z) + (1.0 - beta) * dvzdcp_z;

		hessina.data()[141] = -dvxdcp_z;
		hessina.data()[142] = -dvydcp_z;
		hessina.data()[143] = -dvzdcp_z;
	}

	void MipcSlabSphereConstraint::betaIsZeroHessina(MatrixX & hessina)
	{
		qeal dAdc11_x = 2.0 * sC1.data()[0]; qeal dAdc11_y = 2.0 * sC1.data()[1]; qeal dAdc11_z = 2.0 * sC1.data()[2];
		qeal dAdc13_x = -2.0 * sC1.data()[0]; qeal dAdc13_y = -2.0 * sC1.data()[1]; qeal dAdc13_z = -2.0 * sC1.data()[2];

		qeal dBdc11_x = 2.0 * sC2.data()[0]; qeal dBdc11_y = 2.0 * sC2.data()[1]; qeal dBdc11_z = 2.0 * sC2.data()[2];
		qeal dBdc12_x = 2.0 * sC1.data()[0]; qeal dBdc12_y = 2.0 * sC1.data()[1]; qeal dBdc12_z = 2.0 * sC1.data()[2];
		qeal dBdc13_x = -2.0 * (sC1.data()[0] + sC2.data()[0]); qeal dBdc13_y = -2.0 * (sC1.data()[1] + sC2.data()[1]); qeal dBdc13_z = -2.0 * (sC1.data()[2] + sC2.data()[2]);

		qeal dCdc12_x = 2.0 * sC2.data()[0]; qeal dCdc12_y = 2.0 * sC2.data()[1]; qeal dCdc12_z = 2.0 * sC2.data()[2];
		qeal dCdc13_x = -2.0 * sC2.data()[0]; qeal dCdc13_y = -2.0 * sC2.data()[1]; qeal dCdc13_z = -2.0 * sC2.data()[2];

		qeal dDdc11_x = 2.0 * sC3.data()[0]; qeal dDdc11_y = 2.0 * sC3.data()[1]; qeal dDdc11_z = 2.0 * sC3.data()[2];
		qeal dDdc13_x = 2.0 * (sC1.data()[0] - sC3.data()[0]); qeal dDdc13_y = 2.0 * (sC1.data()[1] - sC3.data()[1]); qeal dDdc13_z = 2.0 * (sC1.data()[2] - sC3.data()[2]);
		qeal dDdcp_x = -2.0 * sC1.data()[0]; qeal dDdcp_y = -2.0 * sC1.data()[1]; qeal dDdcp_z = -2.0 * sC1.data()[2];

		qeal dEdc12_x = 2.0 * sC3.data()[0]; qeal dEdc12_y = 2.0 * sC3.data()[1]; qeal dEdc12_z = 2.0 * sC3.data()[2];
		qeal dEdc13_x = 2.0 * (sC2.data()[0] - sC3.data()[0]); qeal dEdc13_y = 2.0 * (sC2.data()[1] - sC3.data()[1]); qeal dEdc13_z = 2.0 * (sC2.data()[2] - sC3.data()[2]);
		qeal dEdcp_x = -2.0 * sC2.data()[0]; qeal dEdcp_y = -2.0 * sC2.data()[1]; qeal dEdcp_z = -2.0 * sC2.data()[2];

		qeal dFdc13_x = 2.0 * sC3.data()[0]; qeal dFdc13_y = 2.0 * sC3.data()[1]; qeal dFdc13_z = 2.0 * sC3.data()[2];
		qeal dFdcp_x = -2.0 * sC3.data()[0]; qeal dFdcp_y = -2.0 * sC3.data()[1]; qeal dFdcp_z = -2.0 * sC3.data()[2];

		qeal q1 = 0.5 * D / (A * A);  qeal q2 = -0.5 / A;

		qeal dALdc11_x = q1 * dAdc11_x + q2 * dDdc11_x;
		qeal dALdc11_y = q1 * dAdc11_y + q2 * dDdc11_y;
		qeal dALdc11_z = q1 * dAdc11_z + q2 * dDdc11_z;

		qeal dALdc13_x = q1 * dAdc13_x + q2 * dDdc13_x;
		qeal dALdc13_y = q1 * dAdc13_y + q2 * dDdc13_y;
		qeal dALdc13_z = q1 * dAdc13_z + q2 * dDdc13_z;

		qeal dALdcp_x = q2 * dDdcp_x;
		qeal dALdcp_y = q2 * dDdcp_y;
		qeal dALdcp_z = q2 * dDdcp_z;

		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (alpha * C1x + C3x);
		qeal vy = 2.0 * (alpha * C1y + C3y);
		qeal vz = 2.0 * (alpha * C1z + C3z);

		qeal dvxdc11_x = 2.0 * (C1x * dALdc11_x + alpha);
		qeal dvxdc11_y = 2.0 * (C1x * dALdc11_y);
		qeal dvxdc11_z = 2.0 * (C1x * dALdc11_z);

		qeal dvxdc13_x = 2.0 * (C1x * dALdc13_x + 1 - alpha);
		qeal dvxdc13_y = 2.0 * (C1x * dALdc13_y);
		qeal dvxdc13_z = 2.0 * (C1x * dALdc13_z);

		qeal dvxdcp_x = 2.0 * (C1x * dALdcp_x - 1);
		qeal dvxdcp_y = 2.0 * (C1x * dALdcp_y);
		qeal dvxdcp_z = 2.0 * (C1x * dALdcp_z);
		//
		qeal dvydc11_x = 2.0 * (C1y * dALdc11_x);
		qeal dvydc11_y = 2.0 * (C1y * dALdc11_y + alpha);
		qeal dvydc11_z = 2.0 * (C1y * dALdc11_z);

		qeal dvydc13_x = 2.0 * (C1y * dALdc13_x);
		qeal dvydc13_y = 2.0 * (C1y * dALdc13_y + 1 - alpha);
		qeal dvydc13_z = 2.0 * (C1y * dALdc13_z);

		qeal dvydcp_x = 2.0 * (C1y * dALdcp_x);
		qeal dvydcp_y = 2.0 * (C1y * dALdcp_y - 1);
		qeal dvydcp_z = 2.0 * (C1y * dALdcp_z);
		//
		qeal dvzdc11_x = 2.0 * (C1z * dALdc11_x);
		qeal dvzdc11_y = 2.0 * (C1z * dALdc11_y);
		qeal dvzdc11_z = 2.0 * (C1z * dALdc11_z + alpha);

		qeal dvzdc13_x = 2.0 * (C1z * dALdc13_x);
		qeal dvzdc13_y = 2.0 * (C1z * dALdc13_y);
		qeal dvzdc13_z = 2.0 * (C1z * dALdc13_z + 1 - alpha);

		qeal dvzdcp_x = 2.0 * (C1z * dALdcp_x);
		qeal dvzdcp_y = 2.0 * (C1z * dALdcp_y);
		qeal dvzdcp_z = 2.0 * (C1z * dALdcp_z - 1);

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = Check_QEAL_ZERO(vx * dALdc11_x + alpha * dvxdc11_x);
		hessina.data()[1] = Check_QEAL_ZERO(vy * dALdc11_x + alpha * dvydc11_x);
		hessina.data()[2] = Check_QEAL_ZERO(vz * dALdc11_x + alpha * dvzdc11_x);

		hessina.data()[3] = 0;
		hessina.data()[4] = 0;
		hessina.data()[5] = 0;

		hessina.data()[6] = Check_QEAL_ZERO(vx * (-dALdc11_x) + (1.0 - alpha) * dvxdc11_x);
		hessina.data()[7] = Check_QEAL_ZERO(vy * (-dALdc11_x) + (1.0 - alpha) * dvydc11_x);
		hessina.data()[8] = Check_QEAL_ZERO(vz * (-dALdc11_x) + (1.0 - alpha) * dvzdc11_x);

		hessina.data()[9] = Check_QEAL_ZERO(-dvxdc11_x);
		hessina.data()[10] = Check_QEAL_ZERO(-dvydc11_x);
		hessina.data()[11] = Check_QEAL_ZERO(-dvzdc11_x);
		// f / c11y
		hessina.data()[12] = Check_QEAL_ZERO(vx * dALdc11_y + alpha * dvxdc11_y);
		hessina.data()[13] = Check_QEAL_ZERO(vy * dALdc11_y + alpha * dvydc11_y);
		hessina.data()[14] = Check_QEAL_ZERO(vz * dALdc11_y + alpha * dvzdc11_y);

		hessina.data()[15] = 0;
		hessina.data()[16] = 0;
		hessina.data()[17] = 0;

		hessina.data()[18] = Check_QEAL_ZERO(vx * (-dALdc11_y) + (1.0 - alpha) * dvxdc11_y);
		hessina.data()[19] = Check_QEAL_ZERO(vy * (-dALdc11_y) + (1.0 - alpha) * dvydc11_y);
		hessina.data()[20] = Check_QEAL_ZERO(vz * (-dALdc11_y) + (1.0 - alpha) * dvzdc11_y);

		hessina.data()[21] = Check_QEAL_ZERO(-dvxdc11_y);
		hessina.data()[22] = Check_QEAL_ZERO(-dvydc11_y);
		hessina.data()[23] = Check_QEAL_ZERO(-dvzdc11_y);
		// f / c11z
		hessina.data()[24] = Check_QEAL_ZERO(vx * dALdc11_z + alpha * dvxdc11_z);
		hessina.data()[25] = Check_QEAL_ZERO(vy * dALdc11_z + alpha * dvydc11_z);
		hessina.data()[26] = Check_QEAL_ZERO(vz * dALdc11_z + alpha * dvzdc11_z);

		hessina.data()[27] = 0;
		hessina.data()[28] = 0;
		hessina.data()[29] = 0;

		hessina.data()[30] = Check_QEAL_ZERO(vx * (-dALdc11_z) + (1.0 - alpha) * dvxdc11_z);
		hessina.data()[31] = Check_QEAL_ZERO(vy * (-dALdc11_z) + (1.0 - alpha) * dvydc11_z);
		hessina.data()[32] = Check_QEAL_ZERO(vz * (-dALdc11_z) + (1.0 - alpha) * dvzdc11_z);

		hessina.data()[33] = Check_QEAL_ZERO(-dvxdc11_z);
		hessina.data()[34] = Check_QEAL_ZERO(-dvydc11_z);
		hessina.data()[35] = Check_QEAL_ZERO(-dvzdc11_z);
		// f / c12x
		hessina.data()[36] = 0;
		hessina.data()[37] = 0;
		hessina.data()[38] = 0;

		hessina.data()[39] = 0;
		hessina.data()[40] = 0;
		hessina.data()[41] = 0;

		hessina.data()[42] = 0;
		hessina.data()[43] = 0;
		hessina.data()[44] = 0;

		hessina.data()[45] = 0;
		hessina.data()[46] = -0;
		hessina.data()[47] = 0;
		// f / c12y
		hessina.data()[48] = 0;
		hessina.data()[49] = 0;
		hessina.data()[50] = 0;

		hessina.data()[51] = 0;
		hessina.data()[52] = 0;
		hessina.data()[53] = 0;

		hessina.data()[54] = 0;
		hessina.data()[55] = 0;
		hessina.data()[56] = 0;

		hessina.data()[57] = 0;
		hessina.data()[58] = 0;
		hessina.data()[59] = 0;
		// f / c12z
		hessina.data()[60] = 0;
		hessina.data()[61] = 0;
		hessina.data()[62] = 0;

		hessina.data()[63] = 0;
		hessina.data()[64] = 0;
		hessina.data()[65] = 0;

		hessina.data()[66] = 0;
		hessina.data()[67] = 0;
		hessina.data()[68] = 0;

		hessina.data()[69] = 0;
		hessina.data()[70] = 0;
		hessina.data()[71] = 0;
		//	f / c13x
		hessina.data()[72] = Check_QEAL_ZERO(vx * dALdc13_x + alpha * dvxdc13_x);
		hessina.data()[73] = Check_QEAL_ZERO(vy * dALdc13_x + alpha * dvydc13_x);
		hessina.data()[74] = Check_QEAL_ZERO(vz * dALdc13_x + alpha * dvzdc13_x);

		hessina.data()[75] = 0;
		hessina.data()[76] = 0;
		hessina.data()[77] = 0;

		hessina.data()[78] = Check_QEAL_ZERO(vx * (-dALdc13_x) + (1.0 - alpha) * dvxdc13_x);
		hessina.data()[79] = Check_QEAL_ZERO(vy * (-dALdc13_x) + (1.0 - alpha) * dvydc13_x);
		hessina.data()[80] = Check_QEAL_ZERO(vz * (-dALdc13_x) + (1.0 - alpha) * dvzdc13_x);

		hessina.data()[81] = Check_QEAL_ZERO(-dvxdc13_x);
		hessina.data()[82] = Check_QEAL_ZERO(-dvydc13_x);
		hessina.data()[83] = Check_QEAL_ZERO(-dvzdc13_x);
		//	f / c13y
		hessina.data()[84] = Check_QEAL_ZERO(vx * dALdc13_y + alpha * dvxdc13_y);
		hessina.data()[85] = Check_QEAL_ZERO(vy * dALdc13_y + alpha * dvydc13_y);
		hessina.data()[86] = Check_QEAL_ZERO(vz * dALdc13_y + alpha * dvzdc13_y);

		hessina.data()[87] = Check_QEAL_ZERO(beta * dvxdc13_y);
		hessina.data()[88] = Check_QEAL_ZERO(beta * dvydc13_y);
		hessina.data()[89] = Check_QEAL_ZERO(beta * dvzdc13_y);

		hessina.data()[90] = Check_QEAL_ZERO(vx * (-dALdc13_y) + (1.0 - alpha) * dvxdc13_y);
		hessina.data()[91] = Check_QEAL_ZERO(vy * (-dALdc13_y) + (1.0 - alpha) * dvydc13_y);
		hessina.data()[92] = Check_QEAL_ZERO(vz * (-dALdc13_y) + (1.0 - alpha) * dvzdc13_y);

		hessina.data()[93] = Check_QEAL_ZERO(-dvxdc13_y);
		hessina.data()[94] = Check_QEAL_ZERO(-dvydc13_y);
		hessina.data()[95] = Check_QEAL_ZERO(-dvzdc13_y);
		//	f / c13z
		hessina.data()[96] = Check_QEAL_ZERO(vx * dALdc13_z + alpha * dvxdc13_z);
		hessina.data()[97] = Check_QEAL_ZERO(vy * dALdc13_z + alpha * dvydc13_z);
		hessina.data()[98] = Check_QEAL_ZERO(vz * dALdc13_z + alpha * dvzdc13_z);

		hessina.data()[99] = Check_QEAL_ZERO(beta * dvxdc13_z);
		hessina.data()[100] = Check_QEAL_ZERO(beta * dvydc13_z);
		hessina.data()[101] = Check_QEAL_ZERO(beta * dvzdc13_z);

		hessina.data()[102] = Check_QEAL_ZERO(vx * (-dALdc13_z) + (1.0 - alpha) * dvxdc13_z);
		hessina.data()[103] = Check_QEAL_ZERO(vy * (-dALdc13_z) + (1.0 - alpha) * dvydc13_z);
		hessina.data()[104] = Check_QEAL_ZERO(vz * (-dALdc13_z) + (1.0 - alpha) * dvzdc13_z);

		hessina.data()[105] = Check_QEAL_ZERO(-dvxdc13_z);
		hessina.data()[106] = Check_QEAL_ZERO(-dvydc13_z);
		hessina.data()[107] = Check_QEAL_ZERO(-dvzdc13_z);
		//	f / cpx
		hessina.data()[108] = Check_QEAL_ZERO(vx * dALdcp_x + alpha * dvxdcp_x);
		hessina.data()[109] = Check_QEAL_ZERO(vy * dALdcp_x + alpha * dvydcp_x);
		hessina.data()[110] = Check_QEAL_ZERO(vz * dALdcp_x + alpha * dvzdcp_x);

		hessina.data()[111] = Check_QEAL_ZERO(beta * dvxdcp_x);
		hessina.data()[112] = Check_QEAL_ZERO(beta * dvydcp_x);
		hessina.data()[113] = Check_QEAL_ZERO(beta * dvzdcp_x);

		hessina.data()[114] = Check_QEAL_ZERO(vx * (-dALdcp_x) + (1.0 - alpha) * dvxdcp_x);
		hessina.data()[115] = Check_QEAL_ZERO(vy * (-dALdcp_x) + (1.0 - alpha) * dvydcp_x);
		hessina.data()[116] = Check_QEAL_ZERO(vz * (-dALdcp_x) + (1.0 - alpha) * dvzdcp_x);

		hessina.data()[117] = Check_QEAL_ZERO(-dvxdcp_x);
		hessina.data()[118] = Check_QEAL_ZERO(-dvydcp_x);
		hessina.data()[119] = Check_QEAL_ZERO(-dvzdcp_x);
		//	f / cpy
		hessina.data()[120] = Check_QEAL_ZERO(vx * dALdcp_y + alpha * dvxdcp_y);
		hessina.data()[121] = Check_QEAL_ZERO(vy * dALdcp_y + alpha * dvydcp_y);
		hessina.data()[122] = Check_QEAL_ZERO(vz * dALdcp_y + alpha * dvzdcp_y);

		hessina.data()[123] = Check_QEAL_ZERO(beta * dvxdcp_y);
		hessina.data()[124] = Check_QEAL_ZERO(beta * dvydcp_y);
		hessina.data()[125] = Check_QEAL_ZERO(beta * dvzdcp_y);

		hessina.data()[126] = Check_QEAL_ZERO(vx * (-dALdcp_y) + (1.0 - alpha) * dvxdcp_y);
		hessina.data()[127] = Check_QEAL_ZERO(vy * (-dALdcp_y) + (1.0 - alpha) * dvydcp_y);
		hessina.data()[128] = Check_QEAL_ZERO(vz * (-dALdcp_y) + (1.0 - alpha) * dvzdcp_y);

		hessina.data()[129] = Check_QEAL_ZERO(-dvxdcp_y);
		hessina.data()[130] = Check_QEAL_ZERO(-dvydcp_y);
		hessina.data()[131] = Check_QEAL_ZERO(-dvzdcp_y);
		//	f / cpz
		hessina.data()[132] = Check_QEAL_ZERO(vx * dALdcp_z + alpha * dvxdcp_z);
		hessina.data()[133] = Check_QEAL_ZERO(vy * dALdcp_z + alpha * dvydcp_z);
		hessina.data()[134] = Check_QEAL_ZERO(vz * dALdcp_z + alpha * dvzdcp_z);

		hessina.data()[135] = Check_QEAL_ZERO(beta * dvxdcp_z);
		hessina.data()[136] = Check_QEAL_ZERO(beta * dvydcp_z);
		hessina.data()[137] = Check_QEAL_ZERO(beta * dvzdcp_z);

		hessina.data()[138] = Check_QEAL_ZERO(vx * (-dALdcp_z) + (1.0 - alpha) * dvxdcp_z);
		hessina.data()[139] = Check_QEAL_ZERO(vy * (-dALdcp_z) + (1.0 - alpha) * dvydcp_z);
		hessina.data()[140] = Check_QEAL_ZERO(vz * (-dALdcp_z) + (1.0 - alpha) * dvzdcp_z);

		hessina.data()[141] = Check_QEAL_ZERO(-dvxdcp_z);
		hessina.data()[142] = Check_QEAL_ZERO(-dvydcp_z);
		hessina.data()[143] = Check_QEAL_ZERO(-dvzdcp_z);
	}

	void MipcSlabSphereConstraint::alphaBetaPlusOneHessina(MatrixX & hessina)
	{
		qeal dAdc11_x = 2.0 * sC1.data()[0]; qeal dAdc11_y = 2.0 * sC1.data()[1]; qeal dAdc11_z = 2.0 * sC1.data()[2];
		qeal dAdc13_x = -2.0 * sC1.data()[0]; qeal dAdc13_y = -2.0 * sC1.data()[1]; qeal dAdc13_z = -2.0 * sC1.data()[2];

		qeal dBdc11_x = 2.0 * sC2.data()[0]; qeal dBdc11_y = 2.0 * sC2.data()[1]; qeal dBdc11_z = 2.0 * sC2.data()[2];
		qeal dBdc12_x = 2.0 * sC1.data()[0]; qeal dBdc12_y = 2.0 * sC1.data()[1]; qeal dBdc12_z = 2.0 * sC1.data()[2];
		qeal dBdc13_x = -2.0 * (sC1.data()[0] + sC2.data()[0]); qeal dBdc13_y = -2.0 * (sC1.data()[1] + sC2.data()[1]); qeal dBdc13_z = -2.0 * (sC1.data()[2] + sC2.data()[2]);

		qeal dCdc12_x = 2.0 * sC2.data()[0]; qeal dCdc12_y = 2.0 * sC2.data()[1]; qeal dCdc12_z = 2.0 * sC2.data()[2];
		qeal dCdc13_x = -2.0 * sC2.data()[0]; qeal dCdc13_y = -2.0 * sC2.data()[1]; qeal dCdc13_z = -2.0 * sC2.data()[2];

		qeal dDdc11_x = 2.0 * sC3.data()[0]; qeal dDdc11_y = 2.0 * sC3.data()[1]; qeal dDdc11_z = 2.0 * sC3.data()[2];
		qeal dDdc13_x = 2.0 * (sC1.data()[0] - sC3.data()[0]); qeal dDdc13_y = 2.0 * (sC1.data()[1] - sC3.data()[1]); qeal dDdc13_z = 2.0 * (sC1.data()[2] - sC3.data()[2]);
		qeal dDdcp_x = -2.0 * sC1.data()[0]; qeal dDdcp_y = -2.0 * sC1.data()[1]; qeal dDdcp_z = -2.0 * sC1.data()[2];

		qeal dEdc12_x = 2.0 * sC3.data()[0]; qeal dEdc12_y = 2.0 * sC3.data()[1]; qeal dEdc12_z = 2.0 * sC3.data()[2];
		qeal dEdc13_x = 2.0 * (sC2.data()[0] - sC3.data()[0]); qeal dEdc13_y = 2.0 * (sC2.data()[1] - sC3.data()[1]); qeal dEdc13_z = 2.0 * (sC2.data()[2] - sC3.data()[2]);
		qeal dEdcp_x = -2.0 * sC2.data()[0]; qeal dEdcp_y = -2.0 * sC2.data()[1]; qeal dEdcp_z = -2.0 * sC2.data()[2];

		qeal dFdc13_x = 2.0 * sC3.data()[0]; qeal dFdc13_y = 2.0 * sC3.data()[1]; qeal dFdc13_z = 2.0 * sC3.data()[2];
		qeal dFdcp_x = -2.0 * sC3.data()[0]; qeal dFdcp_y = -2.0 * sC3.data()[1]; qeal dFdcp_z = -2.0 * sC3.data()[2];

		//
		qeal q1 = 0.5 / (A + C - B); qeal q2 = alpha / (A + C - B);

		qeal dALdc11_x = q1 * (-(dBdc11_x + dDdc11_x)) - q2 * (dAdc11_x - dBdc11_x);
		qeal dALdc11_y = q1 * (-(dBdc11_y + dDdc11_y)) - q2 * (dAdc11_y - dBdc11_y);
		qeal dALdc11_z = q1 * (-(dBdc11_z + dDdc11_z)) - q2 * (dAdc11_z - dBdc11_z);

		qeal dALdc12_x = q1 * ((2.0 * dCdc12_x + dEdc12_x) - (dBdc12_x)) - q2 * (dCdc12_x - dBdc12_x);
		qeal dALdc12_y = q1 * ((2.0 * dCdc12_y + dEdc12_y) - (dBdc12_y)) - q2 * (dCdc12_y - dBdc12_y);
		qeal dALdc12_z = q1 * ((2.0 * dCdc12_z + dEdc12_z) - (dBdc12_z)) - q2 * (dCdc12_z - dBdc12_z);

		qeal dALdc13_x = 0;
		qeal dALdc13_y = 0;
		qeal dALdc13_z = 0;

		qeal dALdcp_x = q1 * (dEdcp_x - dDdcp_x);
		qeal dALdcp_y = q1 * (dEdcp_y - dDdcp_y);
		qeal dALdcp_z = q1 * (dEdcp_z - dDdcp_z);

		qeal dBTdc11_x = -dALdc11_x;
		qeal dBTdc11_y = -dALdc11_y;
		qeal dBTdc11_z = -dALdc11_z;

		qeal dBTdc12_x = -dALdc12_x;
		qeal dBTdc12_y = -dALdc12_y;
		qeal dBTdc12_z = -dALdc12_z;

		qeal dBTdc13_x = 0;
		qeal dBTdc13_y = 0;
		qeal dBTdc13_z = 0;

		qeal dBTdcp_x = -dALdcp_x;
		qeal dBTdcp_y = -dALdcp_y;
		qeal dBTdcp_z = -dALdcp_z;

		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (alpha * C1x + beta * C2x + C3x);
		qeal vy = 2.0 * (alpha * C1y + beta * C2y + C3y);
		qeal vz = 2.0 * (alpha * C1z + beta * C2z + C3z);

		qeal dvxdc11_x = 2.0 * (C1x * dALdc11_x + C2x * dBTdc11_x + alpha);
		qeal dvxdc11_y = 2.0 * (C1x * dALdc11_y + C2x * dBTdc11_y);
		qeal dvxdc11_z = 2.0 * (C1x * dALdc11_z + C2x * dBTdc11_z);

		qeal dvxdc12_x = 2.0 * (C1x * dALdc12_x + C2x * dBTdc12_x + beta);
		qeal dvxdc12_y = 2.0 * (C1x * dALdc12_y + C2x * dBTdc12_y);
		qeal dvxdc12_z = 2.0 * (C1x * dALdc12_z + C2x * dBTdc12_z);

		qeal dvxdc13_x = 0;
		qeal dvxdc13_y = 0;
		qeal dvxdc13_z = 0;

		qeal dvxdcp_x = 2.0 * (C1x * dALdcp_x + C2x * dBTdcp_x - 1);
		qeal dvxdcp_y = 2.0 * (C1x * dALdcp_y + C2x * dBTdcp_y);
		qeal dvxdcp_z = 2.0 * (C1x * dALdcp_z + C2x * dBTdcp_z);
		//
		qeal dvydc11_x = 2.0 * (C1y * dALdc11_x + C2y * dBTdc11_x);
		qeal dvydc11_y = 2.0 * (C1y * dALdc11_y + C2y * dBTdc11_y + alpha);
		qeal dvydc11_z = 2.0 * (C1y * dALdc11_z + C2y * dBTdc11_z);

		qeal dvydc12_x = 2.0 * (C1y * dALdc12_x + C2y * dBTdc12_x);
		qeal dvydc12_y = 2.0 * (C1y * dALdc12_y + C2y * dBTdc12_y + beta);
		qeal dvydc12_z = 2.0 * (C1y * dALdc12_z + C2y * dBTdc12_z);

		qeal dvydc13_x = 0;
		qeal dvydc13_y = 0;
		qeal dvydc13_z = 0;

		qeal dvydcp_x = 2.0 * (C1y * dALdcp_x + C2y * dBTdcp_x);
		qeal dvydcp_y = 2.0 * (C1y * dALdcp_y + C2y * dBTdcp_y - 1);
		qeal dvydcp_z = 2.0 * (C1y * dALdcp_z + C2y * dBTdcp_z);
		//
		qeal dvzdc11_x = 2.0 * (C1z * dALdc11_x + C2z * dBTdc11_x);
		qeal dvzdc11_y = 2.0 * (C1z * dALdc11_y + C2z * dBTdc11_y);
		qeal dvzdc11_z = 2.0 * (C1z * dALdc11_z + C2z * dBTdc11_z + alpha);

		qeal dvzdc12_x = 2.0 * (C1z * dALdc12_x + C2z * dBTdc12_x);
		qeal dvzdc12_y = 2.0 * (C1z * dALdc12_y + C2z * dBTdc12_y);
		qeal dvzdc12_z = 2.0 * (C1z * dALdc12_z + C2z * dBTdc12_z + beta);

		qeal dvzdc13_x = 0;
		qeal dvzdc13_y = 0;
		qeal dvzdc13_z = 0;

		qeal dvzdcp_x = 2.0 * (C1z * dALdcp_x + C2z * dBTdcp_x);
		qeal dvzdcp_y = 2.0 * (C1z * dALdcp_y + C2z * dBTdcp_y);
		qeal dvzdcp_z = 2.0 * (C1z * dALdcp_z + C2z * dBTdcp_z - 1);

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = Check_QEAL_ZERO(vx * dALdc11_x + alpha * dvxdc11_x);
		hessina.data()[1] = Check_QEAL_ZERO(vy * dALdc11_x + alpha * dvydc11_x);
		hessina.data()[2] = Check_QEAL_ZERO(vz * dALdc11_x + alpha * dvzdc11_x);

		hessina.data()[3] = Check_QEAL_ZERO(vx * dBTdc11_x + beta * dvxdc11_x);
		hessina.data()[4] = Check_QEAL_ZERO(vy * dBTdc11_x + beta * dvydc11_x);
		hessina.data()[5] = Check_QEAL_ZERO(vz * dBTdc11_x + beta * dvzdc11_x);

		hessina.data()[6] = Check_QEAL_ZERO(vx * (-dALdc11_x - dBTdc11_x) + (1.0 - alpha - beta) * dvxdc11_x);
		hessina.data()[7] = Check_QEAL_ZERO(vy * (-dALdc11_x - dBTdc11_x) + (1.0 - alpha - beta) * dvydc11_x);
		hessina.data()[8] = Check_QEAL_ZERO(vz * (-dALdc11_x - dBTdc11_x) + (1.0 - alpha - beta) * dvzdc11_x);

		hessina.data()[9] = Check_QEAL_ZERO(-dvxdc11_x);
		hessina.data()[10] = Check_QEAL_ZERO(-dvydc11_x);
		hessina.data()[11] = Check_QEAL_ZERO(-dvzdc11_x);
		// f / c11y
		hessina.data()[12] = Check_QEAL_ZERO(vx * dALdc11_y + alpha * dvxdc11_y);
		hessina.data()[13] = Check_QEAL_ZERO(vy * dALdc11_y + alpha * dvydc11_y);
		hessina.data()[14] = Check_QEAL_ZERO(vz * dALdc11_y + alpha * dvzdc11_y);

		hessina.data()[15] = Check_QEAL_ZERO(vx * dBTdc11_y + beta * dvxdc11_y);
		hessina.data()[16] = Check_QEAL_ZERO(vy * dBTdc11_y + beta * dvydc11_y);
		hessina.data()[17] = Check_QEAL_ZERO(vz * dBTdc11_y + beta * dvzdc11_y);

		hessina.data()[18] = Check_QEAL_ZERO(vx * (-dALdc11_y - dBTdc11_y) + (1.0 - alpha - beta) * dvxdc11_y);
		hessina.data()[19] = Check_QEAL_ZERO(vy * (-dALdc11_y - dBTdc11_y) + (1.0 - alpha - beta) * dvydc11_y);
		hessina.data()[20] = Check_QEAL_ZERO(vz * (-dALdc11_y - dBTdc11_y) + (1.0 - alpha - beta) * dvzdc11_y);

		hessina.data()[21] = Check_QEAL_ZERO(-dvxdc11_y);
		hessina.data()[22] = Check_QEAL_ZERO(-dvydc11_y);
		hessina.data()[23] = Check_QEAL_ZERO(-dvzdc11_y);
		// f / c11z
		hessina.data()[24] = Check_QEAL_ZERO(vx * dALdc11_z + alpha * dvxdc11_z);
		hessina.data()[25] = Check_QEAL_ZERO(vy * dALdc11_z + alpha * dvydc11_z);
		hessina.data()[26] = Check_QEAL_ZERO(vz * dALdc11_z + alpha * dvzdc11_z);

		hessina.data()[27] = Check_QEAL_ZERO(vx * dBTdc11_z + beta * dvxdc11_z);
		hessina.data()[28] = Check_QEAL_ZERO(vy * dBTdc11_z + beta * dvydc11_z);
		hessina.data()[29] = Check_QEAL_ZERO(vz * dBTdc11_z + beta * dvzdc11_z);

		hessina.data()[30] = Check_QEAL_ZERO(vx * (-dALdc11_z - dBTdc11_z) + (1.0 - alpha - beta) * dvxdc11_z);
		hessina.data()[31] = Check_QEAL_ZERO(vy * (-dALdc11_z - dBTdc11_z) + (1.0 - alpha - beta) * dvydc11_z);
		hessina.data()[32] = Check_QEAL_ZERO(vz * (-dALdc11_z - dBTdc11_z) + (1.0 - alpha - beta) * dvzdc11_z);

		hessina.data()[33] = Check_QEAL_ZERO(-dvxdc11_z);
		hessina.data()[34] = Check_QEAL_ZERO(-dvydc11_z);
		hessina.data()[35] = Check_QEAL_ZERO(-dvzdc11_z);
		// f / c12x
		hessina.data()[36] = Check_QEAL_ZERO(vx * dALdc12_x + alpha * dvxdc12_x);
		hessina.data()[37] = Check_QEAL_ZERO(vy * dALdc12_x + alpha * dvydc12_x);
		hessina.data()[38] = Check_QEAL_ZERO(vz * dALdc12_x + alpha * dvzdc12_x);

		hessina.data()[39] = Check_QEAL_ZERO(vx * dBTdc12_x + beta * dvxdc12_x);
		hessina.data()[40] = Check_QEAL_ZERO(vy * dBTdc12_x + beta * dvydc12_x);
		hessina.data()[41] = Check_QEAL_ZERO(vz * dBTdc12_x + beta * dvzdc12_x);

		hessina.data()[42] = Check_QEAL_ZERO(vx * (-dALdc12_x - dBTdc12_x) + (1.0 - alpha - beta) * dvxdc12_x);
		hessina.data()[43] = Check_QEAL_ZERO(vy * (-dALdc12_x - dBTdc12_x) + (1.0 - alpha - beta) * dvydc12_x);
		hessina.data()[44] = Check_QEAL_ZERO(vz * (-dALdc12_x - dBTdc12_x) + (1.0 - alpha - beta) * dvzdc12_x);

		hessina.data()[45] = Check_QEAL_ZERO(-dvxdc12_x);
		hessina.data()[46] = Check_QEAL_ZERO(-dvydc12_x);
		hessina.data()[47] = Check_QEAL_ZERO(-dvzdc12_x);
		// f / c12y
		hessina.data()[48] = Check_QEAL_ZERO(vx * dALdc12_y + alpha * dvxdc12_y);
		hessina.data()[49] = Check_QEAL_ZERO(vy * dALdc12_y + alpha * dvydc12_y);
		hessina.data()[50] = Check_QEAL_ZERO(vz * dALdc12_y + alpha * dvzdc12_y);

		hessina.data()[51] = Check_QEAL_ZERO(vx * dBTdc12_y + beta * dvxdc12_y);
		hessina.data()[52] = Check_QEAL_ZERO(vy * dBTdc12_y + beta * dvydc12_y);
		hessina.data()[53] = Check_QEAL_ZERO(vz * dBTdc12_y + beta * dvzdc12_y);

		hessina.data()[54] = Check_QEAL_ZERO(vx * (-dALdc12_y - dBTdc12_y) + (1.0 - alpha - beta) * dvxdc12_y);
		hessina.data()[55] = Check_QEAL_ZERO(vy * (-dALdc12_y - dBTdc12_y) + (1.0 - alpha - beta) * dvydc12_y);
		hessina.data()[56] = Check_QEAL_ZERO(vz * (-dALdc12_y - dBTdc12_y) + (1.0 - alpha - beta) * dvzdc12_y);

		hessina.data()[57] = Check_QEAL_ZERO(-dvxdc12_y);
		hessina.data()[58] = Check_QEAL_ZERO(-dvydc12_y);
		hessina.data()[59] = Check_QEAL_ZERO(-dvzdc12_y);
		// f / c12z
		hessina.data()[60] = Check_QEAL_ZERO(vx * dALdc12_z + alpha * dvxdc12_z);
		hessina.data()[61] = Check_QEAL_ZERO(vy * dALdc12_z + alpha * dvydc12_z);
		hessina.data()[62] = Check_QEAL_ZERO(vz * dALdc12_z + alpha * dvzdc12_z);

		hessina.data()[63] = Check_QEAL_ZERO(vx * dBTdc12_z + beta * dvxdc12_z);
		hessina.data()[64] = Check_QEAL_ZERO(vy * dBTdc12_z + beta * dvydc12_z);
		hessina.data()[65] = Check_QEAL_ZERO(vz * dBTdc12_z + beta * dvzdc12_z);

		hessina.data()[66] = Check_QEAL_ZERO(vx * (-dALdc12_z - dBTdc12_z) + (1.0 - alpha - beta) * dvxdc12_z);
		hessina.data()[67] = Check_QEAL_ZERO(vy * (-dALdc12_z - dBTdc12_z) + (1.0 - alpha - beta) * dvydc12_z);
		hessina.data()[68] = Check_QEAL_ZERO(vz * (-dALdc12_z - dBTdc12_z) + (1.0 - alpha - beta) * dvzdc12_z);

		hessina.data()[69] = Check_QEAL_ZERO(-dvxdc12_z);
		hessina.data()[70] = Check_QEAL_ZERO(-dvydc12_z);
		hessina.data()[71] = Check_QEAL_ZERO(-dvzdc12_z);
		//	f / c13x
		hessina.data()[72] = Check_QEAL_ZERO(vx * dALdc13_x + alpha * dvxdc13_x);
		hessina.data()[73] = Check_QEAL_ZERO(vy * dALdc13_x + alpha * dvydc13_x);
		hessina.data()[74] = Check_QEAL_ZERO(vz * dALdc13_x + alpha * dvzdc13_x);

		hessina.data()[75] = Check_QEAL_ZERO(vx * dBTdc13_x + beta * dvxdc13_x);
		hessina.data()[76] = Check_QEAL_ZERO(vy * dBTdc13_x + beta * dvydc13_x);
		hessina.data()[77] = Check_QEAL_ZERO(vz * dBTdc13_x + beta * dvzdc13_x);

		hessina.data()[78] = Check_QEAL_ZERO(vx * (-dALdc13_x - dBTdc13_x) + (1.0 - alpha - beta) * dvxdc13_x);
		hessina.data()[79] = Check_QEAL_ZERO(vy * (-dALdc13_x - dBTdc13_x) + (1.0 - alpha - beta) * dvydc13_x);
		hessina.data()[80] = Check_QEAL_ZERO(vz * (-dALdc13_x - dBTdc13_x) + (1.0 - alpha - beta) * dvzdc13_x);

		hessina.data()[81] = Check_QEAL_ZERO(-dvxdc13_x);
		hessina.data()[82] = Check_QEAL_ZERO(-dvydc13_x);
		hessina.data()[83] = Check_QEAL_ZERO(-dvzdc13_x);
		//	f / c13y
		hessina.data()[84] = Check_QEAL_ZERO(vx * dALdc13_y + alpha * dvxdc13_y);
		hessina.data()[85] = Check_QEAL_ZERO(vy * dALdc13_y + alpha * dvydc13_y);
		hessina.data()[86] = Check_QEAL_ZERO(vz * dALdc13_y + alpha * dvzdc13_y);

		hessina.data()[87] = Check_QEAL_ZERO(vx * dBTdc13_y + beta * dvxdc13_y);
		hessina.data()[88] = Check_QEAL_ZERO(vy * dBTdc13_y + beta * dvydc13_y);
		hessina.data()[89] = Check_QEAL_ZERO(vz * dBTdc13_y + beta * dvzdc13_y);

		hessina.data()[90] = Check_QEAL_ZERO(vx * (-dALdc13_y - dBTdc13_y) + (1.0 - alpha - beta) * dvxdc13_y);
		hessina.data()[91] = Check_QEAL_ZERO(vy * (-dALdc13_y - dBTdc13_y) + (1.0 - alpha - beta) * dvydc13_y);
		hessina.data()[92] = Check_QEAL_ZERO(vz * (-dALdc13_y - dBTdc13_y) + (1.0 - alpha - beta) * dvzdc13_y);

		hessina.data()[93] = Check_QEAL_ZERO(-dvxdc13_y);
		hessina.data()[94] = Check_QEAL_ZERO(-dvydc13_y);
		hessina.data()[95] = Check_QEAL_ZERO(-dvzdc13_y);
		//	f / c13z
		hessina.data()[96] = Check_QEAL_ZERO(vx * dALdc13_z + alpha * dvxdc13_z);
		hessina.data()[97] = Check_QEAL_ZERO(vy * dALdc13_z + alpha * dvydc13_z);
		hessina.data()[98] = Check_QEAL_ZERO(vz * dALdc13_z + alpha * dvzdc13_z);

		hessina.data()[99] = Check_QEAL_ZERO(vx * dBTdc13_z + beta * dvxdc13_z);
		hessina.data()[100] = Check_QEAL_ZERO(vy * dBTdc13_z + beta * dvydc13_z);
		hessina.data()[101] = Check_QEAL_ZERO(vz * dBTdc13_z + beta * dvzdc13_z);

		hessina.data()[102] = Check_QEAL_ZERO(vx * (-dALdc13_z - dBTdc13_z) + (1.0 - alpha - beta) * dvxdc13_z);
		hessina.data()[103] = Check_QEAL_ZERO(vy * (-dALdc13_z - dBTdc13_z) + (1.0 - alpha - beta) * dvydc13_z);
		hessina.data()[104] = Check_QEAL_ZERO(vz * (-dALdc13_z - dBTdc13_z) + (1.0 - alpha - beta) * dvzdc13_z);

		hessina.data()[105] = Check_QEAL_ZERO(-dvxdc13_z);
		hessina.data()[106] = Check_QEAL_ZERO(-dvydc13_z);
		hessina.data()[107] = Check_QEAL_ZERO(-dvzdc13_z);
		//	f / cpx
		hessina.data()[108] = Check_QEAL_ZERO(vx * dALdcp_x + alpha * dvxdcp_x);
		hessina.data()[109] = Check_QEAL_ZERO(vy * dALdcp_x + alpha * dvydcp_x);
		hessina.data()[110] = Check_QEAL_ZERO(vz * dALdcp_x + alpha * dvzdcp_x);

		hessina.data()[111] = Check_QEAL_ZERO(vx * dBTdcp_x + beta * dvxdcp_x);
		hessina.data()[112] = Check_QEAL_ZERO(vy * dBTdcp_x + beta * dvydcp_x);
		hessina.data()[113] = Check_QEAL_ZERO(vz * dBTdcp_x + beta * dvzdcp_x);

		hessina.data()[114] = Check_QEAL_ZERO(vx * (-dALdcp_x - dBTdcp_x) + (1.0 - alpha - beta) * dvxdcp_x);
		hessina.data()[115] = Check_QEAL_ZERO(vy * (-dALdcp_x - dBTdcp_x) + (1.0 - alpha - beta) * dvydcp_x);
		hessina.data()[116] = Check_QEAL_ZERO(vz * (-dALdcp_x - dBTdcp_x) + (1.0 - alpha - beta) * dvzdcp_x);

		hessina.data()[117] = Check_QEAL_ZERO(-dvxdcp_x);
		hessina.data()[118] = Check_QEAL_ZERO(-dvydcp_x);
		hessina.data()[119] = Check_QEAL_ZERO(-dvzdcp_x);
		//	f / cpy
		hessina.data()[120] = Check_QEAL_ZERO(vx * dALdcp_y + alpha * dvxdcp_y);
		hessina.data()[121] = Check_QEAL_ZERO(vy * dALdcp_y + alpha * dvydcp_y);
		hessina.data()[122] = Check_QEAL_ZERO(vz * dALdcp_y + alpha * dvzdcp_y);

		hessina.data()[123] = Check_QEAL_ZERO(vx * dBTdcp_y + beta * dvxdcp_y);
		hessina.data()[124] = Check_QEAL_ZERO(vy * dBTdcp_y + beta * dvydcp_y);
		hessina.data()[125] = Check_QEAL_ZERO(vz * dBTdcp_y + beta * dvzdcp_y);

		hessina.data()[126] = Check_QEAL_ZERO(vx * (-dALdcp_y - dBTdcp_y) + (1.0 - alpha - beta) * dvxdcp_y);
		hessina.data()[127] = Check_QEAL_ZERO(vy * (-dALdcp_y - dBTdcp_y) + (1.0 - alpha - beta) * dvydcp_y);
		hessina.data()[128] = Check_QEAL_ZERO(vz * (-dALdcp_y - dBTdcp_y) + (1.0 - alpha - beta) * dvzdcp_y);

		hessina.data()[129] = Check_QEAL_ZERO(-dvxdcp_y);
		hessina.data()[130] = Check_QEAL_ZERO(-dvydcp_y);
		hessina.data()[131] = Check_QEAL_ZERO(-dvzdcp_y);
		//	f / cpz
		hessina.data()[132] = Check_QEAL_ZERO(vx * dALdcp_z + alpha * dvxdcp_z);
		hessina.data()[133] = Check_QEAL_ZERO(vy * dALdcp_z + alpha * dvydcp_z);
		hessina.data()[134] = Check_QEAL_ZERO(vz * dALdcp_z + alpha * dvzdcp_z);

		hessina.data()[135] = Check_QEAL_ZERO(vx * dBTdcp_z + beta * dvxdcp_z);
		hessina.data()[136] = Check_QEAL_ZERO(vy * dBTdcp_z + beta * dvydcp_z);
		hessina.data()[137] = Check_QEAL_ZERO(vz * dBTdcp_z + beta * dvzdcp_z);

		hessina.data()[138] = Check_QEAL_ZERO(vx * (-dALdcp_z - dBTdcp_z) + (1.0 - alpha - beta) * dvxdcp_z);
		hessina.data()[139] = Check_QEAL_ZERO(vy * (-dALdcp_z - dBTdcp_z) + (1.0 - alpha - beta) * dvydcp_z);
		hessina.data()[140] = Check_QEAL_ZERO(vz * (-dALdcp_z - dBTdcp_z) + (1.0 - alpha - beta) * dvzdcp_z);

		hessina.data()[141] = Check_QEAL_ZERO(-dvxdcp_z);
		hessina.data()[142] = Check_QEAL_ZERO(-dvydcp_z);
		hessina.data()[143] = Check_QEAL_ZERO(-dvzdcp_z);
	}

	void MipcSlabSphereConstraint::alphaBetaHessina(MatrixX & hessina)
	{
		qeal dAdc11_x = 2.0 * sC1.data()[0]; qeal dAdc11_y = 2.0 * sC1.data()[1]; qeal dAdc11_z = 2.0 * sC1.data()[2];
		qeal dAdc13_x = -2.0 * sC1.data()[0]; qeal dAdc13_y = -2.0 * sC1.data()[1]; qeal dAdc13_z = -2.0 * sC1.data()[2];

		qeal dBdc11_x = 2.0 * sC2.data()[0]; qeal dBdc11_y = 2.0 * sC2.data()[1]; qeal dBdc11_z = 2.0 * sC2.data()[2];
		qeal dBdc12_x = 2.0 * sC1.data()[0]; qeal dBdc12_y = 2.0 * sC1.data()[1]; qeal dBdc12_z = 2.0 * sC1.data()[2];
		qeal dBdc13_x = -2.0 * (sC1.data()[0] + sC2.data()[0]); qeal dBdc13_y = -2.0 * (sC1.data()[1] + sC2.data()[1]); qeal dBdc13_z = -2.0 * (sC1.data()[2] + sC2.data()[2]);

		qeal dCdc12_x = 2.0 * sC2.data()[0]; qeal dCdc12_y = 2.0 * sC2.data()[1]; qeal dCdc12_z = 2.0 * sC2.data()[2];
		qeal dCdc13_x = -2.0 * sC2.data()[0]; qeal dCdc13_y = -2.0 * sC2.data()[1]; qeal dCdc13_z = -2.0 * sC2.data()[2];

		qeal dDdc11_x = 2.0 * sC3.data()[0]; qeal dDdc11_y = 2.0 * sC3.data()[1]; qeal dDdc11_z = 2.0 * sC3.data()[2];
		qeal dDdc13_x = 2.0 * (sC1.data()[0] - sC3.data()[0]); qeal dDdc13_y = 2.0 * (sC1.data()[1] - sC3.data()[1]); qeal dDdc13_z = 2.0 * (sC1.data()[2] - sC3.data()[2]);
		qeal dDdcp_x = -2.0 * sC1.data()[0]; qeal dDdcp_y = -2.0 * sC1.data()[1]; qeal dDdcp_z = -2.0 * sC1.data()[2];

		qeal dEdc12_x = 2.0 * sC3.data()[0]; qeal dEdc12_y = 2.0 * sC3.data()[1]; qeal dEdc12_z = 2.0 * sC3.data()[2];
		qeal dEdc13_x = 2.0 * (sC2.data()[0] - sC3.data()[0]); qeal dEdc13_y = 2.0 * (sC2.data()[1] - sC3.data()[1]); qeal dEdc13_z = 2.0 * (sC2.data()[2] - sC3.data()[2]);
		qeal dEdcp_x = -2.0 * sC2.data()[0]; qeal dEdcp_y = -2.0 * sC2.data()[1]; qeal dEdcp_z = -2.0 * sC2.data()[2];

		qeal dFdc13_x = 2.0 * sC3.data()[0]; qeal dFdc13_y = 2.0 * sC3.data()[1]; qeal dFdc13_z = 2.0 * sC3.data()[2];
		qeal dFdcp_x = -2.0 * sC3.data()[0]; qeal dFdcp_y = -2.0 * sC3.data()[1]; qeal dFdcp_z = -2.0 * sC3.data()[2];

		//
		qeal q1 = 1.0 / delta; qeal q2 = alpha / delta; qeal q3 = beta / delta;

		qeal dALdc11_x = q1 * (E * dBdc11_x - 2.0 * C * dDdc11_x) - q2 * (4.0 * C * dAdc11_x - 2.0 * B * dBdc11_x);
		qeal dALdc11_y = q1 * (E * dBdc11_y - 2.0 * C * dDdc11_y) - q2 * (4.0 * C * dAdc11_y - 2.0 * B * dBdc11_y);
		qeal dALdc11_z = q1 * (E * dBdc11_z - 2.0 * C * dDdc11_z) - q2 * (4.0 * C * dAdc11_z - 2.0 * B * dBdc11_z);

		qeal dALdc12_x = q1 * (E * dBdc12_x + B * dEdc12_x - 2.0 * D * dCdc12_x) - q2 * (4.0 * A * dCdc12_x - 2.0 * B * dBdc12_x);
		qeal dALdc12_y = q1 * (E * dBdc12_y + B * dEdc12_y - 2.0 * D * dCdc12_y) - q2 * (4.0 * A * dCdc12_y - 2.0 * B * dBdc12_y);
		qeal dALdc12_z = q1 * (E * dBdc12_z + B * dEdc12_z - 2.0 * D * dCdc12_z) - q2 * (4.0 * A * dCdc12_z - 2.0 * B * dBdc12_z);

		qeal dALdc13_x = q1 * (E * dBdc13_x + B * dEdc13_x - 2.0 * (C * dDdc13_x + D * dCdc13_x)) - q2 * (4.0 * (A * dCdc13_x + C * dAdc13_x) - 2.0 * B * dBdc13_x);
		qeal dALdc13_y = q1 * (E * dBdc13_y + B * dEdc13_y - 2.0 * (C * dDdc13_y + D * dCdc13_y)) - q2 * (4.0 * (A * dCdc13_y + C * dAdc13_y) - 2.0 * B * dBdc13_y);
		qeal dALdc13_z = q1 * (E * dBdc13_z + B * dEdc13_z - 2.0 * (C * dDdc13_z + D * dCdc13_z)) - q2 * (4.0 * (A * dCdc13_z + C * dAdc13_z) - 2.0 * B * dBdc13_z);

		qeal dALdcp_x = q1 * (B * dEdcp_x - 2.0 * C * dDdcp_x);
		qeal dALdcp_y = q1 * (B * dEdcp_y - 2.0 * C * dDdcp_y);
		qeal dALdcp_z = q1 * (B * dEdcp_z - 2.0 * C * dDdcp_z);

		qeal dBTdc11_x = q1 * (D * dBdc11_x + B * dDdc11_x - 2.0 * E * dAdc11_x) - q3 * (4.0 * C * dAdc11_x - 2.0 * B * dBdc11_x);
		qeal dBTdc11_y = q1 * (D * dBdc11_y + B * dDdc11_y - 2.0 * E * dAdc11_y) - q3 * (4.0 * C * dAdc11_y - 2.0 * B * dBdc11_y);
		qeal dBTdc11_z = q1 * (D * dBdc11_z + B * dDdc11_z - 2.0 * E * dAdc11_z) - q3 * (4.0 * C * dAdc11_z - 2.0 * B * dBdc11_z);

		qeal dBTdc12_x = q1 * (D * dBdc12_x - 2.0 * A * dEdc12_x) - q3 * (4.0 * A * dCdc12_x - 2.0 * B * dBdc12_x);
		qeal dBTdc12_y = q1 * (D * dBdc12_y - 2.0 * A * dEdc12_y) - q3 * (4.0 * A * dCdc12_y - 2.0 * B * dBdc12_y);
		qeal dBTdc12_z = q1 * (D * dBdc12_z - 2.0 * A * dEdc12_z) - q3 * (4.0 * A * dCdc12_z - 2.0 * B * dBdc12_z);

		qeal dBTdc13_x = q1 * (D * dBdc13_x + B * dDdc13_x - 2.0 * (A * dEdc13_x + E * dAdc13_x)) - q3 * (4.0 * (A * dCdc13_x + C * dAdc13_x) - 2.0 * B * dBdc13_x);
		qeal dBTdc13_y = q1 * (D * dBdc13_y + B * dDdc13_y - 2.0 * (A * dEdc13_y + E * dAdc13_y)) - q3 * (4.0 * (A * dCdc13_y + C * dAdc13_y) - 2.0 * B * dBdc13_y);
		qeal dBTdc13_z = q1 * (D * dBdc13_z + B * dDdc13_z - 2.0 * (A * dEdc13_z + E * dAdc13_z)) - q3 * (4.0 * (A * dCdc13_z + C * dAdc13_z) - 2.0 * B * dBdc13_z);

		qeal dBTdcp_x = q1 * (B * dDdcp_x - 2.0 * A * dEdcp_x);
		qeal dBTdcp_y = q1 * (B * dDdcp_y - 2.0 * A * dEdcp_y);
		qeal dBTdcp_z = q1 * (B * dDdcp_z - 2.0 * A * dEdcp_z);

		qeal C1x = sC1.data()[0]; qeal C1y = sC1.data()[1]; qeal C1z = sC1.data()[2];
		qeal C2x = sC2.data()[0]; qeal C2y = sC2.data()[1]; qeal C2z = sC2.data()[2];
		qeal C3x = sC3.data()[0]; qeal C3y = sC3.data()[1]; qeal C3z = sC3.data()[2];

		qeal vx = 2.0 * (alpha * C1x + beta * C2x + C3x);
		qeal vy = 2.0 * (alpha * C1y + beta * C2y + C3y);
		qeal vz = 2.0 * (alpha * C1z + beta * C2z + C3z);

		qeal dvxdc11_x = 2.0 * (C1x * dALdc11_x + C2x * dBTdc11_x + alpha);
		qeal dvxdc11_y = 2.0 * (C1x * dALdc11_y + C2x * dBTdc11_y);
		qeal dvxdc11_z = 2.0 * (C1x * dALdc11_z + C2x * dBTdc11_z);

		qeal dvxdc12_x = 2.0 * (C1x * dALdc12_x + C2x * dBTdc12_x + beta);
		qeal dvxdc12_y = 2.0 * (C1x * dALdc12_y + C2x * dBTdc12_y);
		qeal dvxdc12_z = 2.0 * (C1x * dALdc12_z + C2x * dBTdc12_z);

		qeal dvxdc13_x = 2.0 * (C1x * dALdc13_x + C2x * dBTdc13_x + 1 - alpha - beta);
		qeal dvxdc13_y = 2.0 * (C1x * dALdc13_y + C2x * dBTdc13_y);
		qeal dvxdc13_z = 2.0 * (C1x * dALdc13_z + C2x * dBTdc13_z);

		qeal dvxdcp_x = 2.0 * (C1x * dALdcp_x + C2x * dBTdcp_x - 1);
		qeal dvxdcp_y = 2.0 * (C1x * dALdcp_y + C2x * dBTdcp_y);
		qeal dvxdcp_z = 2.0 * (C1x * dALdcp_z + C2x * dBTdcp_z);
		//
		qeal dvydc11_x = 2.0 * (C1y * dALdc11_x + C2y * dBTdc11_x);
		qeal dvydc11_y = 2.0 * (C1y * dALdc11_y + C2y * dBTdc11_y + alpha);
		qeal dvydc11_z = 2.0 * (C1y * dALdc11_z + C2y * dBTdc11_z);

		qeal dvydc12_x = 2.0 * (C1y * dALdc12_x + C2y * dBTdc12_x);
		qeal dvydc12_y = 2.0 * (C1y * dALdc12_y + C2y * dBTdc12_y + beta);
		qeal dvydc12_z = 2.0 * (C1y * dALdc12_z + C2y * dBTdc12_z);

		qeal dvydc13_x = 2.0 * (C1y * dALdc13_x + C2y * dBTdc13_x);
		qeal dvydc13_y = 2.0 * (C1y * dALdc13_y + C2y * dBTdc13_y + 1 - alpha - beta);
		qeal dvydc13_z = 2.0 * (C1y * dALdc13_z + C2y * dBTdc13_z);

		qeal dvydcp_x = 2.0 * (C1y * dALdcp_x + C2y * dBTdcp_x);
		qeal dvydcp_y = 2.0 * (C1y * dALdcp_y + C2y * dBTdcp_y - 1);
		qeal dvydcp_z = 2.0 * (C1y * dALdcp_z + C2y * dBTdcp_z);
		//
		qeal dvzdc11_x = 2.0 * (C1z * dALdc11_x + C2z * dBTdc11_x);
		qeal dvzdc11_y = 2.0 * (C1z * dALdc11_y + C2z * dBTdc11_y);
		qeal dvzdc11_z = 2.0 * (C1z * dALdc11_z + C2z * dBTdc11_z + alpha);

		qeal dvzdc12_x = 2.0 * (C1z * dALdc12_x + C2z * dBTdc12_x);
		qeal dvzdc12_y = 2.0 * (C1z * dALdc12_y + C2z * dBTdc12_y);
		qeal dvzdc12_z = 2.0 * (C1z * dALdc12_z + C2z * dBTdc12_z + beta);

		qeal dvzdc13_x = 2.0 * (C1z * dALdc13_x + C2z * dBTdc13_x);
		qeal dvzdc13_y = 2.0 * (C1z * dALdc13_y + C2z * dBTdc13_y);
		qeal dvzdc13_z = 2.0 * (C1z * dALdc13_z + C2z * dBTdc13_z + 1 - alpha - beta);

		qeal dvzdcp_x = 2.0 * (C1z * dALdcp_x + C2z * dBTdcp_x);
		qeal dvzdcp_y = 2.0 * (C1z * dALdcp_y + C2z * dBTdcp_y);
		qeal dvzdcp_z = 2.0 * (C1z * dALdcp_z + C2z * dBTdcp_z - 1);

		hessina.setZero();
		// f / c11x
		hessina.data()[0] = Check_QEAL_ZERO(vx * dALdc11_x + alpha * dvxdc11_x);
		hessina.data()[1] = Check_QEAL_ZERO(vy * dALdc11_x + alpha * dvydc11_x);
		hessina.data()[2] = Check_QEAL_ZERO(vz * dALdc11_x + alpha * dvzdc11_x);

		hessina.data()[3] = Check_QEAL_ZERO(vx * dBTdc11_x + beta * dvxdc11_x);
		hessina.data()[4] = Check_QEAL_ZERO(vy * dBTdc11_x + beta * dvydc11_x);
		hessina.data()[5] = Check_QEAL_ZERO(vz * dBTdc11_x + beta * dvzdc11_x);

		hessina.data()[6] = Check_QEAL_ZERO(vx * (-dALdc11_x - dBTdc11_x) + (1.0 - alpha - beta) * dvxdc11_x);
		hessina.data()[7] = Check_QEAL_ZERO(vy * (-dALdc11_x - dBTdc11_x) + (1.0 - alpha - beta) * dvydc11_x);
		hessina.data()[8] = Check_QEAL_ZERO(vz * (-dALdc11_x - dBTdc11_x) + (1.0 - alpha - beta) * dvzdc11_x);

		hessina.data()[9] = Check_QEAL_ZERO(-dvxdc11_x);
		hessina.data()[10] = Check_QEAL_ZERO(-dvydc11_x);
		hessina.data()[11] = Check_QEAL_ZERO(-dvzdc11_x);
		// f / c11y
		hessina.data()[12] = Check_QEAL_ZERO(vx * dALdc11_y + alpha * dvxdc11_y);
		hessina.data()[13] = Check_QEAL_ZERO(vy * dALdc11_y + alpha * dvydc11_y);
		hessina.data()[14] = Check_QEAL_ZERO(vz * dALdc11_y + alpha * dvzdc11_y);

		hessina.data()[15] = Check_QEAL_ZERO(vx * dBTdc11_y + beta * dvxdc11_y);
		hessina.data()[16] = Check_QEAL_ZERO(vy * dBTdc11_y + beta * dvydc11_y);
		hessina.data()[17] = Check_QEAL_ZERO(vz * dBTdc11_y + beta * dvzdc11_y);

		hessina.data()[18] = Check_QEAL_ZERO(vx * (-dALdc11_y - dBTdc11_y) + (1.0 - alpha - beta) * dvxdc11_y);
		hessina.data()[19] = Check_QEAL_ZERO(vy * (-dALdc11_y - dBTdc11_y) + (1.0 - alpha - beta) * dvydc11_y);
		hessina.data()[20] = Check_QEAL_ZERO(vz * (-dALdc11_y - dBTdc11_y) + (1.0 - alpha - beta) * dvzdc11_y);

		hessina.data()[21] = Check_QEAL_ZERO(-dvxdc11_y);
		hessina.data()[22] = Check_QEAL_ZERO(-dvydc11_y);
		hessina.data()[23] = Check_QEAL_ZERO(-dvzdc11_y);
		// f / c11z
		hessina.data()[24] = Check_QEAL_ZERO(vx * dALdc11_z + alpha * dvxdc11_z);
		hessina.data()[25] = Check_QEAL_ZERO(vy * dALdc11_z + alpha * dvydc11_z);
		hessina.data()[26] = Check_QEAL_ZERO(vz * dALdc11_z + alpha * dvzdc11_z);

		hessina.data()[27] = Check_QEAL_ZERO(vx * dBTdc11_z + beta * dvxdc11_z);
		hessina.data()[28] = Check_QEAL_ZERO(vy * dBTdc11_z + beta * dvydc11_z);
		hessina.data()[29] = Check_QEAL_ZERO(vz * dBTdc11_z + beta * dvzdc11_z);

		hessina.data()[30] = Check_QEAL_ZERO(vx * (-dALdc11_z - dBTdc11_z) + (1.0 - alpha - beta) * dvxdc11_z);
		hessina.data()[31] = Check_QEAL_ZERO(vy * (-dALdc11_z - dBTdc11_z) + (1.0 - alpha - beta) * dvydc11_z);
		hessina.data()[32] = Check_QEAL_ZERO(vz * (-dALdc11_z - dBTdc11_z) + (1.0 - alpha - beta) * dvzdc11_z);

		hessina.data()[33] = Check_QEAL_ZERO(-dvxdc11_z);
		hessina.data()[34] = Check_QEAL_ZERO(-dvydc11_z);
		hessina.data()[35] = Check_QEAL_ZERO(-dvzdc11_z);
		// f / c12x
		hessina.data()[36] = Check_QEAL_ZERO(vx * dALdc12_x + alpha * dvxdc12_x);
		hessina.data()[37] = Check_QEAL_ZERO(vy * dALdc12_x + alpha * dvydc12_x);
		hessina.data()[38] = Check_QEAL_ZERO(vz * dALdc12_x + alpha * dvzdc12_x);

		hessina.data()[39] = Check_QEAL_ZERO(vx * dBTdc12_x + beta * dvxdc12_x);
		hessina.data()[40] = Check_QEAL_ZERO(vy * dBTdc12_x + beta * dvydc12_x);
		hessina.data()[41] = Check_QEAL_ZERO(vz * dBTdc12_x + beta * dvzdc12_x);

		hessina.data()[42] = Check_QEAL_ZERO(vx * (-dALdc12_x - dBTdc12_x) + (1.0 - alpha - beta) * dvxdc12_x);
		hessina.data()[43] = Check_QEAL_ZERO(vy * (-dALdc12_x - dBTdc12_x) + (1.0 - alpha - beta) * dvydc12_x);
		hessina.data()[44] = Check_QEAL_ZERO(vz * (-dALdc12_x - dBTdc12_x) + (1.0 - alpha - beta) * dvzdc12_x);

		hessina.data()[45] = Check_QEAL_ZERO(-dvxdc12_x);
		hessina.data()[46] = Check_QEAL_ZERO(-dvydc12_x);
		hessina.data()[47] = Check_QEAL_ZERO(-dvzdc12_x);
		// f / c12y
		hessina.data()[48] = Check_QEAL_ZERO(vx * dALdc12_y + alpha * dvxdc12_y);
		hessina.data()[49] = Check_QEAL_ZERO(vy * dALdc12_y + alpha * dvydc12_y);
		hessina.data()[50] = Check_QEAL_ZERO(vz * dALdc12_y + alpha * dvzdc12_y);

		hessina.data()[51] = Check_QEAL_ZERO(vx * dBTdc12_y + beta * dvxdc12_y);
		hessina.data()[52] = Check_QEAL_ZERO(vy * dBTdc12_y + beta * dvydc12_y);
		hessina.data()[53] = Check_QEAL_ZERO(vz * dBTdc12_y + beta * dvzdc12_y);

		hessina.data()[54] = Check_QEAL_ZERO(vx * (-dALdc12_y - dBTdc12_y) + (1.0 - alpha - beta) * dvxdc12_y);
		hessina.data()[55] = Check_QEAL_ZERO(vy * (-dALdc12_y - dBTdc12_y) + (1.0 - alpha - beta) * dvydc12_y);
		hessina.data()[56] = Check_QEAL_ZERO(vz * (-dALdc12_y - dBTdc12_y) + (1.0 - alpha - beta) * dvzdc12_y);

		hessina.data()[57] = Check_QEAL_ZERO(-dvxdc12_y);
		hessina.data()[58] = Check_QEAL_ZERO(-dvydc12_y);
		hessina.data()[59] = Check_QEAL_ZERO(-dvzdc12_y);
		// f / c12z
		hessina.data()[60] = Check_QEAL_ZERO(vx * dALdc12_z + alpha * dvxdc12_z);
		hessina.data()[61] = Check_QEAL_ZERO(vy * dALdc12_z + alpha * dvydc12_z);
		hessina.data()[62] = Check_QEAL_ZERO(vz * dALdc12_z + alpha * dvzdc12_z);

		hessina.data()[63] = Check_QEAL_ZERO(vx * dBTdc12_z + beta * dvxdc12_z);
		hessina.data()[64] = Check_QEAL_ZERO(vy * dBTdc12_z + beta * dvydc12_z);
		hessina.data()[65] = Check_QEAL_ZERO(vz * dBTdc12_z + beta * dvzdc12_z);

		hessina.data()[66] = Check_QEAL_ZERO(vx * (-dALdc12_z - dBTdc12_z) + (1.0 - alpha - beta) * dvxdc12_z);
		hessina.data()[67] = Check_QEAL_ZERO(vy * (-dALdc12_z - dBTdc12_z) + (1.0 - alpha - beta) * dvydc12_z);
		hessina.data()[68] = Check_QEAL_ZERO(vz * (-dALdc12_z - dBTdc12_z) + (1.0 - alpha - beta) * dvzdc12_z);

		hessina.data()[69] = Check_QEAL_ZERO(-dvxdc12_z);
		hessina.data()[70] = Check_QEAL_ZERO(-dvydc12_z);
		hessina.data()[71] = Check_QEAL_ZERO(-dvzdc12_z);
		//	f / c13x
		hessina.data()[72] = Check_QEAL_ZERO(vx * dALdc13_x + alpha * dvxdc13_x);
		hessina.data()[73] = Check_QEAL_ZERO(vy * dALdc13_x + alpha * dvydc13_x);
		hessina.data()[74] = Check_QEAL_ZERO(vz * dALdc13_x + alpha * dvzdc13_x);

		hessina.data()[75] = Check_QEAL_ZERO(vx * dBTdc13_x + beta * dvxdc13_x);
		hessina.data()[76] = Check_QEAL_ZERO(vy * dBTdc13_x + beta * dvydc13_x);
		hessina.data()[77] = Check_QEAL_ZERO(vz * dBTdc13_x + beta * dvzdc13_x);

		hessina.data()[78] = Check_QEAL_ZERO(vx * (-dALdc13_x - dBTdc13_x) + (1.0 - alpha - beta) * dvxdc13_x);
		hessina.data()[79] = Check_QEAL_ZERO(vy * (-dALdc13_x - dBTdc13_x) + (1.0 - alpha - beta) * dvydc13_x);
		hessina.data()[80] = Check_QEAL_ZERO(vz * (-dALdc13_x - dBTdc13_x) + (1.0 - alpha - beta) * dvzdc13_x);

		hessina.data()[81] = Check_QEAL_ZERO(-dvxdc13_x);
		hessina.data()[82] = Check_QEAL_ZERO(-dvydc13_x);
		hessina.data()[83] = Check_QEAL_ZERO(-dvzdc13_x);
		//	f / c13y
		hessina.data()[84] = Check_QEAL_ZERO(vx * dALdc13_y + alpha * dvxdc13_y);
		hessina.data()[85] = Check_QEAL_ZERO(vy * dALdc13_y + alpha * dvydc13_y);
		hessina.data()[86] = Check_QEAL_ZERO(vz * dALdc13_y + alpha * dvzdc13_y);

		hessina.data()[87] = Check_QEAL_ZERO(vx * dBTdc13_y + beta * dvxdc13_y);
		hessina.data()[88] = Check_QEAL_ZERO(vy * dBTdc13_y + beta * dvydc13_y);
		hessina.data()[89] = Check_QEAL_ZERO(vz * dBTdc13_y + beta * dvzdc13_y);

		hessina.data()[90] = Check_QEAL_ZERO(vx * (-dALdc13_y - dBTdc13_y) + (1.0 - alpha - beta) * dvxdc13_y);
		hessina.data()[91] = Check_QEAL_ZERO(vy * (-dALdc13_y - dBTdc13_y) + (1.0 - alpha - beta) * dvydc13_y);
		hessina.data()[92] = Check_QEAL_ZERO(vz * (-dALdc13_y - dBTdc13_y) + (1.0 - alpha - beta) * dvzdc13_y);

		hessina.data()[93] = Check_QEAL_ZERO(-dvxdc13_y);
		hessina.data()[94] = Check_QEAL_ZERO(-dvydc13_y);
		hessina.data()[95] = Check_QEAL_ZERO(-dvzdc13_y);
		//	f / c13z
		hessina.data()[96] = Check_QEAL_ZERO(vx * dALdc13_z + alpha * dvxdc13_z);
		hessina.data()[97] = Check_QEAL_ZERO(vy * dALdc13_z + alpha * dvydc13_z);
		hessina.data()[98] = Check_QEAL_ZERO(vz * dALdc13_z + alpha * dvzdc13_z);

		hessina.data()[99] = Check_QEAL_ZERO(vx * dBTdc13_z + beta * dvxdc13_z);
		hessina.data()[100] = Check_QEAL_ZERO(vy * dBTdc13_z + beta * dvydc13_z);
		hessina.data()[101] = Check_QEAL_ZERO(vz * dBTdc13_z + beta * dvzdc13_z);

		hessina.data()[102] = Check_QEAL_ZERO(vx * (-dALdc13_z - dBTdc13_z) + (1.0 - alpha - beta) * dvxdc13_z);
		hessina.data()[103] = Check_QEAL_ZERO(vy * (-dALdc13_z - dBTdc13_z) + (1.0 - alpha - beta) * dvydc13_z);
		hessina.data()[104] = Check_QEAL_ZERO(vz * (-dALdc13_z - dBTdc13_z) + (1.0 - alpha - beta) * dvzdc13_z);

		hessina.data()[105] = Check_QEAL_ZERO(-dvxdc13_z);
		hessina.data()[106] = Check_QEAL_ZERO(-dvydc13_z);
		hessina.data()[107] = Check_QEAL_ZERO(-dvzdc13_z);
		//	f / cpx
		hessina.data()[108] = Check_QEAL_ZERO(vx * dALdcp_x + alpha * dvxdcp_x);
		hessina.data()[109] = Check_QEAL_ZERO(vy * dALdcp_x + alpha * dvydcp_x);
		hessina.data()[110] = Check_QEAL_ZERO(vz * dALdcp_x + alpha * dvzdcp_x);

		hessina.data()[111] = Check_QEAL_ZERO(vx * dBTdcp_x + beta * dvxdcp_x);
		hessina.data()[112] = Check_QEAL_ZERO(vy * dBTdcp_x + beta * dvydcp_x);
		hessina.data()[113] = Check_QEAL_ZERO(vz * dBTdcp_x + beta * dvzdcp_x);

		hessina.data()[114] = Check_QEAL_ZERO(vx * (-dALdcp_x - dBTdcp_x) + (1.0 - alpha - beta) * dvxdcp_x);
		hessina.data()[115] = Check_QEAL_ZERO(vy * (-dALdcp_x - dBTdcp_x) + (1.0 - alpha - beta) * dvydcp_x);
		hessina.data()[116] = Check_QEAL_ZERO(vz * (-dALdcp_x - dBTdcp_x) + (1.0 - alpha - beta) * dvzdcp_x);

		hessina.data()[117] = Check_QEAL_ZERO(-dvxdcp_x);
		hessina.data()[118] = Check_QEAL_ZERO(-dvydcp_x);
		hessina.data()[119] = Check_QEAL_ZERO(-dvzdcp_x);
		//	f / cpy
		hessina.data()[120] = Check_QEAL_ZERO(vx * dALdcp_y + alpha * dvxdcp_y);
		hessina.data()[121] = Check_QEAL_ZERO(vy * dALdcp_y + alpha * dvydcp_y);
		hessina.data()[122] = Check_QEAL_ZERO(vz * dALdcp_y + alpha * dvzdcp_y);

		hessina.data()[123] = Check_QEAL_ZERO(vx * dBTdcp_y + beta * dvxdcp_y);
		hessina.data()[124] = Check_QEAL_ZERO(vy * dBTdcp_y + beta * dvydcp_y);
		hessina.data()[125] = Check_QEAL_ZERO(vz * dBTdcp_y + beta * dvzdcp_y);

		hessina.data()[126] = Check_QEAL_ZERO(vx * (-dALdcp_y - dBTdcp_y) + (1.0 - alpha - beta) * dvxdcp_y);
		hessina.data()[127] = Check_QEAL_ZERO(vy * (-dALdcp_y - dBTdcp_y) + (1.0 - alpha - beta) * dvydcp_y);
		hessina.data()[128] = Check_QEAL_ZERO(vz * (-dALdcp_y - dBTdcp_y) + (1.0 - alpha - beta) * dvzdcp_y);

		hessina.data()[129] = Check_QEAL_ZERO(-dvxdcp_y);
		hessina.data()[130] = Check_QEAL_ZERO(-dvydcp_y);
		hessina.data()[131] = Check_QEAL_ZERO(-dvzdcp_y);
		//	f / cpz
		hessina.data()[132] = Check_QEAL_ZERO(vx * dALdcp_z + alpha * dvxdcp_z);
		hessina.data()[133] = Check_QEAL_ZERO(vy * dALdcp_z + alpha * dvydcp_z);
		hessina.data()[134] = Check_QEAL_ZERO(vz * dALdcp_z + alpha * dvzdcp_z);

		hessina.data()[135] = Check_QEAL_ZERO(vx * dBTdcp_z + beta * dvxdcp_z);
		hessina.data()[136] = Check_QEAL_ZERO(vy * dBTdcp_z + beta * dvydcp_z);
		hessina.data()[137] = Check_QEAL_ZERO(vz * dBTdcp_z + beta * dvzdcp_z);

		hessina.data()[138] = Check_QEAL_ZERO(vx * (-dALdcp_z - dBTdcp_z) + (1.0 - alpha - beta) * dvxdcp_z);
		hessina.data()[139] = Check_QEAL_ZERO(vy * (-dALdcp_z - dBTdcp_z) + (1.0 - alpha - beta) * dvydcp_z);
		hessina.data()[140] = Check_QEAL_ZERO(vz * (-dALdcp_z - dBTdcp_z) + (1.0 - alpha - beta) * dvzdcp_z);

		hessina.data()[141] = Check_QEAL_ZERO(-dvxdcp_z);
		hessina.data()[142] = Check_QEAL_ZERO(-dvydcp_z);
		hessina.data()[143] = Check_QEAL_ZERO(-dvzdcp_z);
	}

}