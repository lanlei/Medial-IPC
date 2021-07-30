#pragma once
#ifndef MATRIX_CORE_H
#define MATRIX_CORE_H
#ifndef  EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#endif // ! EIGEN_USE_MKL_ALL

#include <Eigen/Eigen>
#include <Eigen/PardisoSupport>
#include "DataCore.h"

typedef Eigen::Matrix<qeal, 2, 1> Vector2;
typedef Eigen::Matrix<qeal, 1, 2> VectorR2;
typedef Eigen::Matrix<qeal, 3, 1> Vector3;
typedef Eigen::Matrix<qeal, 1, 3> VectorR3;
typedef Eigen::Matrix<qeal, 4, 1> Vector4;
typedef Eigen::Matrix<qeal, 1, 4> VectorR4;
typedef Eigen::Matrix<qeal, Eigen::Dynamic, 1> VectorX;
typedef Eigen::Matrix<qeal, 1, Eigen::Dynamic> VectorXR;

typedef Eigen::Matrix<qeal, 2, 2> Matrix2;
typedef Eigen::Matrix<qeal, 3, 3> Matrix3;
typedef Eigen::Matrix<qeal, 4, 4> Matrix4;
typedef Eigen::Matrix<qeal, 2, Eigen::Dynamic> Matrix2X;
typedef Eigen::Matrix<qeal, Eigen::Dynamic, 2> MatrixXR2;
typedef Eigen::Matrix<qeal, 3, Eigen::Dynamic> Matrix3X;
typedef Eigen::Matrix<qeal, Eigen::Dynamic, 3> MatrixXR3;
typedef Eigen::Matrix<qeal, 4, Eigen::Dynamic> Matrix4X;
typedef Eigen::Matrix<qeal, Eigen::Dynamic, 4> MatrixXR4;
typedef Eigen::Matrix<qeal, Eigen::Dynamic, Eigen::Dynamic> MatrixX;

typedef Eigen::Matrix<int, 2, 1> Vector2i;
typedef Eigen::Matrix<int, 1, 2> VectorR2i;
typedef Eigen::Matrix<int, 3, 1> Vector3i;
typedef Eigen::Matrix<int, 1, 3> VectorR3i;
typedef Eigen::Matrix<int, 4, 1> Vector4i;
typedef Eigen::Matrix<int, 1, 4> VectorR4i;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorXi;
typedef Eigen::Matrix<int, 1, Eigen::Dynamic> VectorXRi;

typedef Eigen::Matrix<int, 2, 2> Matrix2i;
typedef Eigen::Matrix<int, 3, 3> Matrix3i;
typedef Eigen::Matrix<int, 4, 4> Matrix4i;
typedef Eigen::Matrix<int, 2, Eigen::Dynamic> Matrix2Xi;
typedef Eigen::Matrix<int, Eigen::Dynamic, 2> MatrixXR2i;
typedef Eigen::Matrix<int, 3, Eigen::Dynamic> Matrix3Xi;
typedef Eigen::Matrix<int, Eigen::Dynamic, 3> MatrixXR3i;
typedef Eigen::Matrix<int, 4, Eigen::Dynamic> Matrix4Xi;
typedef Eigen::Matrix<int, Eigen::Dynamic, 4> MatrixXR4i;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> MatrixXi;

typedef Eigen::SparseMatrix<qeal/*, Eigen::RowMajor*/> SparseMatrix;
typedef Eigen::Triplet<qeal> TripletX;
typedef Eigen::Triplet<int> TripletXi;

#endif