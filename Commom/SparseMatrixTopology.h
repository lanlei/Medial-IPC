#pragma once
#ifndef SPARSE_MATRIX_TOPOLOGY_H
#define SPARSE_MATRIX_TOPOLOGY_H
#include <map>
#include "MatrixCore.h"

class SparseMatrixTopology
{
public:
	SparseMatrixTopology(SparseMatrix* sparseMatrix);
	~SparseMatrixTopology();

	int getValueIndex(int row, int col);
protected:

	SparseMatrix sparse_matrix;
	std::map<int, int> mapOfSparse;
	std::map<std::vector<int>, int> mapOfSparsev2;
	int rows_;
	int cols_;

};

template <class TYPE>
class SparseMatrixTopologyTYPE
{
	typedef Eigen::Triplet<TYPE> EIGEN_TRI;
public:
	SparseMatrixTopologyTYPE(Eigen::SparseMatrix<TYPE>* sparseMatrix);
	~SparseMatrixTopologyTYPE();

	int getValueIndex(int row, int col);

protected:

	Eigen::SparseMatrix<TYPE> sparse_matrix;
	std::map<int, int> mapOfSparse;
	std::map<std::vector<int>, int> mapOfSparsev2;
	int rows_;
	int cols_;
};

template <class TYPE>
int SparseMatrixTopologyTYPE<TYPE>::getValueIndex(int row, int col)
{
	std::vector<int> rowcolindex(2);
	rowcolindex[0] = row;
	rowcolindex[1] = col;

	return mapOfSparsev2[rowcolindex];
}

template <class TYPE>
SparseMatrixTopologyTYPE<TYPE>::~SparseMatrixTopologyTYPE()
{

}

template <class TYPE>
SparseMatrixTopologyTYPE<TYPE>::SparseMatrixTopologyTYPE(Eigen::SparseMatrix<TYPE>* sparseMatrix)
{
	rows_ = sparseMatrix->rows();
	cols_ = sparseMatrix->cols();

	this->sparse_matrix.resize(sparseMatrix->rows(), sparseMatrix->cols());
	std::vector<EIGEN_TRI> resultcoef;

	int i = 0;
	for (int j = 0; j < sparseMatrix->outerSize(); ++j)
		for (typename Eigen::SparseMatrix<TYPE>::InnerIterator it(*sparseMatrix, j); it; ++it)
		{
			resultcoef.push_back(EIGEN_TRI(it.row(), it.col(), 0));
			std::vector<int> rowcolindex(2);
			rowcolindex[0] = it.row();
			rowcolindex[1] = it.col();
			mapOfSparsev2.insert(std::make_pair(rowcolindex, i));
			i++;
		}
	sparse_matrix.setFromTriplets(resultcoef.begin(), resultcoef.end());

}


#endif