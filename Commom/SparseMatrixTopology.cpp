#include "SparseMatrixTopology.h"
#include <iostream>

SparseMatrixTopology::SparseMatrixTopology(SparseMatrix* sparseMatrix)
{
	rows_ = sparseMatrix->rows();
	cols_ = sparseMatrix->cols();
	//copy sparseMatrix topology and set all value to 0 first
	this->sparse_matrix.resize(sparseMatrix->rows(), sparseMatrix->cols());
	std::vector<TripletX> resultcoef;
	int i = 0;
	for (int j = 0; j < sparseMatrix->outerSize(); ++j)
		for (SparseMatrix::InnerIterator it(*sparseMatrix, j); it; ++it, i++)
		{
			resultcoef.push_back(TripletX(it.row(), it.col(), 0));
			mapOfSparse.insert(std::make_pair(it.col()*sparseMatrix->rows() + it.row(), i));
		}
	sparse_matrix.setFromTriplets(resultcoef.begin(), resultcoef.end());
}


SparseMatrixTopology::~SparseMatrixTopology()
{

}

int SparseMatrixTopology::getValueIndex(int row, int col)
{

	return mapOfSparse[col*rows_ + row];

	int i = 0;

	for (int j = 0; j < sparse_matrix.outerSize(); j++)
	{
		for (SparseMatrix::InnerIterator it(sparse_matrix, j); it; ++it, i++)
		{
			if (row == it.row() && col == it.col())
			{
				return i;
			}
		}
	}

	return -1;
}