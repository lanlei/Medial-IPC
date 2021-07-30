#include "SparseMatrixRemoveRows.h"

void subSparseMatrix(SparseMatrix &source, SparseMatrix &result, std::vector<int> &map)
{
	std::vector<TripletX> reusltcoef;
	if (map.size() == source.rows())
	{
		for (int j = 0; j < source.outerSize(); ++j)
			for (SparseMatrix::InnerIterator it(source, j); it; ++it)
			{
				int i_p = map[it.row()];
				int j_p = map[it.col()];
				if (!(i_p == -1 || j_p == -1))
				{
					reusltcoef.push_back(TripletX(i_p, j_p, it.value()));
				}
			}
		result.setFromTriplets(reusltcoef.begin(), reusltcoef.end());
	}
	else
		if (map.size() == result.rows())
		{
			std::cout << "sparse matrix does not support this map" << std::endl;
			return;
		}
}

void SparseMatrixRemoveRows(SparseMatrix* sparseMatrix, SparseMatrix* resultMatrix, std::vector<int>& entry_map, int r, int numConstrainedDOFs_, int* constrainedDOFs_)
{
	int i = 0;
	for (int j = 0; j < sparseMatrix->outerSize(); ++j)
		for (SparseMatrix::InnerIterator it(*sparseMatrix, j); it; ++it, ++i)
		{
			if (entry_map[i] != -1)
			{
				resultMatrix->valuePtr()[entry_map[i]] = it.value();
			}
		}
}

void VectorRemoveRows(std::vector<int> &map, VectorX &target, VectorX &result, int numConstrainedDOFs_, int* constrainedDOFs_)
{
	for (int i = 0; i < target.size(); i++)
	{
		if (map[i] != -1)
		{
			result.data()[map[i]] = target.data()[i];
		}
	}
}

void VectorInsertRows(std::vector<int> &map, VectorX &target, VectorX &result, int numConstrainedDOFs_, int* constrainedDOFs_)
{
	result.setZero();
	for (int i = 0; i < result.size(); i++)
	{
		if (map[i] != -1)
		{
			result.data()[i] = target.data()[map[i]];
		}
	}
}

void MatrixRemoveDofs(std::vector<int> &map, MatrixX& target, MatrixX &result)
{
	result.setZero();
	for (int i = 0; i < target.rows(); i++)
	{
		for (int j = 0; j < target.cols(); j++)
		{
			if (map[i] != -1 && map[j] != -1)
			{
				int row = map[i];
				int col = map[j];
				result.data()[col*result.rows() + row] = target.data()[j*target.rows() + i];
			}
		}
	}
}

void MatrixInsertDofs(std::vector<int> &map, MatrixX& target, MatrixX &result)
{
	result.setZero();
	for (int i = 0; i < result.rows(); i++)
	{
		for (int j = 0; j < result.cols(); j++)
		{
			if (map[i] != -1 && map[j] != -1)
			{
				int row = map[i];
				int col = map[j];
				result.data()[j*result.rows() + i] = target.data()[col*target.rows() + row];
			}
		}
	}
}

void MatrixRemoveCols(std::vector<int> &map, MatrixX& target, MatrixX &result)
{
	result.setZero();
	for (int i = 0; i < result.rows(); i++)
	{
		for (int j = 0; j < result.cols(); j++)
		{
			if (map[j] != -1)
			{
				int row = i;
				int col = map[j];
				result.data()[j*result.rows() + i] = target.data()[col*target.rows() + row];
			}
		}
	}
}


void MatrixRemoveRows(std::vector<int> &map, MatrixX& target, MatrixX &result)
{
	result.setZero();
	for (int i = 0; i < result.rows(); i++)
	{
		for (int j = 0; j < result.cols(); j++)
		{
			if (map[i] != -1)
			{
				int row = map[i];
				int col = j;
				result.data()[j*result.rows() + i] = target.data()[col*target.rows() + row];
			}
		}
	}
}

void createMapByConstrains(std::vector<int> &map, int r, int numConstrainedDOFs_, int* constrainedDOFs_)
{
	map.resize(r);
	std::fill(map.begin(), map.end(), 0);
	for (int i = 0; i < numConstrainedDOFs_; i++)
	{
		map[constrainedDOFs_[i]] = -1;
	}

	int ofset = 0;
	for (int i = 0; i < r; i++)
	{
		if (map[i] != -1)
		{
			map[i] = ofset;
			ofset++;
		}
	}
}

/// <summary>
/// Creates the sparse map by topology.
/// </summary>
/// <param name="sparseMatrix">The sparse matrix.</param>
/// <param name="subsparseMatrix">The subsparse matrix.</param>
/// <param name="entry_map">The entry_map.</param>
/// <param name="rowmap">The rowmap.</param>
/// <param name="r">The r.</param>
/// <param name="numConstrainedDOFs_">The number constrained do FS_.</param>
/// <param name="constrainedDOFs_">The constrained do FS_.</param>
void createSparseMapbyTopology(SparseMatrix* sparseMatrix, SparseMatrix* subsparseMatrix, std::vector<int>& entry_map, std::vector<int>& rowmap, int r, int numConstrainedDOFs_, int* constrainedDOFs_)
{
	subsparseMatrix->resize(r - numConstrainedDOFs_, r - numConstrainedDOFs_);
	std::vector<TripletX> subcoef;

	for (int j = 0; j < sparseMatrix->outerSize(); ++j)
		for (SparseMatrix::InnerIterator it(*sparseMatrix, j); it; ++it)
		{
			int i_p = rowmap[it.row()];
			int j_p = rowmap[it.col()];
			if (!(i_p == -1 || j_p == -1))
			{
				subcoef.push_back(TripletX(i_p, j_p, it.value()));
			}
		}
	subsparseMatrix->setFromTriplets(subcoef.begin(), subcoef.end());

	int supersize = sparseMatrix->nonZeros();
	int subsize = subsparseMatrix->nonZeros();

	entry_map.resize(supersize);
	std::fill(entry_map.begin(), entry_map.end(), -1);

	SparseMatrixTopology supermatrix(sparseMatrix);
	SparseMatrixTopology submatrix(subsparseMatrix);

	for (int j = 0; j < sparseMatrix->outerSize(); ++j)
		for (SparseMatrix::InnerIterator it(*sparseMatrix, j); it; ++it)
		{
			int i_p = rowmap[it.row()];
			int j_p = rowmap[it.col()];
			if (!(i_p == -1 || j_p == -1))
			{
				entry_map[supermatrix.getValueIndex(it.row(), it.col())] =
					submatrix.getValueIndex(i_p, j_p);
			}
		}
}