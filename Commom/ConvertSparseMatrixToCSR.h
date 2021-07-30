#pragma once
#ifndef CONVERT_SPARSE_MATRIX_TO_CSR_H
#define CONVERT_SPARSE_MATRIX_TO_CSR_H
#include "MatrixCore.h"

static void convertSparseMatrixToCSR(SparseMatrix& sm, std::vector<int>& coeff_csrRowPtr, std::vector<int>& coeff_csrColInd, std::vector<qeal>& coeff_csrVal, int& coeff_non_zero_num)
{
	int rows = sm.rows();
	int cols = sm.cols();
	coeff_non_zero_num = sm.nonZeros();

	coeff_csrRowPtr.resize(rows + 1);
	coeff_csrColInd.resize(coeff_non_zero_num);
	coeff_csrVal.resize(coeff_non_zero_num);

	std::vector<std::vector<int>> row_nnz_index;
	row_nnz_index.resize(rows);
	std::vector<std::vector <qeal>> spMat;
	spMat.resize(rows);

	for (int k = 0; k < sm.outerSize(); ++k)
	{
		for (SparseMatrix::InnerIterator it(sm, k); it; ++it)
		{
			qeal value = it.value();
			int row_id = it.row();
			int col_id = it.col();

			row_nnz_index[row_id].push_back(col_id);
			spMat[row_id].push_back(value);
		}
	}

	int count = 0;
	for (int i = 0; i < rows; i++)
	{
		int size = row_nnz_index[i].size();
		for (int j = 0; j < size; j++)
		{
			coeff_csrVal[count] = spMat[i][j];
			coeff_csrColInd[count] = row_nnz_index[i][j];

			int row_id = i;
			int col_id = row_nnz_index[i][j];
			count++;
		}
	}

	int perv_size = 0;
	for (int i = 0; i < rows; i++)
	{
		if (i == 0)
		{
			coeff_csrRowPtr[0] = 0;
		}
		else
		{
			coeff_csrRowPtr[i] = perv_size;
		}
		perv_size += row_nnz_index[i].size();
	}

	coeff_csrRowPtr[rows] = coeff_non_zero_num;
}

static void convertSparseMatrixToCSR(SparseMatrix& sm, std::vector<int>& coeff_csrRowPtr, std::vector<int>& coeff_csrColInd, std::vector<qeal>& coeff_csrVal, int& coeff_non_zero_num, std::vector<int>& smCsrIndexMap)
{
	int rows = sm.rows();
	int cols = sm.cols();
	coeff_non_zero_num = sm.nonZeros();

	coeff_csrRowPtr.resize(rows + 1);
	coeff_csrColInd.resize(coeff_non_zero_num);
	coeff_csrVal.resize(coeff_non_zero_num);

	std::vector<std::vector<int>> row_nnz_index;
	row_nnz_index.resize(rows);
	std::vector<std::vector <qeal>> spMat;
	spMat.resize(rows);

	for (int k = 0; k < sm.outerSize(); ++k)
	{
		for (SparseMatrix::InnerIterator it(sm, k); it; ++it)
		{
			qeal value = it.value();
			int row_id = it.row();
			int col_id = it.col();

			row_nnz_index[row_id].push_back(col_id);
			spMat[row_id].push_back(value);
		}
	}

	int count = 0;
	for (int i = 0; i < rows; i++)
	{
		int size = row_nnz_index[i].size();
		for (int j = 0; j < size; j++)
		{
			coeff_csrVal[count] = spMat[i][j];
			coeff_csrColInd[count] = row_nnz_index[i][j];

			int row_id = i;
			int col_id = row_nnz_index[i][j];

			smCsrIndexMap.push_back(col_id * rows + row_id);
			count++;
		}
	}

	int perv_size = 0;
	for (int i = 0; i < rows; i++)
	{
		if (i == 0)
		{
			coeff_csrRowPtr[0] = 0;
		}
		else
		{
			coeff_csrRowPtr[i] = perv_size;
		}
		perv_size += row_nnz_index[i].size();
	}

	coeff_csrRowPtr[rows] = coeff_non_zero_num;
}

static void getSpareMatrixNeighborInfo(SparseMatrix& sm, std::vector<std::vector<int>>& neighborIndices, std::vector<qeal>& mainCoeff, std::vector<std::vector<qeal>>& neighborCoeff)
{
	int dim = sm.rows();
	neighborIndices.resize(dim);
	mainCoeff.resize(dim);
	neighborCoeff.resize(dim);

	for (int k = 0; k < sm.outerSize(); ++k)
	{
		for (SparseMatrix::InnerIterator it(sm, k); it; ++it)
		{
			qeal value = it.value();
			int row_id = it.row();
			int col_id = it.col();

			if (row_id == col_id)
			{
				mainCoeff[row_id] = value;
				continue;
			}
			else
			{
				neighborIndices[row_id].push_back(col_id);
				neighborCoeff[row_id].push_back(value);
			}
		}
	}
}

static void convertSparseMatrixFromCSR(int dim, SparseMatrix& sm, std::vector<int>& coeff_csrRowPtr, std::vector<int>& coeff_csrColInd, std::vector<qeal>& coeff_csrVal, int& coeff_non_zero_num)
{
	sm.resize(dim, dim);
	
	std::vector<TripletX> triplet;
	for (int i = 0; i < dim; i++) // rows
	{
		
		int size = 0;
		int buffer = coeff_csrRowPtr[i];
		size = coeff_csrRowPtr[i + 1] - coeff_csrRowPtr[i];

		for (int j = 0; j < size; j++)
		{			
			int col_id = coeff_csrColInd[buffer + j];
			qeal val = coeff_csrVal[buffer + j];
			triplet.push_back(TripletX(i, col_id, val));
		}
	}
	sm.setFromTriplets(triplet.begin(), triplet.end());
}

#endif