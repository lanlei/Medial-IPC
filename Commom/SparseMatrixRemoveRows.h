#pragma once
#ifndef SPARSE_MATRIX_REMOVE_ROWS_H
#define SPARSE_MATRIX_REMOVE_ROWS_H
#include "SparseMatrixTopology.h"
#include <iostream>

void subSparseMatrix(SparseMatrix &source, SparseMatrix &result, std::vector<int> &map);

void SparseMatrixRemoveRows(SparseMatrix* sparseMatrix, SparseMatrix* resultMatrix, std::vector<int>& entry_map, int r, int numConstrainedDOFs_, int* constrainedDOFs_);

void VectorRemoveRows(std::vector<int> &map, VectorX &target, VectorX &result, int numConstrainedDOFs_, int* constrainedDOFs_);

void VectorInsertRows(std::vector<int> &map, VectorX &target, VectorX &result, int numConstrainedDOFs_, int* constrainedDOFs_);

void MatrixRemoveDofs(std::vector<int> &map, MatrixX& target, MatrixX &result);

void MatrixInsertDofs(std::vector<int> &map, MatrixX& target, MatrixX &result);

void MatrixRemoveCols(std::vector<int> &map, MatrixX& target, MatrixX &result);

void MatrixRemoveRows(std::vector<int> &map, MatrixX& target, MatrixX &result);

void createMapByConstrains(std::vector<int> &map, int r, int numConstrainedDOFs_, int* constrainedDOFs_);

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
void createSparseMapbyTopology(SparseMatrix* sparseMatrix, SparseMatrix* subsparseMatrix, std::vector<int>& entry_map, std::vector<int>& rowmap, int r, int numConstrainedDOFs_, int* constrainedDOFs_);

#endif