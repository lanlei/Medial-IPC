#pragma once
#pragma once
#define outputprecision 64
#include <fstream>
#include <iostream>
#include "MatrixCore.h"

namespace EigenMatrixIO {
	using namespace Eigen;
	template<class Matrix>
	void write_binary(const char* filename, const Matrix& matrix) {
		std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
		out.precision(outputprecision);
		typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
		out.write((char*)(&rows), sizeof(typename Matrix::Index));
		out.write((char*)(&cols), sizeof(typename Matrix::Index));
		out.write((char*)matrix.data(), rows*cols * sizeof(typename Matrix::Scalar));
		out.close();
	};

	template<class Matrix>
	void write_binary(std::ofstream & out, const Matrix& matrix) {
		typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
		out.write((char*)(&rows), sizeof(typename Matrix::Index));
		out.write((char*)(&cols), sizeof(typename Matrix::Index));
		out.write((char*)matrix.data(), rows*cols * sizeof(typename Matrix::Scalar));
	};

	template<class Matrix>
	bool read_binary(const char* filename, Matrix& matrix) {
		std::ifstream in(filename, std::ios::in | std::ios::binary);
		if (!in.good())
		{
			std::cout << "file not open" << std::endl;
			return false;
		}
		typename Matrix::Index rows = 0, cols = 0;
		in.read((char*)(&rows), sizeof(typename Matrix::Index));
		in.read((char*)(&cols), sizeof(typename Matrix::Index));
		matrix.resize(rows, cols);
		in.read((char *)matrix.data(), rows*cols * sizeof(typename Matrix::Scalar));
		in.close();
		return true;
	};

	template<class Matrix>
	void read_binary(std::ifstream &in, Matrix& matrix) {
		/*if (!in.good())
		{
		std::cout << "file not open" << std::endl;
		return;
		}*/
		typename Matrix::Index rows = 0, cols = 0;
		in.read((char*)(&rows), sizeof(typename Matrix::Index));
		in.read((char*)(&cols), sizeof(typename Matrix::Index));
		matrix.resize(rows, cols);
		in.read((char *)matrix.data(), rows*cols * sizeof(typename Matrix::Scalar));
	};

	template<typename T>
	void write_sp_binary(std::ofstream &out, Eigen::SparseMatrix<T>& sp)
	{
		int rows = sp.rows();
		int cols = sp.cols();
		out.write((char*)&rows, sizeof(int));
		out.write((char*)&cols, sizeof(int));
		int nnz = sp.nonZeros();
		out.write((char*)&nnz, sizeof(int));

		for (int k = 0; k < sp.outerSize(); ++k)
			for (Eigen::SparseMatrix<T>::InnerIterator it(sp, k); it; ++it)
			{
				int row = it.row();
				int col = it.col();
				T value = it.value();
				out.write((char*)&row, sizeof(int));
				out.write((char*)&col, sizeof(int));
				out.write((char*)&value, sizeof(T));
			}
	}

	template<typename T>
	void read_sp_binary(std::ifstream &in, Eigen::SparseMatrix<T>& sp)
	{
		int rows;
		int cols;
		in.read((char*)&rows, sizeof(int));
		in.read((char*)&cols, sizeof(int));
		sp.resize(rows, cols);
		int nnz;
		in.read((char*)&nnz, sizeof(int));
		std::vector<Eigen::Triplet<T>> triplet;
		for (int i = 0; i < nnz; i++)
		{
			int row;
			int col;
			T value;
			in.read((char*)&row, sizeof(int));
			in.read((char*)&col, sizeof(int));
			in.read((char*)&value, sizeof(T));
			triplet.push_back(Eigen::Triplet<T>(row, col, value));
		}

		sp.setFromTriplets(triplet.begin(), triplet.end());
	}
};
