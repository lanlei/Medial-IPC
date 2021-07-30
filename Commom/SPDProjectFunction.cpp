#include "SPDProjectFunction.h"

void makePD(MatrixX& symMtr)
{
	Eigen::SelfAdjointEigenSolver<MatrixX> eigenSolver(symMtr);
	if (eigenSolver.eigenvalues()[0] >= 0) {
		return;
	}
	int size = symMtr.rows();
	MatrixX D = eigenSolver.eigenvalues().asDiagonal();

	for (int i = 0; i < size; ++i) {
		if (D.diagonal()[i] < 0) D.diagonal()[i] = 0;
		else break;
	}

	symMtr = eigenSolver.eigenvectors() * D * eigenSolver.eigenvectors().transpose();
}