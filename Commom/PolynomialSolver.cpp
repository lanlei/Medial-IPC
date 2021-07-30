#include "PolynomialSolver.h"

namespace PolynoimalSolver
{

	qeal computeQuadricEquation(qeal & x, qeal & a, qeal & b, qeal & c)
	{
		return a * x*x + b * x + c;
	}

	qeal computeQuadricEquation(qeal & x1, qeal & x2, qeal & a, qeal & b, qeal & c, qeal & d, qeal & e, qeal & f)
	{
		return a * x1 * x1 + b * x1 * x2 + c * x2 * x2 + d * x1 + e * x2 + f;
	}

	qeal computeQuadricEquation(qeal & x1, qeal & x2, qeal & x3, qeal & a, qeal & b, qeal & c, qeal & d, qeal & e, qeal & f, qeal & g, qeal & h, qeal& i, qeal & j)
	{
		return a * x1 * x1 + b * x2 * x2 + c * x3 * x2 + d * x1 * x2 + e * x1 * x3 + f * x2 * x3 + g * x1 + h * x2 + i * x3 + j;
	}

	qeal computeQuadricEquation(qeal & x1, qeal & x2, qeal & x3, qeal& x4, qeal & a, qeal & b, qeal & c, qeal & d, qeal & e, qeal & f, qeal & g, qeal & h, qeal & i, qeal & j, qeal & k, qeal & l, qeal & m, qeal & n, qeal & o)
	{
		return a * x1 * x1 + b * x2 * x2 + c * x3 * x3 + d * x4 * x4 + e * x1 * x2 + f * x1 * x3 + g * x1 * x4 + h * x2 * x3 + i * x2 * x4 + j * x3 * x4 + k * x1 + l * x2 + m * x3 + n * x4 + o;
	}

	qeal computeCubicEquation(qeal & x, qeal & a, qeal & b, qeal & c, qeal & d)
	{
		return a * x * x * x + b * x * x + c * x + d;
	}

	qeal computeQuarticEquation(qeal & x, qeal & a, qeal & b, qeal & c, qeal & d, qeal & e)
	{
		return a * x * x * x * x + b * x * x * x + c * x * x + d * x + e;
	}

	qeal computeQuarticEquation(qeal & x, qeal & A_2, qeal & A_1, qeal& A_0, qeal & B_2, qeal & B_1, qeal & B_0, qeal & C_2, qeal & C_1, qeal & C_0, qeal & D_2, qeal & D_1, qeal & D_0)
	{
		return (A_2 * x * x + A_1 * x + A_0) *  (B_2 * x * x + B_1 * x + B_0) - (C_2 * x * x + C_1 * x + C_0) *  (D_2 * x * x + D_1 * x + D_0);
	}

	bool solveQuadricEquation(qeal& a, qeal& b, qeal& c, qeal* x, int& countRoot, bool banSame, qeal mini, qeal maxi)
	{
		countRoot = 0;
		qeal temp_x[2];
		qeal root;
		if (IS_QEAL_ZERO(a))
		{
			if (IS_QEAL_ZERO(b))
			{
				if (IS_QEAL_ZERO(c))
				{
					temp_x[0] = 0;
					temp_x[1] = 1;
					countRoot = -1; // All real numbers are roots
				}
			}
			else
			{
				root = -c / b;
				temp_x[countRoot++] = root;
				if (!banSame)
					temp_x[countRoot++] = root;
			}
		}
		else
		{
			qeal a_ = 1.0;
			qeal b_ = b / a;
			qeal c_ = c / a;
			
			qeal delta = b_ * b_ - 4.0 * a_ * c_;
			if (IS_QEAL_ZERO(delta))
			{
				root = -b_ / (2.0 * a_);
				temp_x[countRoot++] = root;
				if (!banSame)
					temp_x[countRoot++] = root;
			}
			else if (delta > 0)
			{
				if (IS_QEAL_ZERO(c_))
				{
					root = 0;
					temp_x[countRoot++] = root;
					root = -b_ / a_;
					temp_x[countRoot++] = root;
				}
				else
				{
					qeal sdelta = sqrt(delta);
					if (b_ > 0.0)
					{
						root = (-2.0 * c_) / (b_ + sdelta);
						temp_x[countRoot++] = root;
						root = (-b_ - sdelta) / (2.0 * a_);
						temp_x[countRoot++] = root;
					}
					else
					{
						root = (-b_ + sdelta) / (2.0 * a_);
						temp_x[countRoot++] = root;
						root = (-2.0 * c_) / (b_ - sdelta);
						temp_x[countRoot++] = root;
					}
				}
			}
		}

		if (countRoot == 0)
			return false;
		if (countRoot == 2)
			if (temp_x[0] > temp_x[1])
			{
				qeal temp = temp_x[0];
				temp_x[0] = temp_x[1];
				temp_x[1] = temp;
			}
		
		int rc = 0;
		for (int i = 0; i < countRoot; i++)
		{
			if (isInValidRange(temp_x[i], mini, maxi))
			{
				x[rc++] = temp_x[i];
			}
		}
		countRoot = rc;
		if (countRoot == 0) return false;
		return true;
	}

	qeal correctQuadricEquationSolution(qeal & x, qeal & a, qeal & b, qeal & c)
	{
		int maxIter = 200;
		int iter = 0;
		qeal fx = a * x * x + b * x + c;
		fx = POLYNOMIAL_SOLVER_PRECISION(fx);
		while (fx != 0.0 && iter < maxIter)
		{
			qeal fd = 2.0 * a * x + b;
			fd = POLYNOMIAL_SOLVER_PRECISION(fd);
			if (fd == 0.0) break;
			x = x - fx / fd;
			fx = a * x * x + b * x + c;
			fx = POLYNOMIAL_SOLVER_PRECISION(fx);
			iter++;
		}
		return fx;
	}

	bool solveCubicEquation(qeal & a, qeal & b, qeal & c, qeal & d, qeal * x, int & countRoot, bool banSame, qeal mini, qeal maxi)
	{
		countRoot = 0;
		qeal temp_x[3];
		qeal root;

		if (IS_QEAL_ZERO(a))
			return solveQuadricEquation(b, c, d, x, countRoot, banSame, mini, maxi);
		else
		{
			qeal a_ = 1.0;
			qeal b_ = b / a;
			qeal c_ = c / a;
			qeal d_ = d / a;
			
			qeal A = b_ * b_ - 3.0 * a_ * c_;
			qeal B = b_ * c_ - 9.0 * a_ * d_;
			qeal C = c_ * c_ - 3.0 * b_ * d_;

			qeal delta = B * B - 4.0 * A * C;

			if (IS_QEAL_ZERO(A) && IS_QEAL_ZERO(B))
			{
				root = -b_ / (3.0 * a_);
				temp_x[countRoot++] = root;
				if (!banSame)
					temp_x[countRoot++] = root;
			}
			else if (delta > QEAL_MIN)
			{
				qeal sdelta = sqrt(delta);
				qeal Z1 = (-B - sdelta) * 0.5;
				qeal Z2 = (-B + sdelta) * 0.5;
				qeal Y1 = A * b_ + 3.0 * a_ * Z1;
				qeal Y2 = A * b_ + 3.0 * a_ * Z2;

				if (Y1 < 0.0 && Y2 < 0.0)
					root = (-b_ + pow(-Y1, 1.0 / 3.0) + pow(-Y2, 1.0 / 3.0)) / (3.0 * a_);
				else if (Y1 < 0.0 && Y2 > 0.0)
					root = (-b_ + pow(-Y1, 1.0 / 3.0) - pow(Y2, 1.0 / 3.0)) / (3.0 * a_);
				else if (Y1 > 0.0 && Y2 < 0.0)
					root = (-b_ - pow(Y1, 1.0 / 3.0) + pow(-Y2, 1.0 / 3.0)) / (3.0 * a_);
				else
					root = (-b_ - pow(Y1, 1.0 / 3.0) - pow(Y2, 1.0 / 3.0)) / (3.0 * a_);
				temp_x[countRoot++] = root;
			}
			else if (IS_QEAL_ZERO(delta))
			{
				if (!IS_QEAL_ZERO(A))
				{
					qeal K = B / A;
					root = -b_ / a_ + K;

					temp_x[countRoot++] = root;
					root = -0.5 * K;
					temp_x[countRoot++] = root;
					if (!banSame)
						temp_x[countRoot++] = root;
				}
			}
			else
			{
				const qeal CERROR = MIN_VALUE;
				if (A > 0.0)
				{
					qeal T = (2.0 * A * b_ - 3.0 * a_ * B) / (2.0 * pow(A, 3.0 / 2.0));

					if (T > 1.0)
					{
						if (T < 1.0 + CERROR)
							T = 1.0;
						else return false;
					}
					else if (T < -1.0)
					{
						if (T > -1.0 - CERROR)
							T = -1.0;
						else return false;
					}

					qeal theta = acos(T);
					root = (-b_ - 2.0 * sqrt(A) * cos(theta / 3.0)) / (3.0 * a_);
					temp_x[countRoot++] = root;
					root = (-b_ + sqrt(A) * (cos(theta / 3.0) + sqrt(3.0) * sin(theta / 3.0))) / (3.0 * a_);
					temp_x[countRoot++] = root;
					root = (-b_ + sqrt(A) * (cos(theta / 3.0) - sqrt(3.0) * sin(theta / 3.0))) / (3.0 * a_);
					temp_x[countRoot++] = root;
				}
			}
		}
		if (countRoot == 0) return false;
		if (countRoot >= 2)
			std::sort(temp_x, temp_x + countRoot);

		int rc = 0;
		for (int i = 0; i < countRoot; i++)
		{
			if (isInValidRange(temp_x[i], mini, maxi))
			{
				x[rc++] = temp_x[i];
			}
		}
		countRoot = rc;
		if (countRoot == 0) return false;
		return true;
	}

	qeal correctCubicEquationSolution(qeal & x, qeal & a, qeal & b, qeal & c, qeal & d)
	{
		int maxIter = 200;
		int iter = 0;
		x = 0;
		qeal fx = a * x * x * x + b * x * x+ c * x + d;
		while (fabs(fx) > 1e-6 && iter < maxIter)
		{
			qeal fd = 3.0 * a * x * x+ 2.0 * b * x + c;
			fd = POLYNOMIAL_SOLVER_PRECISION(fd);
			if (fd == 0.0) break;
			x = x - fx / fd;
			fx = a * x * x * x + b * x * x + c * x + d;
			fx = POLYNOMIAL_SOLVER_PRECISION(fx);
			iter++;
		}
		return fx;
	}

	bool solveQuarticEquation(qeal & a, qeal & b, qeal & c, qeal & d, qeal & e, qeal * x, int & countRoot, bool banSame, qeal mini, qeal maxi)
	{
		countRoot = 0;
		qeal temp_x[3];
		qeal root;

		if (IS_QEAL_ZERO(a))
			return solveCubicEquation(b, c, d, e, x, countRoot, banSame, mini, maxi);
		else
		{
			qeal b_ = b / a;
			qeal c_ = c / a;
			qeal d_ = d / a;
			qeal e_ = e / a;
			
			qeal J0 = 8.0;
			qeal J1 = -4.0 * c_;
			qeal J2 = 2.0 * (b_ * d_ - 4.0 * e_);
			qeal J3 = -e_ * (b_ * b_ - 4.0 * c_) - d_ * d_;
			qeal Jx[3];
			int Jc = 0;
			if (!solveCubicEquation(J0, J1, J2, J3, Jx, Jc, false))
				return false;
			qeal y, M, N;
			qeal x1[2], x2[2];
			int countRoot1, countRoot2, i;
			qeal MSquareTemp, MSquare, yTemp;
			if (Jc != 0)
			{
				y = Jx[0];

				MSquare = 8.0 * y + b_ * b_ - 4.0 * c_;
				for (i = 1; i < Jc; i++)
				{
					yTemp = Jx[i];
					MSquareTemp = 8.0 * yTemp + b_ * b_ - 4.0 * c_;
					if (MSquareTemp > MSquare)
					{
						MSquare = MSquareTemp;
						y = yTemp;
					}
				}

				if (MSquare > 0.0)
				{
					M = sqrt(MSquare);
					N = b_ * y - d_;

					qeal K0 = 2.0;
					qeal K1 = b_ + M;
					qeal K2 = 2.0 * (y + N / M);
					solveQuadricEquation(K0, K1, K2, x1, countRoot1, true, mini, maxi);
					qeal L0 = 2.0;
					qeal L1 = b_ - M;
					qeal L2 = 2.0 * (y - N / M);
					solveQuadricEquation(L0, L1, L2, x2, countRoot2, true, mini, maxi);
 					for (int j = 0; j < countRoot1; j++)
						temp_x[countRoot++] = x1[j];
					for (int j = 0; j < countRoot2; j++)
						temp_x[countRoot++] = x2[j];
				}
			}
		}

		if (countRoot == 0) return false;
		if (countRoot >= 2)
			std::sort(temp_x, temp_x + countRoot);

		int rc = 0;
		for (int i = 0; i < countRoot; i++)
		{
			if (isInValidRange(temp_x[i], mini, maxi))
			{
				x[rc++] = temp_x[i];
			}
		}
		countRoot = rc;
		if (countRoot == 0) return false;
		return true;
	}

	qeal correctQuarticEquationSolution(qeal & x, qeal & a, qeal & b, qeal & c, qeal & d, qeal & e)
	{
		int maxIter = 200;
		int iter = 0;
		qeal fx = a * x * x * x * x + b * x * x * x+ c * x * x+ d * x + e;
		while (fabs(fx) > 1e-6 && iter < maxIter)
		{
			qeal fd = 4.0 * a * x * x * x + 3.0 * b * x * x + 2.0 * c * x + d;
			fd = POLYNOMIAL_SOLVER_PRECISION(fd);
			if (fd == 0.0) break;
			x = x - fx / fd;
			fx = a * x * x * x * x + b * x * x * x + c * x * x + d * x + e;
			fx = POLYNOMIAL_SOLVER_PRECISION(fx);
			iter++;
		}
		return fx;
	}

	void solveQuadricNEq(qeal & a, qeal & b, qeal & c, qeal * x, int & countRoot, ValidRange & posRange, ValidRange & negRange, qeal mini, qeal maxi)
	{
		countRoot = 0;
		if (IS_QEAL_ZERO(a) && IS_QEAL_ZERO(b))
		{
			if (c > 0.0)
				negRange.setNull();
			else if (c < 0.0)
				posRange.setNull();
			else
			{
				x[0] = mini; x[1] = maxi;
				countRoot = 2;
			}
			return;
		}

		bool hasSolution = solveQuadricEquation(a, b, c, x, countRoot, true, mini, maxi);

		if (hasSolution)
		{
			bool posInsert = false;
			bool negInsert = false;

			qeal lx = mini, hx = maxi;
			if (lx != x[0])
			{
				hx = x[0];
				qeal midx = lx + (hx - lx) / 2;
				qeal v = computeQuadricEquation(midx, a, b, c);
				v = POLYNOMIAL_SOLVER_PRECISION(v);
				if (v >= 0.0)
				{
					posRange.insert(Vector2(lx, hx));
					posInsert = true;
				}
				else
				{
					negRange.insert(Vector2(lx, hx));
					negInsert = true;
				}
					
			}
			lx = x[0];
			if (countRoot > 1)
				hx = x[1];
			else hx = maxi;
			qeal midx = lx + (hx - lx) / 2;
			qeal v = computeQuadricEquation(midx, a, b, c);
			v = POLYNOMIAL_SOLVER_PRECISION(v);
			if (v >= 0.0)
			{
				posRange.insert(Vector2(lx, hx));
				posInsert = true;
			}
			else
			{
				negRange.insert(Vector2(lx, hx));
				negInsert = true;
			}
			if (!posInsert)posRange.setNull();
			if (!negInsert)negRange.setNull();
		}
		else
		{
			qeal midx = mini + (maxi - mini) / 2.0;
			qeal v = computeQuadricEquation(midx, a, b, c);

			v = POLYNOMIAL_SOLVER_PRECISION(v);
			if (v >= 0.0)
			{
				posRange.insert(Vector2(mini, maxi));
				negRange.setNull();
			}
			else
			{
				negRange.insert(Vector2(mini, maxi));
				posRange.setNull();
			}
		}
	}

	void solveQuarticNEq(qeal & a, qeal & b, qeal & c, qeal & d, qeal & e, qeal * x, int & countRoot, ValidRange & posRange, ValidRange & negRange, qeal mini, qeal maxi)
	{
		countRoot = 0;
		if (IS_QEAL_ZERO(a) && IS_QEAL_ZERO(b))
		{
			solveQuadricNEq(c, d, e, x, countRoot, posRange, negRange, mini, maxi);
			return;
		}

		bool hasSolution = solveQuarticEquation(a, b, c, d, e, x, countRoot, true, mini, maxi);

		if (hasSolution)
		{
			bool posInsert = false;
			bool negInsert = false;
			qeal lx = mini, hx = maxi;
			if (lx != x[0])
			{
				hx = x[0];
				qeal midx = lx + (hx - lx) / 2;
				qeal v = computeQuarticEquation(midx, a, b, c, d, e);
				v = POLYNOMIAL_SOLVER_PRECISION(v);
				if (v >= 0.0)
				{
					posRange.insert(Vector2(lx, hx));
					posInsert = true;
				}					
				else
				{
					negRange.insert(Vector2(lx, hx));
					negInsert = true;
				}					
			}

			//
			lx = x[0];
			if (countRoot > 1)
				hx = x[1];
			else hx = maxi;
			qeal midx = lx + (hx - lx) / 2;
			qeal v = computeQuarticEquation(midx, a, b, c, d, e);
			v = POLYNOMIAL_SOLVER_PRECISION(v);
			if (v >= 0.0)
			{
				posRange.insert(Vector2(lx, hx));
				posInsert = true;
			}
			else
			{
				negRange.insert(Vector2(lx, hx));
				negInsert = true;
			}
			
			if (!posInsert)posRange.setNull();
			if (!negInsert)negRange.setNull();
		}
		else
		{
			qeal midx = mini + (maxi - mini) / 2.0;
			qeal v = computeQuarticEquation(midx, a, b, c, d, e);
			v = POLYNOMIAL_SOLVER_PRECISION(v);
			if (v >= 0.0)
			{
				posRange.insert(Vector2(mini, maxi));
				negRange.setNull();
			}
			else
			{
				negRange.insert(Vector2(mini, maxi));
				posRange.setNull();
			}
		}
	}
}
