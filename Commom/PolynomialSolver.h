#pragma once
#ifndef POLYNOMIAL_SOLVER_H
#include "MatrixCore.h"
#include <stdio.h>
#include <sstream>
#include <queue>
#include<iostream>

namespace PolynoimalSolver
{
#define POLYNOMIAL_SOLVER_PRECISION(d) (abs(d) > 1e-12 ? d : 0)
	
	// f(x) = ax^2 + bx + c
	qeal computeQuadricEquation(qeal& x, qeal& a, qeal& b, qeal& c);

	//f(x1,x2) = ax1 ^ 2 + bx1x2 + cx2^2 + dx1 + ex2 + f
	qeal computeQuadricEquation(qeal& x1, qeal& x2, qeal& a, qeal& b, qeal& c, qeal& d, qeal& e, qeal& f);

	//f(x1,x2, x3) = ax1 ^ 2 + bx2^2 + cx3^2 + dx1x2 +ex1x3 + fx2x3 + gx1 + hx2 +ix3 + j
	qeal computeQuadricEquation(qeal & x1, qeal & x2, qeal & x3, qeal & a, qeal & b, qeal & c, qeal & d, qeal & e, qeal & f, qeal & g, qeal & h, qeal& i, qeal & j);

	//f(x1,x2, x3, x4) = ax1 ^ 2 + bx2^2 + cx3^2 + dx4^2 + ex1x2 + fx1x3 + gx1x4 + hx2x3 + ix2x4 + jx3x4 + kx1 + lx2 + mx3 + nx4 + o
	qeal computeQuadricEquation(qeal & x1, qeal & x2, qeal & x3, qeal& x4, qeal & a, qeal & b, qeal & c, qeal & d, qeal & e, qeal & f, qeal & g, qeal & h, qeal & i, qeal & j, qeal & k, qeal & l, qeal & m, qeal & n, qeal & o);

	//f(x) = ax^3 + bx^2 + cx + d;
	qeal computeCubicEquation(qeal& x, qeal& a, qeal& b, qeal& c, qeal& d);

	//f(x) = ax^4 + bx^3 + cx^2 + dx + e;
	qeal computeQuarticEquation(qeal& x, qeal& a, qeal& b, qeal& c, qeal& d, qeal& e);

	qeal computeQuarticEquation(qeal & x, qeal & A_2, qeal & A_1, qeal& A_0, qeal & B_2, qeal & B_1, qeal & B_0, qeal & C_2, qeal & C_1, qeal & C_0, qeal & D_2, qeal & D_1, qeal & D_0);

	// get real roots of the equation ax ^ 2 + bx + c = 0;
	// countRoot is the number of return roots; 
	// if banSame = true, remove the same roots;
	// all roots must be on the range [mini, maxi];
	bool solveQuadricEquation(qeal& a, qeal& b, qeal& c, qeal* x, int& countRoot, bool banSame = false, qeal mini = -QEAL_MAX, qeal maxi = QEAL_MAX);

	qeal correctQuadricEquationSolution(qeal& x, qeal& a, qeal& b, qeal& c);

	// get real roots of the equation ax ^ 3 + bx^2 + cx + d = 0;
	// countRoot is the number of return roots; 
	// if banSame = true, remove the same roots;
	// all roots must be on the range [mini, maxi];
	bool solveCubicEquation(qeal& a, qeal& b, qeal& c, qeal& d, qeal* x, int& countRoot, bool banSame = false, qeal mini = -QEAL_MAX, qeal maxi = QEAL_MAX);

	qeal correctCubicEquationSolution(qeal& x, qeal& a, qeal& b, qeal& c, qeal& d);

	// get real roots of the equation ax^4 + bx^3 + cx^2 + dx + e = 0;
	// countRoot is the number of the return roots; 
	// if banSame = true, remove the same roots;
	// all roots must be on the range [mini, maxi];
	bool solveQuarticEquation(qeal& a, qeal& b, qeal& c, qeal& d, qeal& e, qeal* x, int& countRoot, bool banSame = false, qeal mini = -QEAL_MAX, qeal maxi = QEAL_MAX);

	qeal correctQuarticEquationSolution(qeal& x, qeal& a, qeal& b, qeal& c, qeal& d, qeal& e);

	struct ValidRange
	{
		ValidRange(Vector2 valid = Vector2(0.0, 1.0))
		{
			_nullSet = false;
			_valid = valid;
			set.push(valid);
		}

		bool isNull() { return _nullSet; }
		void setNull() {
			_nullSet = true;
			while (!set.empty()) set.pop();
		}

		bool insert(Vector2& s0)
		{
			if (_nullSet)
				return false;
			int size = set.size();
			while (size > 0)
			{
				Vector2 s1 = set.front();
				set.pop();
				if (s0.data()[1] < s1.data()[0] || s0.data()[0] > s1.data()[1])
				{
					size--;
					continue;
				}
				Vector2 newOne;

				if (s0.data()[0] >= s1.data()[0])
					newOne.data()[0] = s0.data()[0];
				else newOne.data()[0] = s1.data()[0];
				if (s0.data()[1] <= s1.data()[1])
					newOne.data()[1] = s0.data()[1];
				else newOne.data()[1] = s1.data()[1];
				set.push(newOne);
				size--;
			}
			if (set.empty())
			{
				setNull();
				return false;
			}
			return true;
		}

		bool insert(std::vector<Vector2>& outter)
		{
			if (_nullSet)
				return false;
			int size = set.size();

			while (size > 0)
			{
				Vector2 s1 = set.front();
				set.pop();

				std::vector<Vector2> overlap;
				for (int i = 0; i < outter.size(); i++)
				{
					Vector2 s0 = outter[i];
					if (s0.data()[1] < s1.data()[0] || s0.data()[0] > s1.data()[1])
						continue;
					Vector2 newOne;
					if (s0.data()[0] >= s1.data()[0])
						newOne.data()[0] = s0.data()[0];
					else newOne.data()[0] = s1.data()[0];
					if (s0.data()[1] <= s1.data()[1])
						newOne.data()[1] = s0.data()[1];
					else newOne.data()[1] = s1.data()[1];
					overlap.push_back(newOne);
				}
				if (overlap.size() == 0)
				{
					size--;
					continue;
				}
				else
				{
					for(int i = 0; i < overlap.size(); i++)
						set.push(overlap[i]);
					size--;
				}

			}
			if (set.empty())
			{
				setNull();
				return false;
			}
			return true;
		}

		bool check(qeal& p)
		{
			if (_nullSet)
				return false;
			int size = set.size();
			while (size > 0)
			{
				Vector2 s1 = set.front();
				set.pop();
				set.push(s1);
				size--;
				if (p <= s1.data()[0] && p >= s1.data()[1])
					return true;
			}
			return false;
		}

		std::queue<Vector2> set;
	protected:
		bool _nullSet;
		Vector2 _valid;
	};

	void solveQuadricNEq(qeal& a, qeal& b, qeal& c, qeal* x, int& countRoot, ValidRange& posRange, ValidRange& negRange, qeal mini = 0.0, qeal maxi = 1.0);

	void solveQuarticNEq(qeal& a, qeal& b, qeal& c, qeal& d, qeal& e, qeal* x, int& countRoot, ValidRange& posRange, ValidRange& negRange, qeal mini = 0.0, qeal maxi = 1.0);
}



#endif