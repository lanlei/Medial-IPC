#include "Simulator\CollisionDetection\CollisionDetectionMedialMesh.h"

namespace CDMM
{
#define NO_COLLIDE 0
#define INTERSECTION 1
#define PENETRATION 2

	qeal valueOfQuadircSurface2D(const qeal x, const qeal y, const qeal A, const qeal B, const qeal C, const qeal D, const qeal E, const qeal F)
	{
		return A * x * x + B * x * y + C * y * y + D * x + E * y + F;
	}

	void genQuarticCoeffs(qeal & P_4, qeal & P_3, qeal & P_2, qeal & P_1, qeal & P_0, qeal A_2, qeal A_1, qeal A_0, qeal B_2, qeal B_1, qeal B_0, qeal C_2, qeal C_1, qeal C_0, qeal D_2, qeal D_1, qeal D_0, qeal S)
	{
		P_4 += (A_2 * B_2 - C_2 * D_2) * S;
		P_3 += (A_2 * B_1 + A_1 * B_2 - C_2 * D_1 - C_1 * D_2) * S;
		P_2 += (A_0 * B_2 + A_2 * B_0 + A_1 * B_1 - C_0 * D_2 - C_2 * D_0 - C_1 * D_1) * S;
		P_1 += (A_0 * B_1 + A_1 * B_0 - C_0 * D_1 - C_1 * D_0) * S;
		P_0 += (A_0 * B_0 - C_0 * D_0) * S;
	}

	void genSexticCoeffs(qeal & P_6, qeal & P_5, qeal & P_4, qeal & P_3, qeal & P_2, qeal & P_1, qeal & P_0, qeal A_2, qeal A_1, qeal A_0, qeal B_2, qeal B_1, qeal B_0, qeal C_2, qeal C_1, qeal C_0, qeal S)
	{
		P_6 += A_2 * B_2 * C_2 * S;
		P_5 += (A_2 * B_1 * C_2 + A_1 * B_2 * C_2 + A_2 * B_2 * C_1) * S;
		P_4 += (A_2 * B_0 * C_2 + A_2 * B_2 * C_0 + A_0 * B_2 * C_2 + A_2 * B_1 * C_1 + A_1 * B_2*C_1 + A_1 * B_1 * C_2) * S;
		P_3 += (A_2 * B_1 * C_0 + A_1 * B_2 * C_0 + A_0 * B_2 * C_1 + A_2 * B_0 * C_1 + A_1 * B_1 * C_1 + A_0 * B_1*C_2 + A_1 * B_0*C_2) * S;
		P_2 += (A_2 * B_0 * C_0 + A_0 * B_2 * C_0 + A_0 * B_0 * C_2 + A_1 * B_1*C_0 + A_0 * B_1 * C_1 + A_1 * B_0 * C_1) * S;
		P_1 += (A_1 * B_0 * C_0 + A_0 * B_1 * C_0 + A_0 * B_0 * C_1) * S;
		P_0 += A_0 * B_0 * C_0 * S;
	}

	qeal SolveSexticEqForFTC(qeal & W_6, qeal & W_5, qeal & W_4, qeal & W_3, qeal & W_2, qeal & W_1, qeal & W_0, int banRootsNum, qeal* banRoots, RootRange & range, qeal norrow)
	{
		qeal ftc = 1.0;
		while (!range.set.empty())
		{
			qeal ft = 1.0;
			Vector2 r = range.set.front();
			range.set.pop();
			int count = 0;
			qeal x[6];
			if (W_6 == 0.0 && W_5 == 0.0)
			{				
				if (solveQuarticEquation(W_4, W_3, W_2, W_1, W_0, x, count, true, r.data()[0], r.data()[1]))
				{					
					qeal fx = W_6 * x[0] * x[0] * x[0] * x[0] * x[0] * x[0] + W_5 * x[0] * x[0] * x[0] * x[0] * x[0] + W_4 * x[0] * x[0] * x[0] * x[0] + W_3 * x[0] * x[0] * x[0] + W_2 * x[0] * x[0] + W_1 * x[0] + W_0;
					
					if (fx >= 0.0 && fx <= std::max(norrow, CDMM_NORROW))
						ft = x[0];
					else if (fx < 0.0)
					{
						ft = x[0];  // _distance still excess 0.0
					}

					if(ft == 1.0)
						ft = NewtonSolverSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, r, norrow);
				}
				else ft = NewtonSolverSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, r, norrow);
			}
			else ft = NewtonSolverSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, r, norrow);

		//	ft = NewtonSolverSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, r, norrow);

			bool isBan = false;
			for (int j = 0; j < banRootsNum; j++)
			{
				if (IS_QEAL_ZERO(ft - banRoots[j]))
				{
					isBan = true;
					break;
				}
			}
			if (isBan) continue;
			if(ft < ftc)
				ftc = ft;
		}
		return ftc;
	}

	qeal NewtonSolverSexticEqForFTC(qeal & W_6, qeal & W_5, qeal & W_4, qeal & W_3, qeal & W_2, qeal & W_1, qeal & W_0, Vector2& r, qeal norrow)
	{		
		qeal ftc = 1.0;
		int maxIter = 200;
		int iter = 0;
		qeal x = r.data()[0];
		qeal fx, fd;

		if (x > ftc) return ftc;

		if (x > 0.0)
		{
			fx = W_6 * x * x * x * x  * x * x + W_5 * x * x * x  * x * x + W_4 * x * x  * x * x + W_3 * x  * x * x + W_2 * x * x + W_1 * x + W_0;
			if (fx >= 0.0 && fx <= std::max(norrow, CDMM_NORROW))
			{
				ftc = x;
				return ftc;
			}
			else if (fx < 0.0)
			{
				ftc = x;  // _distance still excess 0.0	
				return ftc;
			}
		}


		do
		{
			fd = 6.0 * W_6 * x * x * x * x  * x + 5.0 * W_5 * x * x * x  * x + 4.0 * W_4 * x * x  * x + 3.0 * W_3 * x  * x + 2.0 * W_2 * x + W_1;
			fd = POLYNOMIAL_SOLVER_PRECISION(fd);
			if (fd == 0.0)
				break;
			x = x - fx / fd;
			if (x > r.data()[1] || x < r.data()[0])
				x = r.data()[1];
			fx = W_6 * x * x * x * x  * x * x + W_5 * x * x * x  * x * x + W_4 * x * x  * x * x + W_3 * x  * x * x + W_2 * x * x + W_1 * x + W_0;
			iter++;
		} while (abs(fx) > CDMM_NewtonSolverPrecision && iter < maxIter);
		if (x >= r.data()[0] && x <= r.data()[1])
		{
			if (fx >= 0.0 && fx <= std::max(norrow, CDMM_NORROW))
				ftc = x;
			else if (fx < 0.0)
			{
				ftc = x;  // _distance still excess 0.0
			}
		}
		return ftc;
	}

	bool dcd(Vector3 c11, qeal r11, Vector3 c12, qeal r12, Vector3 c21, qeal r21, Vector3 c22, qeal r22, bool space_triangle, qeal norrow)
	{
		qeal R1 = r11 - r12;
		qeal R2 = r21 - r22;
		qeal R3 = r12 + r22;
		Vector3 C1 = c11 - c12;
		Vector3 C2 = c22 - c21;
		Vector3 C3 = c12 - c22;
		return dcd(C1, C2, C3, R1, R2, R3, space_triangle, norrow);
	}

	bool dcd(Vector3 C1, Vector3 C2, Vector3 C3, qeal & R1, qeal & R2, qeal & R3, bool space_triangle, qeal norrow)
	{
		qeal A = C1.dot(C1) - R1 * R1;
		qeal B = 2.0 * (C1.dot(C2) - R1 * R2);
		qeal C = C2.dot(C2) - R2 * R2;
		qeal D = 2.0 * (C1.dot(C3) - R1 * R3);
		qeal E = 2.0 * (C2.dot(C3) - R2 * R3);
		qeal F = C3.dot(C3) - R3 * R3;

		return dcd(A, B, C, D, E, F, space_triangle, norrow);
	}

	bool dcd(qeal A, qeal B, qeal C, qeal D, qeal E, qeal F, bool space_triangle, qeal norrow)
	{
		qeal min = 1;
		//case 1 x = 0, y = 0
		if (D >= 0 && E >= 0)
		{
			qeal f = valueOfQuadircSurface2D(0, 0, A, B, C, D, E, F);
			if (f < min)
			{
				min = f;
			}
		}
		//case 2: x = 0, y != 0,1
		{
			qeal E2C = -1 * (E / (2 * C));
			if (E2C > 0 && E2C < 1)
			{
				qeal DB = B * E2C + D;
				if (DB >= 0)
				{
					qeal f = valueOfQuadircSurface2D(0, E2C, A, B, C, D, E, F);
					if (f < min)
					{
						min = f;
					}
				}
			}
		}
		// case 3 x = 0, y = 1;
		if ((B + D) >= 0 && (2 * C + E) <= 0)
		{
			qeal f = valueOfQuadircSurface2D(0, 1, A, B, C, D, E, F);
			if (f < min)
			{
				min = f;
			}
		}
		//case 4 x != 0, 1, y = 0
		{
			if (!IS_QEAL_ZERO(A))
			{
				qeal D2A = -1 * (D / (2 * A));
				if (D2A > 0 && D2A < 1)
				{
					qeal EB = B * D2A + E;
					if (EB >= 0)
					{
						qeal f = valueOfQuadircSurface2D(D2A, 0, A, B, C, D, E, F);
						if (f < min)
						{
							min = f;
						}
					}
				}
			}
		}
		// case 5 x != 0,1 y != 0,1
		{
			if (!IS_QEAL_ZERO(4 * A*C - B * B))
			{
				qeal x = (B*E - 2 * C*D) / (4 * A*C - B * B);
				qeal y = (B*D - 2 * A*E) / (4 * A*C - B * B);
				if (x > 0 && x < 1 && y > 0 && y < 1)
				{
					if (!space_triangle)
					{
						qeal f = valueOfQuadircSurface2D(x, y, A, B, C, D, E, F);
						if (f < min)
						{
							min = f;
						}
					}
					else
					{
						if ((x + y) <= 1)
						{
							qeal f = valueOfQuadircSurface2D(x, y, A, B, C, D, E, F);
							if (f < min)
							{
								min = f;
							}
						}
					}
				}
			}
		}
		// case 6 x != 0,1 y = 1
		{
			if (!space_triangle)
			{
				if (!IS_QEAL_ZERO(A))
				{
					qeal x = -1 * ((B + D) / (2 * A));
					qeal CBE = 2 * C - (B*B + B * D) / (2 * A) + E;
					if (x > 0 && x < 1 && CBE <= 0)
					{
						qeal f = valueOfQuadircSurface2D(x, 1, A, B, C, D, E, F);
						if (f < min)
						{
							min = f;
						}
					}
				}
			}
		}

		// case 7 x =1 y = 0
		{
			if ((-1 * (2 * A + D)) >= 0 && (B + E) >= 0)
			{
				qeal f = valueOfQuadircSurface2D(1, 0, A, B, C, D, E, F);
				if (f < min)
				{
					min = f;
				}
			}
		}

		// case 8 x =1 y != 0,1
		{
			if (!space_triangle)
			{
				qeal y = -1 * ((B + E) / (2 * C));
				qeal ABD = 2 * A - (B*B + B * E) / (2 * C) + D;
				if (y > 0 && y < 1 && ABD <= 0)
				{
					qeal f = valueOfQuadircSurface2D(1, y, A, B, C, D, E, F);
					if (f < min)
					{
						min = f;
					}
				}
			}
		}

		// case 9 x =1, y = 1
		{
			if (!space_triangle)
			{
				qeal ABD = -1 * (2 * A + B + D);
				qeal CBE = -1 * (2 * C + B + E);
				if (ABD >= 0 && CBE >= 0)
				{
					qeal f = valueOfQuadircSurface2D(1, 1, A, B, C, D, E, F);
					if (f < min)
					{
						min = f;
					}
				}
			}
		}

		min -= norrow;
		min = Check_QEAL_ZERO(min);

		if (min > 0.0)
			return false;

		return true;
	}

	bool ccd(Vector3 & C1, Vector3 & V1, Vector3 & C2, Vector3 & V2, Vector3 & C3, Vector3 & V3, qeal & R1, qeal & R2, qeal & R3, qeal & ftc, bool space_triangle, qeal norrow)
	{
		bool collide = false;
		qeal pftc = 1.0; // potential ftc
		
		qeal JA_2, JA_1, JA_0;
		JA_2 = V1.dot(V1); JA_1 = 2.0 * V1.dot(C1); JA_0 = C1.dot(C1) - R1 * R1;
		qeal JB_2, JB_1, JB_0;
		JB_2 = 2.0 * V1.dot(V2); JB_1 = 2.0 * (V1.dot(C2) + V2.dot(C1)); JB_0 = 2.0 * (C1.dot(C2) - R1 * R2);
		qeal JC_2, JC_1, JC_0;
		JC_2 = V2.dot(V2); JC_1 = 2.0 * V2.dot(C2); JC_0 = C2.dot(C2) - R2 * R2;
		qeal JD_2, JD_1, JD_0;
		JD_2 = 2.0 * V1.dot(V3); JD_1 = 2.0 * (V1.dot(C3) + V3.dot(C1)); JD_0 = 2.0 * (C1.dot(C3) - R1 * R3);
		qeal JE_2, JE_1, JE_0;
		JE_2 = 2.0 * V2.dot(V3); JE_1 = 2.0 * (V2.dot(C3) + V3.dot(C2)); JE_0 = 2.0 * (C2.dot(C3) - R2 * R3);
		qeal JF_2, JF_1, JF_0;
		JF_2 = V3.dot(V3); JF_1 = 2.0 * V3.dot(C3); JF_0 = C3.dot(C3) - R3 * R3 - norrow;


		JA_2 = Check_QEAL_ZERO(JA_2); JA_1 = Check_QEAL_ZERO(JA_1); JA_0 = Check_QEAL_ZERO(JA_0);
		JB_2 = Check_QEAL_ZERO(JB_2); JB_1 = Check_QEAL_ZERO(JB_1); JB_0 = Check_QEAL_ZERO(JB_0);
		JC_2 = Check_QEAL_ZERO(JC_2); JC_1 = Check_QEAL_ZERO(JC_1); JC_0 = Check_QEAL_ZERO(JC_0);
		JD_2 = Check_QEAL_ZERO(JD_2); JD_1 = Check_QEAL_ZERO(JD_1); JD_0 = Check_QEAL_ZERO(JD_0);
		JE_2 = Check_QEAL_ZERO(JE_2); JE_1 = Check_QEAL_ZERO(JE_1); JE_0 = Check_QEAL_ZERO(JE_0);
		JF_2 = Check_QEAL_ZERO(JF_2); JF_1 = Check_QEAL_ZERO(JF_1); JF_0 = Check_QEAL_ZERO(JF_0);

#ifdef CHECK_MEDIAL_PRIMITIVE_VALID
		qeal Aroots[2];
		int ArootCount = 0;
		assert(solveQuadricEquation(JA_2, JA_1, JA_0, Aroots, ArootCount, true, 0.0, 1.0) == false);
		qeal Croots[2];
		int CrootCount = 0;
		assert(solveQuadricEquation(JC_2, JC_1, JC_0, Croots, CrootCount, true, 0.0, 1.0) == false);
#endif 
		// dcd when t = 1.0;
		if (dcd(C1 + V1, C2 + V2, C3 + V3, R1, R2, R3, space_triangle, norrow))
			collide = true; // check all cases;
		//
		qeal ftc_endpoints = 1.0;
		bool collide_endpoints = checkEndpointAlphaBetaCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_endpoints, space_triangle);
		if (collide_endpoints && ftc_endpoints <= ftc)
		{
			collide = true; 
			ftc = ftc_endpoints;
		}

		qeal ftc_a0 = 1.0;
		bool collide_a0 = checkAlphaIsZeroCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_a0);
		if (collide_a0 && ftc_a0 <= ftc)
		{
			collide = true;
			ftc = ftc_a0;
		}

		qeal ftc_b0 = 1.0;
		bool collide_b0 = checkBetaIsZeroCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_b0);
		if (collide_b0 && ftc_b0 <= ftc)
		{
			collide = true;
			ftc = ftc_b0;
		}

		if (space_triangle)
		{
			qeal ftc_ab1 = 1.0;
			bool collide_ab1 = checkAlphaPlusBetaIsOneCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_ab1);
			if (collide_ab1 && ftc_ab1 <= ftc)
			{
				collide = true;
				ftc = ftc_ab1;
			}
		}
		else
		{
			qeal ftc_a1 = 1.0;
			bool collide_a1 = checkAlphaIsOneCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_a1);
			if (collide_a1 && ftc_a1 <= ftc)
			{
				collide = true;
				ftc = ftc_a1;
			}

			qeal ftc_b1 = 1.0;
			bool collide_b1 = checkBetaIsOneCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_b1);
			if (collide_b1 && ftc_b1 <= ftc)
			{
				collide = true;
				ftc = ftc_b1;
			}
		}

		qeal ftc_ab = 1.0;
		bool collide_ab = checkAlphaBetaCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_ab, space_triangle, norrow);
		if (collide_ab && ftc_ab <= ftc)
		{
			collide = true;
			ftc = ftc_ab;
		}

		if (ftc == 0.0)
			return dcd(C1, C2, C3, R1, R2, R3, space_triangle, norrow);
		return collide;
	}

	bool ccd2(Vector3 & C1, Vector3 & V1, Vector3 & C2, Vector3 & V2, Vector3 & C3, Vector3 & V3, qeal & R1, qeal & R2, qeal & R3, qeal & ftc, bool space_triangle, qeal norrow)
	{
		qeal JA_2, JA_1, JA_0;
		JA_2 = V1.dot(V1); JA_1 = 2.0 * V1.dot(C1); JA_0 = C1.dot(C1) - R1 * R1;
		qeal JB_2, JB_1, JB_0;
		JB_2 = 2.0 * V1.dot(V2); JB_1 = 2.0 * (V1.dot(C2) + V2.dot(C1)); JB_0 = 2.0 * (C1.dot(C2) - R1 * R2);
		qeal JC_2, JC_1, JC_0;
		JC_2 = V2.dot(V2); JC_1 = 2.0 * V2.dot(C2); JC_0 = C2.dot(C2) - R2 * R2;
		qeal JD_2, JD_1, JD_0;
		JD_2 = 2.0 * V1.dot(V3); JD_1 = 2.0 * (V1.dot(C3) + V3.dot(C1)); JD_0 = 2.0 * (C1.dot(C3) - R1 * R3);
		qeal JE_2, JE_1, JE_0;
		JE_2 = 2.0 * V2.dot(V3); JE_1 = 2.0 * (V2.dot(C3) + V3.dot(C2)); JE_0 = 2.0 * (C2.dot(C3) - R2 * R3);
		qeal JF_2, JF_1, JF_0;
		JF_2 = V3.dot(V3); JF_1 = 2.0 * V3.dot(C3); JF_0 = C3.dot(C3) - R3 * R3 - norrow;

		JA_2 = Check_QEAL_ZERO(JA_2); JA_1 = Check_QEAL_ZERO(JA_1); JA_0 = Check_QEAL_ZERO(JA_0);
		JB_2 = Check_QEAL_ZERO(JB_2); JB_1 = Check_QEAL_ZERO(JB_1); JB_0 = Check_QEAL_ZERO(JB_0);
		JC_2 = Check_QEAL_ZERO(JC_2); JC_1 = Check_QEAL_ZERO(JC_1); JC_0 = Check_QEAL_ZERO(JC_0);
		JD_2 = Check_QEAL_ZERO(JD_2); JD_1 = Check_QEAL_ZERO(JD_1); JD_0 = Check_QEAL_ZERO(JD_0);
		JE_2 = Check_QEAL_ZERO(JE_2); JE_1 = Check_QEAL_ZERO(JE_1); JE_0 = Check_QEAL_ZERO(JE_0);
		JF_2 = Check_QEAL_ZERO(JF_2); JF_1 = Check_QEAL_ZERO(JF_1); JF_0 = Check_QEAL_ZERO(JF_0);

#ifdef CHECK_MEDIAL_PRIMITIVE_VALID
		qeal Aroots[2];
		int ArootCount = 0;
		assert(solveQuadricEquation(JA_2, JA_1, JA_0, Aroots, ArootCount, true, 0.0, 1.0) == false);
		qeal Croots[2];
		int CrootCount = 0;
		assert(solveQuadricEquation(JC_2, JC_1, JC_0, Croots, CrootCount, true, 0.0, 1.0) == false);
#endif 
		bool collide = false;
		int collideType = NO_COLLIDE;
		qeal pftc = 1.0; // potential ftc

		// dcd when t = 1.0;
		if (dcd(C1 + V1, C2 + V2, C3 + V3, R1, R2, R3, space_triangle, norrow))
			collideType = INTERSECTION; // check all cases;


		//
		qeal ftc_endpoints = ftc;
		bool collide_endpoints = checkEndpointAlphaBetaCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_endpoints, space_triangle);
		if (collide_endpoints && ftc_endpoints < ftc)
		{
			collide = true;
			ftc = ftc_endpoints;
		}

		qeal ftc_a0 = ftc;
		bool collide_a0 = checkAlphaIsZeroCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_a0);
		if (collide_a0 && ftc_a0 < ftc)
		{
			collide = true;
			ftc = ftc_a0;
		}

		qeal ftc_b0 = ftc;
		bool collide_b0 = checkBetaIsZeroCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_b0);
		if (collide_b0 && ftc_b0 < ftc)
		{
			collide = true;
			ftc = ftc_b0;
		}

		if (space_triangle)
		{
			qeal ftc_ab1 = ftc;
			bool collide_ab1 = checkAlphaPlusBetaIsOneCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_ab1);
			if (collide_ab1 && ftc_ab1 < ftc)
			{
				collide = true;
				ftc = ftc_ab1;
			}
		}
		else
		{
			qeal ftc_a1 = ftc;
			bool collide_a1 = checkAlphaIsOneCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_a1);
			if (collide_a1 && ftc_a1 < ftc)
			{
				collide = true;
				ftc = ftc_a1;
			}

			qeal ftc_b1 = ftc;
			bool collide_b1 = checkBetaIsOneCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_b1);
			if (collide_b1 && ftc_b1 < ftc)
			{
				collide = true;
				ftc = ftc_b1;
			}
		}

		qeal ftc_ab = ftc;
		bool collide_ab = checkAlphaBetaCCD2(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc_ab, collideType, space_triangle, norrow);
		if (collide_ab && ftc_ab < ftc)
		{
			collide = true;
			ftc = ftc_ab;
		}

		//if (collide && (ftc >= 1.0 || ftc <= 0.0))
		//{
		//	qeal t = 1e-10;
		//	while (t <= 1.0)
		//	{
		//		qeal At = JA_2 * t * t + JA_1 * t + JA_0;
		//		qeal Bt = JB_2 * t * t + JB_1 * t + JB_0;
		//		qeal Ct = JC_2 * t * t + JC_1 * t + JC_0;
		//		qeal Dt = JD_2 * t * t + JD_1 * t + JD_0;
		//		qeal Et = JE_2 * t * t + JE_1 * t + JE_0;
		//		qeal Ft = JF_2 * t * t + JF_1 * t + JF_0;
		//		qeal delta = 4.0 * At * Ct - Bt * Bt;
		//		qeal alpha = (Bt * Et - 2.0 * Ct * Dt) / delta;
		//		qeal beta = (Bt * Dt - 2.0 * At * Et) / delta;

		//		if (alpha < 0.0)alpha = 0.0;
		//		if (alpha > 1.0)alpha = 1.0;

		//		if (beta < 0.0)beta = 0.0;
		//		if (beta > 1.0)beta = 1.0;

		//		if (space_triangle && alpha + beta > 1.0)
		//		{
		//			alpha = -(Bt + Dt - 2.0 * Ct - Et) / (2.0 * (At - Bt + Ct));
		//			if (alpha < 0.0)alpha = 0.0;
		//			if (alpha > 1.0)alpha = 1.0;
		//			beta = 1.0 - alpha;
		//		}

		//		qeal fx = At * alpha * alpha + Bt * alpha * beta + Ct * beta * beta + Dt * alpha + Et * beta + Ft;
		//		t += 0.001;
		//		if (fx < 0)
		//			break;
		//	}
		//	ftc = t - 0.001;
		//	std::cout << "while " << std::endl;
		//}

		return collide;
	}

	bool checkEndpointAlphaBetaCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, bool space_triangle, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal J4_2, J4_1, J4_0;
		bool collide = false;
		// set alpha = 0.0, beta = 0.0
		//F(t) <= 0
		RootRange pos_00, neg_00;
		qeal roots_00[2];
		int rootCount_00 = 0;
		J1_2 = JF_2; J1_1 = JF_1; J1_0 = JF_0;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_00, rootCount_00, pos_00, neg_00);
		qeal ft00 = -1;
		if (rootCount_00 > 0)
			ft00 = roots_00[0];
		else if (!neg_00.isNull())
			ft00 = neg_00.set.front()[0];
		if (ft00 >= 0.0 && ft00 <= 1.0)
		{
			collide = true;
			if (ft00 < ftc) ftc = ft00;
		}

		// set alpha = 0.0, beta = 1.0
		//C(t) + E(t) + F(t) <= 0
		RootRange pos_01, neg_01;
		qeal roots_01[2];
		int rootCount_01 = 0;
		J1_2 = JC_2 + JE_2 + JF_2;
		J1_1 = JC_1 + JE_1 + JF_1;
		J1_0 = JC_0 + JE_0 + JF_0;

		solveQuadricNEq(J1_2, J1_1, J1_0, roots_01, rootCount_01, pos_01, neg_01);
		qeal ft01 = -1;
		if (rootCount_01 > 0)
			ft01 = roots_01[0];
		else if (!neg_01.isNull())
			ft01 = neg_01.set.front()[0];
		if (ft01 >= 0.0 && ft01 <= 1.0)
		{
			collide = true;
			if (ft01 < ftc) ftc = ft01;
		}

		// set alpha = 1.0, beta = 0.0
		//A(t) + D(t) + F(t) <= 0
		RootRange pos_10, neg_10;
		qeal roots_10[2];
		int rootCount_10 = 0;
		J1_2 = JA_2 + JD_2 + JF_2;
		J1_1 = JA_1 + JD_1 + JF_1;
		J1_0 = JA_0 + JD_0 + JF_0;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_10, rootCount_10, pos_10, neg_10);

		qeal ft10 = -1;
		if (rootCount_10 > 0)
			ft10 = roots_10[0];
		else if (!neg_10.isNull())
			ft10 = neg_10.set.front()[0];
		if (ft10 >= 0.0 && ft10 <= 1.0)
		{
			collide = true;
			if (ft10 < ftc) ftc = ft10;
		}

		if (!space_triangle)
		{
			// set alpha = 1.0, beta = 1.0
			//A(t) + B(t) + C(t) + D(t) + E(t) + F(t) <= 0
			RootRange pos_11, neg_11;
			qeal roots_11[2];
			int rootCount_11 = 0;
			J1_2 = JA_2 + JB_2 + JC_2 + JD_2 + JE_2 + JF_2;
			J1_1 = JA_1 + JB_1 + JC_1 + JD_1 + JE_1 + JF_1;
			J1_0 = JA_0 + JB_0 + JC_0 + JD_0 + JE_0 + JF_0;
			solveQuadricNEq(J1_2, J1_1, J1_0, roots_11, rootCount_11, pos_11, neg_11);
			qeal ft11 = -1;
			if (rootCount_11 > 0)
				ft11 = roots_11[0];
			else if (!neg_11.isNull())
				ft11 = neg_11.set.front()[0];
			if (ft11 >= 0.0 && ft11 <= 1.0)
			{
				collide = true;
				if (ft11 < ftc) ftc = ft11;
			}
		}
		return collide;
	}

	bool checkAlphaIsZeroCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;

		RootRange pos_0b, neg_0b;
		qeal roots_0b[6];
		int rootCount_0b = 0;
		// -E > 0
		J1_2 = -JE_2;
		J1_1 = -JE_1;
		J1_0 = -JE_0;
		if (IS_QEAL_ZERO(J1_2) && IS_QEAL_ZERO(J1_1) && IS_QEAL_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_0b, rootCount_0b, pos_0b, neg_0b);
		if (pos_0b.isNull())
			return false;

		//2C + E > 0
		J1_2 = 2.0 * JC_2 + JE_2;
		J1_1 = 2.0 * JC_1 + JE_1;
		J1_0 = 2.0 * JC_0 + JE_0;
		if (IS_QEAL_ZERO(J1_2) && IS_QEAL_ZERO(J1_1) && IS_QEAL_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_0b, rootCount_0b, pos_0b, neg_0b);
		if (pos_0b.isNull())
			return false;

		//E^2 - 4CF > 0;
		J1_2 = JE_2;
		J1_1 = JE_1;
		J1_0 = JE_0;

		J2_2 = J1_2;
		J2_1 = J1_1;
		J2_0 = J1_0;

		J3_2 = 4.0 * JC_2;
		J3_1 = 4.0 * JC_1;
		J3_0 = 4.0 * JC_0;

		qeal W_4 = 0, W_3 = 0, W_2 = 0, W_1 = 0, W_0 = 0;
		genQuarticCoeffs(W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, JF_2, JF_1, JF_0);
		solveQuarticNEq(W_4, W_3, W_2, W_1, W_0, roots_0b, rootCount_0b, pos_0b, neg_0b);
		if (pos_0b.isNull())
			return false;

		qeal At, Bt, Ct, Dt, Et, Ft;
		qeal beta;
		qeal fx;
		while (!pos_0b.set.empty())
		{
			Vector2 range = pos_0b.set.front();
			pos_0b.set.pop();
			qeal ft = range.data()[0];
			if (ft >= 0.0 && ft < ftc)
			{
				Ct = JC_2 * ft * ft + JC_1 * ft + JC_0;
				Et = JE_2 * ft * ft + JE_1 * ft + JE_0;
				Ft = JF_2 * ft * ft + JF_1 * ft + JF_0;

				beta = -(Et) / (2.0 * Ct);
				if (beta <= 0.0 || beta >= 1.0) return false;

				fx = Ct * beta * beta + Et * beta + Ft;
				if (fx >= 0.0 && fx <= 200.0 * std::max(norrow, CDMM_NORROW))
					ftc = ft;
				else if (fx < 0.0)
				{
					if (fx > -norrow) 
						ftc = ft;  // _distance still excess 0.0
					else // rock back
					{
						qeal scale = 0.9;
						qeal xp = ft;
						while (fx < -norrow)
						{
							xp = ft * scale;
							Ct = JC_2 * xp * xp + JC_1 * xp + JC_0;
							Et = JE_2 * xp * xp + JE_1 * xp + JE_0;
							Ft = JF_2 * xp * xp + JF_1 * xp + JF_0;
							beta = -(Et) / (2.0 * Ct);
							if (beta > 1.0 || beta < 0.0)
								break;
							fx = Ct * beta * beta + Et * beta + Ft;
							scale -= 0.05;
						}
						ftc = xp;
					}
				}
			}
		}
		return true;
	}

	bool checkAlphaIsOneCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal J4_2, J4_1, J4_0;

		RootRange pos_1b, neg_1b;
		qeal roots_1b[6];
		int rootCount_1b = 0;

		// -(B+ E) > 0
		J1_2 = -(JB_2 + JE_2);
		J1_1 = -(JB_1 + JE_1);
		J1_0 = -(JB_0 + JE_0);
		if (IS_QEAL_ZERO(J1_2) && IS_QEAL_ZERO(J1_1) && IS_QEAL_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_1b, rootCount_1b, pos_1b, neg_1b);
		if (pos_1b.isNull())
			return false;

		// 2C + B + E > 0
		J1_2 = 2.0 * JC_2 + JB_2 + JE_2;
		J1_1 = 2.0 * JC_1 + JB_1 + JE_1;
		J1_0 = 2.0 * JC_0 + JB_0 + JE_0;
		if (IS_QEAL_ZERO(J1_2) && IS_QEAL_ZERO(J1_1) && IS_QEAL_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_1b, rootCount_1b, pos_1b, neg_1b);
		if (pos_1b.isNull())
			return false;

		J1_2 = JB_2 + JE_2;
		J1_1 = JB_1 + JE_1;
		J1_0 = JB_0 + JE_0;

		J2_2 = J1_2;
		J2_1 = J1_1;
		J2_0 = J1_0;

		J3_2 = 4.0 * JC_2;
		J3_1 = 4.0 * JC_1;
		J3_0 = 4.0 * JC_0;

		J4_2 = JA_2 + JD_2 + JF_2;
		J4_1 = JA_1 + JD_1 + JF_1;
		J4_0 = JA_0 + JD_0 + JF_0;
		qeal W_4 = 0, W_3 = 0, W_2 = 0, W_1 = 0, W_0 = 0;
		genQuarticCoeffs(W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, J4_2, J4_1, J4_0);
		solveQuarticNEq(W_4, W_3, W_2, W_1, W_0, roots_1b, rootCount_1b, pos_1b, neg_1b);
		if (pos_1b.isNull())
			return false;

		qeal At, Bt, Ct, Dt, Et, Ft;
		qeal beta;
		qeal fx;
		while (!pos_1b.set.empty())
		{
			Vector2 range = pos_1b.set.front();
			pos_1b.set.pop();
			qeal ft = range.data()[0];
			if (ft >= 0.0 && ft < ftc)
			{
				At = JA_2 * ft * ft + JA_1 * ft + JA_0;
				Bt = JB_2 * ft * ft + JB_1 * ft + JB_0;
				Ct = JC_2 * ft * ft + JC_1 * ft + JC_0;
				Dt = JD_2 * ft * ft + JD_1 * ft + JD_0;
				Et = JE_2 * ft * ft + JE_1 * ft + JE_0;
				Ft = JF_2 * ft * ft + JF_1 * ft + JF_0;
				beta = -(Bt + Et) / (2.0 * Ct);
				if (beta <= 0.0 || beta >= 1.0) return false;

				fx = Ct * beta * beta + (Bt + Et) * beta + (At + Dt + Ft);
	
				if (fx >= 0.0 && fx <= 200.0 * std::max(norrow, CDMM_NORROW))
					ftc = ft;
				else if (fx < 0.0)
				{
					if (fx > -norrow) ftc = ft;  // _distance still excess 0.0
					else // rock back
					{
						qeal scale = 0.9;
						qeal xp = ft;
						while (fx < -norrow)
						{
							xp = ft * scale;
							At = JA_2 * xp * xp + JA_1 * xp + JA_0;
							Bt = JB_2 * xp * xp + JB_1 * xp + JB_0;
							Ct = JC_2 * xp * xp + JC_1 * xp + JC_0;
							Dt = JD_2 * xp * xp + JD_1 * xp + JD_0;
							Et = JE_2 * xp * xp + JE_1 * xp + JE_0;
							Ft = JF_2 * xp * xp + JF_1 * xp + JF_0;
							beta = -(Bt + Et) / (2.0 * Ct);
							if (beta > 1.0 || beta < 0.0)
								break;
							fx = Ct * beta * beta + (Bt + Et) * beta + (At + Dt + Ft);
							scale -= 0.05;
						}
						ftc = xp;
					}
				}

			}
		}
		return true;
	}

	bool checkBetaIsZeroCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		RootRange pos_a0, neg_a0;
		qeal roots_a0[6];
		int rootCount_a0 = 0;		
		// -D > 0
		J1_2 = -JD_2;
		J1_1 = -JD_1;
		J1_0 = -JD_0;
		if (IS_QEAL_ZERO(J1_2) && IS_QEAL_ZERO(J1_1) && IS_QEAL_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_a0, rootCount_a0, pos_a0, neg_a0);
		if (pos_a0.isNull())
			return false;
		
		//2A + D > 0
		J1_2 = 2.0 * JA_2 + JD_2;
		J1_1 = 2.0 * JA_1 + JD_1;
		J1_0 = 2.0 * JA_0 + JD_0;
		if (IS_QEAL_ZERO(J1_2) && IS_QEAL_ZERO(J1_1) && IS_QEAL_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_a0, rootCount_a0, pos_a0, neg_a0);
		if (pos_a0.isNull())
			return false;

		//D^2 - 4AF > 0;
		J1_2 = JD_2;
		J1_1 = JD_1;
		J1_0 = JD_0;

		J2_2 = J1_2;
		J2_1 = J1_1;
		J2_0 = J1_0;

		J3_2 = 4.0 * JA_2;
		J3_1 = 4.0 * JA_1;
		J3_0 = 4.0 * JA_0;

		qeal W_4 = 0, W_3 = 0, W_2 = 0, W_1 = 0, W_0 = 0;
		genQuarticCoeffs(W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, JF_2, JF_1, JF_0);
		solveQuarticNEq(W_4, W_3, W_2, W_1, W_0, roots_a0, rootCount_a0, pos_a0, neg_a0);
		if (pos_a0.isNull())
			return false;

		qeal At, Bt, Ct, Dt, Et, Ft;
		qeal alpha;
		qeal fx;
		while (!pos_a0.set.empty())
		{
			Vector2 range = pos_a0.set.front();
			pos_a0.set.pop();
			qeal ft = range.data()[0];
			if (ft > 0.0 && ft < ftc)
			{
				ftc = range.data()[0];
				At = JA_2 * ft * ft + JA_1 * ft + JA_0;
				Dt = JD_2 * ft * ft + JD_1 * ft + JD_0;
				Ft = JF_2 * ft * ft + JF_1 * ft + JF_0;
				alpha = -(Dt) / (2.0 * At);
				if (alpha <= 0.0 || alpha >= 1.0) return false;
				fx = At * alpha * alpha + Dt * alpha + Ft;

				if (fx >= 0.0 && fx <= 200.0 * std::max(norrow, CDMM_NORROW))
					ftc = ft;
				else if (fx < 0.0)
				{
					if (fx > -norrow) ftc = ft;  // _distance still excess 0.0
					else // rock back
					{
						qeal scale = 0.9;
						qeal xp = ft;
						while (fx < -norrow)
						{
							xp = ft * scale;
							At = JA_2 * xp * xp + JA_1 * xp + JA_0;
							Dt = JD_2 * xp * xp + JD_1 * xp + JD_0;
							Ft = JF_2 * xp * xp + JF_1 * xp + JF_0;
							alpha = -(Dt) / (2.0 * At);
							if (alpha > 1.0 || alpha < 0.0)
								break;
							fx = At * alpha * alpha + Dt * alpha + Ft;
							scale -= 0.05;
						}
						ftc = xp;
					}
				}

			}
		}
		return true;
	}

	bool checkBetaIsOneCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal J4_2, J4_1, J4_0;

		RootRange pos_a1, neg_a1;
		qeal roots_a1[6];
		int rootCount_a1 = 0;
		
		// -(B+ D) > 0
		J1_2 = -(JB_2 + JD_2);
		J1_1 = -(JB_1 + JD_1);
		J1_0 = -(JB_0 + JD_0);
		if (IS_QEAL_ZERO(J1_2) && IS_QEAL_ZERO(J1_1) && IS_QEAL_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_a1, rootCount_a1, pos_a1, neg_a1);
		if (pos_a1.isNull())
			return false;

		// 2A + B + D > 0
		J1_2 = 2.0 * JA_2 + JB_2 + JD_2;
		J1_1 = 2.0 * JA_1 + JB_1 + JD_1;
		J1_0 = 2.0 * JA_0 + JB_0 + JD_0;
		if (IS_QEAL_ZERO(J1_2) && IS_QEAL_ZERO(J1_1) && IS_QEAL_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_a1, rootCount_a1, pos_a1, neg_a1);
		if (pos_a1.isNull())
			return false;
		// (B+ D)^2 - 4.0 * A * (C + E + F) > 0
		J1_2 = JB_2 + JD_2;
		J1_1 = JB_1 + JD_1;
		J1_0 = JB_0 + JD_0;

		J2_2 = J1_2;
		J2_1 = J1_1;
		J2_0 = J1_0;

		J3_2 = 4.0 * JA_2;
		J3_1 = 4.0 * JA_1;
		J3_0 = 4.0 * JA_0;

		J4_2 = JC_2 + JE_2 + JF_2;
		J4_1 = JC_1 + JE_1 + JF_1;
		J4_0 = JC_0 + JE_0 + JF_0;
		qeal W_4 = 0, W_3 = 0, W_2 = 0, W_1 = 0, W_0 = 0;
		genQuarticCoeffs(W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, J4_2, J4_1, J4_0);
		solveQuarticNEq(W_4, W_3, W_2, W_1, W_0, roots_a1, rootCount_a1, pos_a1, neg_a1);
		if (pos_a1.isNull())
			return false;

		qeal At, Bt, Ct, Dt, Et, Ft;
		qeal alpha;
		qeal fx;
		while (!pos_a1.set.empty())
		{
			Vector2 range = pos_a1.set.front();
			pos_a1.set.pop();
			qeal ft = range.data()[0];
			if (ft > 0.0 && ft < ftc)
			{
				At = JA_2 * ft * ft + JA_1 * ft + JA_0;
				Bt = JB_2 * ft * ft + JB_1 * ft + JB_0;
				Ct = JC_2 * ft * ft + JC_1 * ft + JC_0;
				Dt = JD_2 * ft * ft + JD_1 * ft + JD_0;
				Et = JE_2 * ft * ft + JE_1 * ft + JE_0;
				Ft = JF_2 * ft * ft + JF_1 * ft + JF_0;
				alpha = -(Bt + Dt) / (2.0 * At);
				if (alpha <= 0.0 || alpha >= 1.0) return false;
				fx = At * alpha * alpha + (Bt + Dt) * alpha + (Ct + Et + Ft);

				if (fx >= 0.0 && fx <= 200.0 * std::max(norrow, CDMM_NORROW))
					ftc = ft;
				else if (fx < 0.0)
				{
					if (fx > -norrow) ftc = ft;  // _distance still excess 0.0
					else // rock back
					{
						qeal scale = 0.9;
						qeal xp = ft;
						while (fx < -norrow)
						{
							xp = ft * scale;
							At = JA_2 * xp * xp + JA_1 * xp + JA_0;
							Bt = JB_2 * xp * xp + JB_1 * xp + JB_0;
							Ct = JC_2 * xp * xp + JC_1 * xp + JC_0;
							Dt = JD_2 * xp * xp + JD_1 * xp + JD_0;
							Et = JE_2 * xp * xp + JE_1 * xp + JE_0;
							Ft = JF_2 * xp * xp + JF_1 * xp + JF_0;
							alpha = -(Bt + Dt) / (2.0 * At);
							if (alpha > 1.0 || alpha < 0.0)
								break;
							fx = At * alpha * alpha + (Bt + Dt) * alpha + (Ct + Et + Ft);
							scale -= 0.05;
						}
						ftc = xp;
					}
				}

			}
		}
		return true;
	}

	bool checkAlphaPlusBetaIsOneCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal J4_2, J4_1, J4_0;

		//
		RootRange pos_ab1, neg_ab1;
		qeal roots_ab1[6];
		int rootCount_ab1 = 0;

		// -(B+ D - 2C - E) > 0
		J1_2 = -(JB_2 + JD_2 - 2.0 * JC_2 - JE_2);
		J1_1 = -(JB_1 + JD_1 - 2.0 * JC_1 - JE_1);
		J1_0 = -(JB_0 + JD_0 - 2.0 * JC_0 - JE_0);
		if (IS_QEAL_ZERO(J1_2) && IS_QEAL_ZERO(J1_1) && IS_QEAL_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_ab1, rootCount_ab1, pos_ab1, neg_ab1);
		if (pos_ab1.isNull())
			return false;

		// 2.0 *  (A + C - B) + (B + D - 2.0 * C - E) > 0
		J1_2 = (2.0 * JA_2 + JD_2) - (JB_2 + JE_2);
		J1_1 = (2.0 * JA_1 + JD_1) - (JB_1 + JE_1);
		J1_0 = (2.0* JA_0 + JD_0) - (JB_0 + JE_0);

		if (IS_QEAL_ZERO(J1_2) && IS_QEAL_ZERO(J1_1) && IS_QEAL_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_ab1, rootCount_ab1, pos_ab1, neg_ab1);
		if (pos_ab1.isNull())
			return false;

		// (B+ D - 2C - E)^2 - 4.0 * (A + C - B) * (C + E + F) > 0
		J1_2 = -(JB_2 + JD_2 - 2.0 * JC_2 - JE_2);
		J1_1 = -(JB_1 + JD_1 - 2.0 * JC_1 - JE_1);
		J1_0 = -(JB_0 + JD_0 - 2.0 * JC_0 - JE_0);
		J2_2 = J1_2;
		J2_1 = J1_1;
		J2_0 = J1_0;

		J3_2 = 4.0 * (JA_2 + JC_2 - JB_2);
		J3_1 = 4.0 * (JA_1 + JC_1 - JB_1);
		J3_0 = 4.0 * (JA_0 + JC_0 - JB_0);

		J4_2 = JC_2 + JE_2 + JF_2;
		J4_1 = JC_1 + JE_1 + JF_1;
		J4_0 = JC_0 + JE_0 + JF_0;

		qeal W_4 = 0, W_3 = 0, W_2 = 0, W_1 = 0, W_0 = 0;
		genQuarticCoeffs(W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, J4_2, J4_1, J4_0);
		solveQuarticNEq(W_4, W_3, W_2, W_1, W_0, roots_ab1, rootCount_ab1, pos_ab1, neg_ab1);
		if (pos_ab1.isNull())
			return false;

		qeal At, Bt, Ct, Dt, Et, Ft;
		qeal alpha, beta;
		qeal fx;
		while (!pos_ab1.set.empty())
		{
			Vector2 range = pos_ab1.set.front();
			pos_ab1.set.pop();
			qeal ft = range.data()[0];
			if (ft > 0.0 && ft < 1.0 && ft < ftc)
			{
				At = JA_2 * ft * ft + JA_1 * ft + JA_0;
				Bt = JB_2 * ft * ft + JB_1 * ft + JB_0;
				Ct = JC_2 * ft * ft + JC_1 * ft + JC_0;
				Dt = JD_2 * ft * ft + JD_1 * ft + JD_0;
				Et = JE_2 * ft * ft + JE_1 * ft + JE_0;
				Ft = JF_2 * ft * ft + JF_1 * ft + JF_0;

				alpha = -(Bt + Dt - 2.0 * Ct - Et) / (2.0 * (At - Bt + Ct));
				if (alpha <= 0.0 || alpha >= 1.0) return false;
				beta = 1.0 - alpha;
				fx = At * alpha * alpha + Bt * alpha * beta + Ct * beta  * beta + Dt * alpha + Et * beta + Ft;

				if (fx >= 0.0 && fx <= 200.0 * std::max(norrow, CDMM_NORROW))
					ftc = ft;
				else if (fx < 0.0)
				{
					if (fx > -norrow) ftc = ft;  // _distance still excess 0.0
					else // rock back
					{
						qeal scale = 0.9;
						qeal xp = ft;
						while (fx < -norrow)
						{
							xp = ft * scale;
							At = JA_2 * xp * xp + JA_1 * xp + JA_0;
							Bt = JB_2 * xp * xp + JB_1 * xp + JB_0;
							Ct = JC_2 * xp * xp + JC_1 * xp + JC_0;
							Dt = JD_2 * xp * xp + JD_1 * xp + JD_0;
							Et = JE_2 * xp * xp + JE_1 * xp + JE_0;
							Ft = JF_2 * xp * xp + JF_1 * xp + JF_0;
							alpha = -(Bt + Dt - 2.0 * Ct - Et) / (2.0 * (At - Bt + Ct));
							beta = 1.0 - alpha;
							if (alpha > 1.0 || alpha < 0.0 || beta > 1.0 || beta < 0.0)
								break;
							fx = At * alpha * alpha + Bt * alpha * beta + Ct * beta  * beta + Dt * alpha + Et * beta + Ft;
							scale -= 0.05;
						}
						ftc = xp;
					}
				}

			}
		}
		return true;
	}

	bool checkAlphaBetaCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, bool space_triangle, qeal norrow)
	{
		RootRange pos, neg;
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal J4_2, J4_1, J4_0;

		qeal L_4 = 0, L_3 = 0, L_2 = 0, L_1 = 0, L_0 = 0;
		qeal K_4 = 0, K_3 = 0, K_2 = 0, K_1 = 0, K_0 = 0;
		qeal H_4 = 0, H_3 = 0, H_2 = 0, H_1 = 0, H_0 = 0;
		qeal P_4 = 0, P_3 = 0, P_2 = 0, P_1 = 0, P_0 = 0;
		qeal Q_4 = 0, Q_3 = 0, Q_2 = 0, Q_1 = 0, Q_0 = 0;

		qeal roots_delta[6];
		int rootCount_delta = 0;
		// 4AC-B^2 != 0
		J1_2 = 4.0 * JA_2;
		J1_1 = 4.0 * JA_1;
		J1_0 = 4.0 * JA_0;

		J2_2 = JC_2;
		J2_1 = JC_1;
		J2_0 = JC_0;

		J3_2 = JB_2;
		J3_1 = JB_1;
		J3_0 = JB_0;

		J4_2 = JB_2;
		J4_1 = JB_1;
		J4_0 = JB_0;

		//4AC-B^2 = H4t^4 + H3t^3 + H2t^2 + H1t + H0;
		genQuarticCoeffs(H_4, H_3, H_2, H_1, H_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, J4_2, J4_1, J4_0);
		solveQuarticNEq(H_4, H_3, H_2, H_1, H_0, roots_delta, rootCount_delta, pos, neg);
		//
		bool collide = false;
		for (int i = 0; i < rootCount_delta; i++)
		{
			qeal ft = roots_delta[i];
			qeal At = JA_2 * ft * ft + JA_1 * ft + JA_0;
			qeal Bt = JB_2 * ft * ft + JB_1 * ft + JB_0;
			qeal Ct = JC_2 * ft * ft + JC_1 * ft + JC_0;
			qeal Dt = JD_2 * ft * ft + JD_1 * ft + JD_0;
			qeal Et = JE_2 * ft * ft + JE_1 * ft + JE_0;
			qeal Ft = JF_2 * ft * ft + JF_1 * ft + JF_0;
			if (dcd(At, Bt, Ct, Dt, Et, Ft, false) && ft < ftc)
			{
				collide = true;
				ftc = ft;
			}
		}

		//
		// BE - 2CD
		qeal roots_ga0[6];
		int rootCount_ga0 = 0;
		J1_2 = JB_2;
		J1_1 = JB_1;
		J1_0 = JB_0;

		J2_2 = JE_2;
		J2_1 = JE_1;
		J2_0 = JE_0;

		J3_2 = 2.0 * JC_2;
		J3_1 = 2.0 * JC_1;
		J3_0 = 2.0 * JC_0;

		J4_2 = JD_2;
		J4_1 = JD_1;
		J4_0 = JD_0;
		genQuarticCoeffs(L_4, L_3, L_2, L_1, L_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, J4_2, J4_1, J4_0);
		solveQuarticNEq(L_4, L_3, L_2, L_1, L_0, roots_ga0, rootCount_ga0, pos, neg);
		if (pos.isNull())
			return collide;

		// BD - 2AE
		qeal roots_gb0[6];
		int rootCount_gb0 = 0;
		J2_2 = JD_2;
		J2_1 = JD_1;
		J2_0 = JD_0;

		J3_2 = 2.0 * JA_2;
		J3_1 = 2.0 * JA_1;
		J3_0 = 2.0 * JA_0;

		J4_2 = JE_2;
		J4_1 = JE_1;
		J4_0 = JE_0;
		genQuarticCoeffs(K_4, K_3, K_2, K_1, K_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, J4_2, J4_1, J4_0);
		solveQuarticNEq(K_4, K_3, K_2, K_1, K_0, roots_gb0, rootCount_gb0, pos, neg);
		if (pos.isNull())
			return collide;

		if (!space_triangle)
		{
			// (4AC - B^2) -  (BE - 2CD)
			qeal roots_ga1[6];
			int rootCount_ga1 = 0;
			P_4 = H_4 - L_4;
			P_3 = H_3 - L_3;
			P_2 = H_2 - L_2;
			P_1 = H_1 - L_1;
			P_0 = H_0 - L_0;
			solveQuarticNEq(P_4, P_3, P_2, P_1, P_0, roots_ga1, rootCount_ga1, pos, neg);
			if (pos.isNull())
				return collide;

			// (4AC - B^2) -  (BD - 2AE)
			qeal roots_gb1[6];
			int rootCount_gb1 = 0;
			Q_4 = H_4 - K_4;
			Q_3 = H_3 - K_3;
			Q_2 = H_2 - K_2;
			Q_1 = H_1 - K_1;
			Q_0 = H_0 - K_0;
			solveQuarticNEq(Q_4, Q_3, Q_2, Q_1, Q_0, roots_gb1, rootCount_gb1, pos, neg);
			if (pos.isNull())
				return collide;
		}
		else
		{
			qeal roots_gab1[6];
			int rootCount_gab1 = 0;
			P_4 = H_4 - L_4 - K_4;
			P_3 = H_3 - L_3 - K_3;
			P_2 = H_2 - L_2 - K_2;
			P_1 = H_1 - L_1 - K_1;
			P_0 = H_0 - L_0 - K_0;
			solveQuarticNEq(P_4, P_3, P_2, P_1, P_0, roots_gab1, rootCount_gab1, pos, neg);
			if (pos.isNull())
				return collide;
		}

		collide = checkAlphaBetaCCD(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc, rootCount_delta, roots_delta, pos, neg, space_triangle);
		return collide;
	}

	bool checkAlphaBetaCCD2(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, int collideType, bool space_triangle, qeal norrow)
	{
		RootRange pos, neg;
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal J4_2, J4_1, J4_0;

		qeal L_4 = 0, L_3 = 0, L_2 = 0, L_1 = 0, L_0 = 0;
		qeal K_4 = 0, K_3 = 0, K_2 = 0, K_1 = 0, K_0 = 0;
		qeal H_4 = 0, H_3 = 0, H_2 = 0, H_1 = 0, H_0 = 0;
		qeal P_4 = 0, P_3 = 0, P_2 = 0, P_1 = 0, P_0 = 0;
		qeal Q_4 = 0, Q_3 = 0, Q_2 = 0, Q_1 = 0, Q_0 = 0;

		qeal roots_delta[6];
		int rootCount_delta = 0;
		// 4AC-B^2 != 0
		J1_2 = 4.0 * JA_2;
		J1_1 = 4.0 * JA_1;
		J1_0 = 4.0 * JA_0;

		J2_2 = JC_2;
		J2_1 = JC_1;
		J2_0 = JC_0;

		J3_2 = JB_2;
		J3_1 = JB_1;
		J3_0 = JB_0;

		J4_2 = JB_2;
		J4_1 = JB_1;
		J4_0 = JB_0;

		//4AC-B^2 = H4t^4 + H3t^3 + H2t^2 + H1t + H0;
		genQuarticCoeffs(H_4, H_3, H_2, H_1, H_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, J4_2, J4_1, J4_0);
		solveQuarticNEq(H_4, H_3, H_2, H_1, H_0, roots_delta, rootCount_delta, pos, neg);
		//
		bool collide = false;
		for (int i = 0; i < rootCount_delta; i++)
		{
			qeal ft = roots_delta[i];
			qeal At = JA_2 * ft * ft + JA_1 * ft + JA_0;
			qeal Bt = JB_2 * ft * ft + JB_1 * ft + JB_0;
			qeal Ct = JC_2 * ft * ft + JC_1 * ft + JC_0;
			qeal Dt = JD_2 * ft * ft + JD_1 * ft + JD_0;
			qeal Et = JE_2 * ft * ft + JE_1 * ft + JE_0;
			qeal Ft = JF_2 * ft * ft + JF_1 * ft + JF_0;
			if (dcd(At, Bt, Ct, Dt, Et, Ft, false) && ft < ftc)
			{
				collide = true;
				ftc = ft;
			}
		}

		//
		// BE - 2CD
		qeal roots_ga0[6];
		int rootCount_ga0 = 0;
		J1_2 = JB_2;
		J1_1 = JB_1;
		J1_0 = JB_0;

		J2_2 = JE_2;
		J2_1 = JE_1;
		J2_0 = JE_0;

		J3_2 = 2.0 * JC_2;
		J3_1 = 2.0 * JC_1;
		J3_0 = 2.0 * JC_0;

		J4_2 = JD_2;
		J4_1 = JD_1;
		J4_0 = JD_0;
		genQuarticCoeffs(L_4, L_3, L_2, L_1, L_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, J4_2, J4_1, J4_0);
		solveQuarticNEq(L_4, L_3, L_2, L_1, L_0, roots_ga0, rootCount_ga0, pos, neg);
		if (pos.isNull())
			return collide;

		// BD - 2AE
		qeal roots_gb0[6];
		int rootCount_gb0 = 0;
		J2_2 = JD_2;
		J2_1 = JD_1;
		J2_0 = JD_0;

		J3_2 = 2.0 * JA_2;
		J3_1 = 2.0 * JA_1;
		J3_0 = 2.0 * JA_0;

		J4_2 = JE_2;
		J4_1 = JE_1;
		J4_0 = JE_0;
		genQuarticCoeffs(K_4, K_3, K_2, K_1, K_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, J4_2, J4_1, J4_0);
		solveQuarticNEq(K_4, K_3, K_2, K_1, K_0, roots_gb0, rootCount_gb0, pos, neg);
		if (pos.isNull())
			return collide;

		if (!space_triangle)
		{
			// (4AC - B^2) -  (BE - 2CD)
			qeal roots_ga1[6];
			int rootCount_ga1 = 0;
			P_4 = H_4 - L_4;
			P_3 = H_3 - L_3;
			P_2 = H_2 - L_2;
			P_1 = H_1 - L_1;
			P_0 = H_0 - L_0;
			solveQuarticNEq(P_4, P_3, P_2, P_1, P_0, roots_ga1, rootCount_ga1, pos, neg);
			if (pos.isNull())
				return collide;

			// (4AC - B^2) -  (BD - 2AE)
			qeal roots_gb1[6];
			int rootCount_gb1 = 0;
			Q_4 = H_4 - K_4;
			Q_3 = H_3 - K_3;
			Q_2 = H_2 - K_2;
			Q_1 = H_1 - K_1;
			Q_0 = H_0 - K_0;
			solveQuarticNEq(Q_4, Q_3, Q_2, Q_1, Q_0, roots_gb1, rootCount_gb1, pos, neg);
			if (pos.isNull())
				return collide;
		}
		else
		{
			qeal roots_gab1[6];
			int rootCount_gab1 = 0;
			P_4 = H_4 - L_4 - K_4;
			P_3 = H_3 - L_3 - K_3;
			P_2 = H_2 - L_2 - K_2;
			P_1 = H_1 - L_1 - K_1;
			P_0 = H_0 - L_0 - K_0;
			solveQuarticNEq(P_4, P_3, P_2, P_1, P_0, roots_gab1, rootCount_gab1, pos, neg);
			if (pos.isNull())
				return collide;
		}

		collide = checkAlphaBetaCCD2(JA_2, JA_1, JA_0, JB_2, JB_1, JB_0, JC_2, JC_1, JC_0, JD_2, JD_1, JD_0, JE_2, JE_1, JE_0, JF_2, JF_1, JF_0, ftc, rootCount_delta, roots_delta, pos, neg, collideType, space_triangle);
		return collide;
	}

	bool checkAlphaBetaCCD(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, int banRootsNum, qeal* banRoots, RootRange & pos, RootRange & neg, bool space_triangle, qeal norrow)
	{
		bool collide = false;
		if (pos.isNull() && neg.isNull()) return false;

		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;

		qeal ft = 1.0;
		qeal W_6 = 0, W_5 = 0, W_4 = 0, W_3 = 0, W_2 = 0, W_1 = 0, W_0 = 0;
		//4A(t)C(t)F(t) + B(t)D(t)E(t) - A(t)E(t)^2 - C(t)D(t)^2 - B(t)^2F(t) = 0
		J1_2 = JA_2; 		J1_1 = JA_1; 		J1_0 = JA_0;
		J2_2 = JC_2; 		J2_1 = JC_1; 		J2_0 = JC_0;
		J3_2 = JF_2; 		J3_1 = JF_1; 		J3_0 = JF_0;
		genSexticCoeffs(W_6, W_5, W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, 4.0);

		J1_2 = JB_2; 		J1_1 = JB_1; 		J1_0 = JB_0;
		J2_2 = JD_2; 		J2_1 = JD_1; 		J2_0 = JD_0;
		J3_2 = JE_2; 		J3_1 = JE_1; 		J3_0 = JE_0;
		genSexticCoeffs(W_6, W_5, W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, 1.0);

		J1_2 = JA_2; 		J1_1 = JA_1; 		J1_0 = JA_0;
		J2_2 = JE_2; 		J2_1 = JE_1; 		J2_0 = JE_0;
		genSexticCoeffs(W_6, W_5, W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, -1.0);

		J1_2 = JC_2; 		J1_1 = JC_1; 		J1_0 = JC_0;
		J2_2 = JD_2; 		J2_1 = JD_1; 		J2_0 = JD_0;
		J3_2 = JD_2; 		J3_1 = JD_1; 		J3_0 = JD_0;
		genSexticCoeffs(W_6, W_5, W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, -1.0);

		J1_2 = JB_2; 		J1_1 = JB_1; 		J1_0 = JB_0;
		J2_2 = JB_2; 		J2_1 = JB_1; 		J2_0 = JB_0;
		J3_2 = JF_2; 		J3_1 = JF_1; 		J3_0 = JF_0;
		genSexticCoeffs(W_6, W_5, W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, -1.0);

		// W_6 t^6 + W_5 t^5 + W_4 t^4 + W_3 t^3 + W_2 t^2 + W_1 t + W_0 = 0
		//W_6 = Check_QEAL_ZERO(W_6);
		//W_5 = Check_QEAL_ZERO(W_5);
		//W_4 = Check_QEAL_ZERO(W_4);
		//W_3 = Check_QEAL_ZERO(W_3);
		//W_2 = Check_QEAL_ZERO(W_2);
		//W_1 = Check_QEAL_ZERO(W_1);
		//W_0 = Check_QEAL_ZERO(W_0);

		if (!pos.isNull())
			ft = SolveSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, banRootsNum, banRoots, pos, norrow);
		if (!neg.isNull())
		{
			qeal temp = SolveSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, banRootsNum, banRoots, neg, norrow);
			if (ft > temp) ft = temp;
		}

		//double check
		qeal At = JA_2 * ft * ft + JA_1 * ft + JA_0;
		qeal Bt = JB_2 * ft * ft + JB_1 * ft + JB_0;
		qeal Ct = JC_2 * ft * ft + JC_1 * ft + JC_0;
		qeal Dt = JD_2 * ft * ft + JD_1 * ft + JD_0;
		qeal Et = JE_2 * ft * ft + JE_1 * ft + JE_0;
		qeal Ft = JF_2 * ft * ft + JF_1 * ft + JF_0;

		qeal delta = 4.0 * At * Ct - Bt * Bt;
		if (!IS_QEAL_ZERO(delta))
		{
			qeal alpha = (Bt * Et - 2.0 * Ct * Dt) / delta;
			qeal beta = (Bt * Dt - 2.0 * At * Et) / delta;
			bool valid = false;

			if(!space_triangle)
			{
				if (alpha > 0.0 && alpha < 1.0 && beta > 0.0 && beta < 1.0) valid = true;
			}
			else
			{
				if (alpha > 0.0 && alpha < 1.0 && beta > 0.0 && beta < 1.0 && (alpha + beta)) valid = true;
			}
			if (valid)
			{
				qeal fx = At * alpha * alpha + Bt * alpha * beta + Ct * beta * beta + Dt * alpha + Et * beta + Ft;

				if (ft > 0.0 && ft < 1.0)
				{
					collide = true;
				}
				else if (ft == 1.0 && IS_QEAL_ZERO(fx))
				{
					collide = true;
					ft = 0.95;
				}
			}
		}
		else
		{			
			if (dcd(At, Bt, Ct, Dt, Et, Ft, false)) // ensure non-intersection
			{
				collide = true;
			}
		}
		if (collide && ft < ftc) ftc = ft;
		return collide;
	}

	bool checkAlphaBetaCCD2(qeal & JA_2, qeal & JA_1, qeal & JA_0, qeal & JB_2, qeal & JB_1, qeal & JB_0, qeal & JC_2, qeal & JC_1, qeal & JC_0, qeal & JD_2, qeal & JD_1, qeal & JD_0, qeal & JE_2, qeal & JE_1, qeal & JE_0, qeal & JF_2, qeal & JF_1, qeal & JF_0, qeal & ftc, int banRootsNum, qeal * banRoots, RootRange & pos, RootRange & neg, int collideType, bool space_triangle, qeal norrow)
	{
		bool collide = false;
		if (pos.isNull() && neg.isNull()) return false;

		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;

		qeal ft = 1.0;
		qeal W_6 = 0, W_5 = 0, W_4 = 0, W_3 = 0, W_2 = 0, W_1 = 0, W_0 = 0;
		//4A(t)C(t)F(t) + B(t)D(t)E(t) - A(t)E(t)^2 - C(t)D(t)^2 - B(t)^2F(t) = 0
		J1_2 = JA_2; 		J1_1 = JA_1; 		J1_0 = JA_0;
		J2_2 = JC_2; 		J2_1 = JC_1; 		J2_0 = JC_0;
		J3_2 = JF_2; 		J3_1 = JF_1; 		J3_0 = JF_0;
		genSexticCoeffs(W_6, W_5, W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, 4.0);

		J1_2 = JB_2; 		J1_1 = JB_1; 		J1_0 = JB_0;
		J2_2 = JD_2; 		J2_1 = JD_1; 		J2_0 = JD_0;
		J3_2 = JE_2; 		J3_1 = JE_1; 		J3_0 = JE_0;
		genSexticCoeffs(W_6, W_5, W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, 1.0);

		J1_2 = JA_2; 		J1_1 = JA_1; 		J1_0 = JA_0;
		J2_2 = JE_2; 		J2_1 = JE_1; 		J2_0 = JE_0;
		genSexticCoeffs(W_6, W_5, W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, -1.0);

		J1_2 = JC_2; 		J1_1 = JC_1; 		J1_0 = JC_0;
		J2_2 = JD_2; 		J2_1 = JD_1; 		J2_0 = JD_0;
		J3_2 = JD_2; 		J3_1 = JD_1; 		J3_0 = JD_0;
		genSexticCoeffs(W_6, W_5, W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, -1.0);

		J1_2 = JB_2; 		J1_1 = JB_1; 		J1_0 = JB_0;
		J2_2 = JB_2; 		J2_1 = JB_1; 		J2_0 = JB_0;
		J3_2 = JF_2; 		J3_1 = JF_1; 		J3_0 = JF_0;
		genSexticCoeffs(W_6, W_5, W_4, W_3, W_2, W_1, W_0, J1_2, J1_1, J1_0, J2_2, J2_1, J2_0, J3_2, J3_1, J3_0, -1.0);

		W_6 = Check_QEAL_ZERO(W_6);
		W_5 = Check_QEAL_ZERO(W_5);
		W_4 = Check_QEAL_ZERO(W_4);
		W_3 = Check_QEAL_ZERO(W_3);
		W_2 = Check_QEAL_ZERO(W_2);
		W_1 = Check_QEAL_ZERO(W_1);
		W_0 = Check_QEAL_ZERO(W_0);
		//if (collideType == NO_COLLIDE)
		//{
		//	qeal fd0 = W_1;
		//	qeal fd1 = 6 * W_6 + 5 * W_5 + 4 * W_4 + 3 * W_3 + 2 * W_2 + W_1;

 	//		if (fd0 * fd1 > 0) return false;
		//}

		if (!pos.isNull())
			ft = SolveSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, banRootsNum, banRoots, pos, norrow);
		if (!neg.isNull())
		{
			qeal temp = SolveSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, banRootsNum, banRoots, neg, norrow);
			if (ft > temp) ft = temp;
		}

		//double check
		qeal At = JA_2 * ft * ft + JA_1 * ft + JA_0;
		qeal Bt = JB_2 * ft * ft + JB_1 * ft + JB_0;
		qeal Ct = JC_2 * ft * ft + JC_1 * ft + JC_0;
		qeal Dt = JD_2 * ft * ft + JD_1 * ft + JD_0;
		qeal Et = JE_2 * ft * ft + JE_1 * ft + JE_0;
		qeal Ft = JF_2 * ft * ft + JF_1 * ft + JF_0;

		qeal delta = 4.0 * At * Ct - Bt * Bt;
		if (!IS_QEAL_ZERO(delta))
		{
			qeal alpha = (Bt * Et - 2.0 * Ct * Dt) / delta;
			qeal beta = (Bt * Dt - 2.0 * At * Et) / delta;
			bool valid = false;

			if (!space_triangle)
			{
				if (alpha > 0.0 && alpha < 1.0 && beta > 0.0 && beta < 1.0) valid = true;
			}
			else
			{
				if (alpha > 0.0 && alpha < 1.0 && beta > 0.0 && beta < 1.0 && (alpha + beta)) valid = true;
			}
			if (valid)
			{
				qeal fx = At * alpha * alpha + Bt * alpha * beta + Ct * beta * beta + Dt * alpha + Et * beta + Ft;

				if (ft > 0.0 && ft < 1.0 && fx <= 200.0 * std::max(norrow, CDMM_NORROW))
				{
					collide = true;
				}
				else if (ft == 1.0 && IS_QEAL_ZERO(fx))
				{
					collide = true;
					ft = 0.95;
				}
			}
		}
		else
		{
			if (dcd(At, Bt, Ct, Dt, Et, Ft, false)) // ensure non-intersection
			{
				collide = true;
			}
		}
		if (collide && ft < ftc) ftc = ft;
		return collide;
	}

}