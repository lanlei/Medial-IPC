#include "Simulator\CollisionDetection\CollisionDetectionMedialMesh.h"

namespace CDMM
{
#define NO_COLLIDE 0
#define INTERSECTION 1
#define PENETRATION 2

	void computeCollisionDistance(Vector3 C1, Vector3 C2, Vector3 C3, qeal& R1, qeal& R2, qeal& R3, bool ss, qeal& dist, qeal& alpha, qeal& beta)
	{
		qeal A = C1.dot(C1) - R1 * R1;
		qeal B = 2.0 * (C1.dot(C2) - R1 * R2);
		qeal C = C2.dot(C2) - R2 * R2;
		qeal D = 2.0 * (C1.dot(C3) - R1 * R3);
		qeal E = 2.0 * (C2.dot(C3) - R2 * R3);
		qeal F = C3.dot(C3) - R3 * R3;

		qeal delta = 4 * A * C - B * B;

		alpha = 0.0; beta = 0.0;

		dist = valueOfQuadircSurface2D(alpha, beta, A, B, C, D, E, F);

		qeal temp_dist;
		qeal temp_alpha, temp_beta;

		temp_alpha = 1.0, temp_beta = 0.0;
		temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
		if (dist > temp_dist)
		{
			dist = temp_dist;
			alpha = 1.0; beta = 0.0;
		}

		temp_alpha = 0.0, temp_beta = 1.0;
		temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
		if (dist > temp_dist)
		{
			dist = temp_dist;
			alpha = 0.0; beta = 1.0;
		}

		temp_alpha = 0.0; temp_beta = -E / (2.0 *C);
		if (temp_beta > 0.0 && temp_beta < 1.0)
		{
			temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
			if (dist > temp_dist)
			{
				dist = temp_dist;
				alpha = temp_alpha; beta = temp_beta;
			}
		}

		temp_alpha = -D / (2.0 *A); temp_beta = 0.0;
		if (temp_alpha > 0.0 && temp_alpha < 1.0)
		{
			temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
			if (dist > temp_dist)
			{
				dist = temp_dist;
				alpha = temp_alpha; beta = temp_beta;
			}
		}

		if (ss)
		{
			temp_alpha = 0.5 * (2.0 * C + E - B - D) / (A - B + C); temp_beta = 1.0 - temp_alpha;
			if (temp_alpha > 0.0 && temp_alpha < 1.0)
			{
				temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
				if (dist > temp_dist)
				{
					dist = temp_dist;
					alpha = temp_alpha; beta = temp_beta;
				}
			}
		}
		else
		{
			temp_alpha = 1.0, temp_beta = 1.0;
			temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
			if (dist > temp_dist)
			{
				dist = temp_dist;
				alpha = 1.0; beta = 1.0;
			}

			temp_alpha = 1.0; temp_beta = -(B + E) / (2.0 *C);
			if (temp_beta > 0.0 && temp_beta < 1.0)
			{
				temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
				if (dist > temp_dist)
				{
					dist = temp_dist;
					alpha = temp_alpha; beta = temp_beta;
				}
			}

			temp_alpha = -(B + D) / (2.0 *A); temp_beta = 1.0;
			if (temp_alpha > 0.0 && temp_alpha < 1.0)
			{
				temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
				if (dist > temp_dist)
				{
					dist = temp_dist;
					alpha = temp_alpha; beta = temp_beta;
				}
			}
		}

		if (!IS_QEAL_ZERO(delta))
		{
			temp_alpha = (B * E - 2.0 * C * D) / delta; temp_beta = (B * D - 2.0 * A * E) / delta;

			if (temp_alpha > 0.0 && temp_alpha < 1.0 && temp_beta> 0.0 && temp_beta < 1.0 && temp_alpha + temp_beta < 1.0)
			{
				temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
				if (dist > temp_dist)
				{
					dist = temp_dist;
					alpha = temp_alpha; beta = temp_beta;
				}
			}
		}

	}


	qeal valueOfQuadircSurface2D(const qeal x, const qeal y, const qeal A, const qeal B, const qeal C, const qeal D, const qeal E, const qeal F)
	{
		return A * x * x + B * x * y + C * y * y + D * x + E * y + F;
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

}