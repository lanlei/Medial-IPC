#include "MPsCCD.cuh"
# 
#define IS_Numerical_ZERO(d) (CUDA_ABS(d) < 1e-8)

namespace MIPC
{

	__host__ void MPsCCD
	(
		int hostCollisionEventNum,
		int* devCollisionEventNum,
		qeal* devMedialPointPosition,
		qeal* devMedialPointRadius,
		qeal* devStaticMedialPointPosition,
		qeal* devStaticMedialPointRadius,
		qeal* devMedialPointMovingDir,
		int* devCollisionEventList,
		qeal* devCCD
	)
	{
		dim3 blockSize(128);
		uint32_t num_block = (hostCollisionEventNum + (128 - 1)) / 128;
		dim3 gridSize(num_block);

		MPsCCD << <gridSize, blockSize >> >
			(
				devCollisionEventNum,
				devMedialPointPosition,
				devMedialPointRadius,
				devStaticMedialPointPosition,
				devStaticMedialPointRadius,
				devMedialPointMovingDir,
				devCollisionEventList,
				devCCD
				);
		cudaDeviceSynchronize();
	}


	__global__ void MPsCCD
	(
		int* devCollisionEventNum,
		qeal* devMedialPointPosition,
		qeal* devMedialPointRadius,
		qeal* devStaticMedialPointPosition,
		qeal* devStaticMedialPointRadius,
		qeal* devMedialPointMovingDir,
		int* devCollisionEventList,
		qeal* devCCD
	)
	{
		__shared__ int collisionEventNum;
		__shared__ qeal blockC1[384]; // 3 * 128
		__shared__ qeal blockC2[384];
		__shared__ qeal blockC3[384];
		__shared__ qeal blockV1[384];
		__shared__ qeal blockV2[384];
		__shared__ qeal blockV3[384];
		__shared__ qeal blockR1[128]; // 1 * 128
		__shared__ qeal blockR2[128];
		__shared__ qeal blockR3[128];

		__shared__ qeal J_ABCDE[2304]; // 18 *128

		const int length = gridDim.x *  blockDim.x;
		int eventId = (blockIdx.x  * blockDim.x) + threadIdx.x;
		if (threadIdx.x == 0)
			collisionEventNum = *devCollisionEventNum;
		__syncthreads();
		if (eventId >= collisionEventNum)
			return;

		int offset1 = threadIdx.x;
		int offset3 = 3 * threadIdx.x;
		int offset18 = 18 * threadIdx.x;

		qeal* C1, *C2, *C3;
		qeal* V1, *V2, *V3;
		qeal* R1, *R2, *R3;

		C1 = blockC1 + offset3;
		C2 = blockC2 + offset3;
		C3 = blockC3 + offset3;
		V1 = blockV1 + offset3;
		V2 = blockV2 + offset3;
		V3 = blockV3 + offset3;
		R1 = blockR1 + offset1;
		R2 = blockR2 + offset1;
		R3 = blockR3 + offset1;

		qeal *JA_2, *JA_1, *JA_0;
		qeal *JB_2, *JB_1, *JB_0;
		qeal *JC_2, *JC_1, *JC_0;
		qeal *JD_2, *JD_1, *JD_0;
		qeal *JE_2, *JE_1, *JE_0;
		qeal *JF_2, *JF_1, *JF_0;

		JA_2 = J_ABCDE + offset18;
		JA_1 = J_ABCDE + offset18 + 1;
		JA_0 = J_ABCDE + offset18 + 2;

		JB_2 = J_ABCDE + offset18 + 3;
		JB_1 = J_ABCDE + offset18 + 4;
		JB_0 = J_ABCDE + offset18 + 5;

		JC_2 = J_ABCDE + offset18 + 6;
		JC_1 = J_ABCDE + offset18 + 7;
		JC_0 = J_ABCDE + offset18 + 8;

		JD_2 = J_ABCDE + offset18 + 9;
		JD_1 = J_ABCDE + offset18 + 10;
		JD_0 = J_ABCDE + offset18 + 11;

		JE_2 = J_ABCDE + offset18 + 12;
		JE_1 = J_ABCDE + offset18 + 13;
		JE_0 = J_ABCDE + offset18 + 14;

		JF_2 = J_ABCDE + offset18 + 15;
		JF_1 = J_ABCDE + offset18 + 16;
		JF_0 = J_ABCDE + offset18 + 17;
		__syncthreads();

		for (; eventId < collisionEventNum; eventId += length)
		{
			int flag = devCollisionEventList[5 * eventId];
			int mid[4];
			for (int i = 0; i < 4; i++)
				mid[i] = devCollisionEventList[5 * eventId + 1 + i];  //medial vertex index
			bool ss = false;

			if (flag == COLLISION_CC) // cc obj vs obj
			{
				for (int i = 0; i < 3; i++)
				{
					C1[i] = devMedialPointPosition[3 * mid[0] + i] - devMedialPointPosition[3 * mid[1] + i];
					C2[i] = devMedialPointPosition[3 * mid[3] + i] - devMedialPointPosition[3 * mid[2] + i];
					C3[i] = devMedialPointPosition[3 * mid[1] + i] - devMedialPointPosition[3 * mid[3] + i];

					V1[i] = devMedialPointMovingDir[3 * mid[0] + i] - devMedialPointMovingDir[3 * mid[1] + i];
					V2[i] = devMedialPointMovingDir[3 * mid[3] + i] - devMedialPointMovingDir[3 * mid[2] + i];
					V3[i] = devMedialPointMovingDir[3 * mid[1] + i] - devMedialPointMovingDir[3 * mid[3] + i];
				}
				R1[0] = devMedialPointRadius[mid[0]] - devMedialPointRadius[mid[1]];
				R2[0] = devMedialPointRadius[mid[2]] - devMedialPointRadius[mid[3]];
				R3[0] = devMedialPointRadius[mid[1]] + devMedialPointRadius[mid[3]];
			}
			else if (flag == COLLISION_SS)// ss obj vs obj
			{
				ss = true;
				for (int i = 0; i < 3; i++)
				{
					C1[i] = devMedialPointPosition[3 * mid[0] + i] - devMedialPointPosition[3 * mid[2] + i];
					C2[i] = devMedialPointPosition[3 * mid[1] + i] - devMedialPointPosition[3 * mid[2] + i];
					C3[i] = devMedialPointPosition[3 * mid[2] + i] - devMedialPointPosition[3 * mid[3] + i];

					V1[i] = devMedialPointMovingDir[3 * mid[0] + i] - devMedialPointMovingDir[3 * mid[2] + i];
					V2[i] = devMedialPointMovingDir[3 * mid[1] + i] - devMedialPointMovingDir[3 * mid[2] + i];
					V3[i] = devMedialPointMovingDir[3 * mid[2] + i] - devMedialPointMovingDir[3 * mid[3] + i];
				}
				R1[0] = devMedialPointRadius[mid[0]] - devMedialPointRadius[mid[2]];
				R2[0] = devMedialPointRadius[mid[1]] - devMedialPointRadius[mid[2]];
				R3[0] = devMedialPointRadius[mid[2]] + devMedialPointRadius[mid[3]];
			}
			else if (flag == COLLISION_DEFORMABLE_WITH_STATIC_CC) // cc obj vs static obj
			{
				for (int i = 0; i < 3; i++)
				{
					C1[i] = devMedialPointPosition[3 * mid[0] + i] - devMedialPointPosition[3 * mid[1] + i];
					C2[i] = devStaticMedialPointPosition[3 * mid[3] + i] - devStaticMedialPointPosition[3 * mid[2] + i];
					C3[i] = devMedialPointPosition[3 * mid[1] + i] - devStaticMedialPointPosition[3 * mid[3] + i];

					V1[i] = devMedialPointMovingDir[3 * mid[0] + i] - devMedialPointMovingDir[3 * mid[1] + i];
					V2[i] = 0.0;
					V3[i] = devMedialPointMovingDir[3 * mid[1] + i];
				}
				R1[0] = devMedialPointRadius[mid[0]] - devMedialPointRadius[mid[1]];
				R2[0] = devStaticMedialPointRadius[mid[2]] - devStaticMedialPointRadius[mid[3]];
				R3[0] = devMedialPointRadius[mid[1]] + devStaticMedialPointRadius[mid[3]];
			}
			else if (flag == COLLISION_DEFORMABLE_WITH_STATIC_SS) // ss obj vs static obj
			{
				ss = true;
				for (int i = 0; i < 3; i++)
				{
					C1[i] = devMedialPointPosition[3 * mid[0] + i] - devMedialPointPosition[3 * mid[2] + i];
					C2[i] = devMedialPointPosition[3 * mid[1] + i] - devMedialPointPosition[3 * mid[2] + i];
					C3[i] = devMedialPointPosition[3 * mid[2] + i] - devStaticMedialPointPosition[3 * mid[3] + i];

					V1[i] = devMedialPointMovingDir[3 * mid[0] + i] - devMedialPointMovingDir[3 * mid[2] + i];
					V2[i] = devMedialPointMovingDir[3 * mid[1] + i] - devMedialPointMovingDir[3 * mid[2] + i];
					V3[i] = devMedialPointMovingDir[3 * mid[2] + i];
				}
				R1[0] = devMedialPointRadius[mid[0]] - devMedialPointRadius[mid[2]];
				R2[0] = devMedialPointRadius[mid[1]] - devMedialPointRadius[mid[2]];
				R3[0] = devMedialPointRadius[mid[2]] + devStaticMedialPointRadius[mid[3]];
			}
			else if (flag == COLLISION_STATIC_WITH_DEFORMABLE_CC) // cc static obj vs obj
			{
				for (int i = 0; i < 3; i++)
				{
					C1[i] = devStaticMedialPointPosition[3 * mid[0] + i] - devStaticMedialPointPosition[3 * mid[1] + i];
					C2[i] = devMedialPointPosition[3 * mid[3] + i] - devMedialPointPosition[3 * mid[2] + i];
					C3[i] = devStaticMedialPointPosition[3 * mid[1] + i] - devMedialPointPosition[3 * mid[3] + i];

					V1[i] = 0.0;
					V2[i] = devMedialPointMovingDir[3 * mid[3] + i] - devMedialPointMovingDir[3 * mid[2] + i];
					V3[i] = -devMedialPointMovingDir[3 * mid[3] + i];
				}
				R1[0] = devStaticMedialPointRadius[mid[0]] - devStaticMedialPointRadius[mid[1]];
				R2[0] = devMedialPointRadius[mid[2]] - devMedialPointRadius[mid[3]];
				R3[0] = devStaticMedialPointRadius[mid[1]] + devMedialPointRadius[mid[3]];
			}
			else if (flag == COLLISION_STATIC_WITH_DEFORMABLE_SS) // ss static obj vs obj
			{
				ss = true;
				for (int i = 0; i < 3; i++)
				{
					C1[i] = devStaticMedialPointPosition[3 * mid[0] + i] - devStaticMedialPointPosition[3 * mid[2] + i];
					C2[i] = devStaticMedialPointPosition[3 * mid[1] + i] - devStaticMedialPointPosition[3 * mid[2] + i];
					C3[i] = devStaticMedialPointPosition[3 * mid[2] + i] - devMedialPointPosition[3 * mid[3] + i];

					V1[i] = 0.0;
					V2[i] = 0.0;
					V3[i] = -devMedialPointMovingDir[3 * mid[3] + i];
				}
				R1[0] = devStaticMedialPointRadius[mid[0]] - devStaticMedialPointRadius[mid[2]];
				R2[0] = devStaticMedialPointRadius[mid[1]] - devStaticMedialPointRadius[mid[2]];
				R3[0] = devStaticMedialPointRadius[mid[2]] + devMedialPointRadius[mid[3]];
			}

			qeal norrow = QEAL_ZERO;

			JA_2[0] = VectorDot(V1, V1);
			JA_1[0] = 2.0 * VectorDot(V1, C1);
			JA_0[0] = VectorDot(C1, C1) - R1[0] * R1[0];

			JB_2[0] = 2.0 *  VectorDot(V1, V2);
			JB_1[0] = 2.0 * (VectorDot(V1, C2) + VectorDot(V2, C1));
			JB_0[0] = 2.0 * (VectorDot(C1, C2) - R1[0] * R2[0]);

			JC_2[0] = VectorDot(V2, V2);
			JC_1[0] = 2.0 * VectorDot(V2, C2);
			JC_0[0] = VectorDot(C2, C2) - R2[0] * R2[0];

			JD_2[0] = 2.0 *  VectorDot(V1, V3);
			JD_1[0] = 2.0 * (VectorDot(V1, C3) + VectorDot(V3, C1));
			JD_0[0] = 2.0 * (VectorDot(C1, C3) - R1[0] * R3[0]);

			JE_2[0] = 2.0 *  VectorDot(V2, V3);
			JE_1[0] = 2.0 * (VectorDot(V2, C3) + VectorDot(V3, C2));
			JE_0[0] = 2.0 * (VectorDot(C2, C3) - R2[0] * R3[0]);

			JF_2[0] = VectorDot(V3, V3);
			JF_1[0] = 2.0 * VectorDot(V3, C3);
			JF_0[0] = VectorDot(C3, C3) - R3[0] * R3[0] - norrow;

			//JA_2[0] = Check_CUDA_ZERO(JA_2[0]); JA_1[0] = Check_CUDA_ZERO(JA_1[0]); JA_0[0] = Check_CUDA_ZERO(JA_0[0]);
			//JB_2[0] = Check_CUDA_ZERO(JB_2[0]); JB_1[0] = Check_CUDA_ZERO(JB_1[0]); JB_0[0] = Check_CUDA_ZERO(JB_0[0]);
			//JC_2[0] = Check_CUDA_ZERO(JC_2[0]); JC_1[0] = Check_CUDA_ZERO(JC_1[0]); JC_0[0] = Check_CUDA_ZERO(JC_0[0]);
			//JD_2[0] = Check_CUDA_ZERO(JD_2[0]); JD_1[0] = Check_CUDA_ZERO(JD_1[0]); JD_0[0] = Check_CUDA_ZERO(JD_0[0]);
			//JE_2[0] = Check_CUDA_ZERO(JE_2[0]); JE_1[0] = Check_CUDA_ZERO(JE_1[0]); JE_0[0] = Check_CUDA_ZERO(JE_0[0]);
			//JF_2[0] = Check_CUDA_ZERO(JF_2[0]); JF_1[0] = Check_CUDA_ZERO(JF_1[0]); JF_0[0] = Check_CUDA_ZERO(JF_0[0]);

			bool isCollide = false;
			qeal ftc = 1.0;
			int collideType = 0; // no colliding
			if (dcd(C1, V1, C2, V2, C3, V3, R1[0], R2[0], R3[0], ss, norrow))
			{
				//devCCD[eventId] = 1.0; // ill-condition
				//return;
			}
			//
			qeal ftcEndpoints = 1.0;
			bool collideEndpoints = checkEndpointAlphaBetaCCD(eventId, *JA_2, *JA_1, *JA_0, *JB_2, *JB_1, *JB_0, *JC_2, *JC_1, *JC_0, *JD_2, *JD_1, *JD_0, *JE_2, *JE_1, *JE_0, *JF_2, *JF_1, *JF_0, ftcEndpoints, ss, norrow);
			if (collideEndpoints && ftcEndpoints < ftc)
			{
				isCollide = true;
				ftc = ftcEndpoints;
			}
			//
			qeal ftcA0 = 1.0;
			bool collideA0 = checkAlphaIsZeroCCD(eventId, *JA_2, *JA_1, *JA_0, *JB_2, *JB_1, *JB_0, *JC_2, *JC_1, *JC_0, *JD_2, *JD_1, *JD_0, *JE_2, *JE_1, *JE_0, *JF_2, *JF_1, *JF_0, ftcA0, ss, norrow);
			if (collideA0 && ftcA0 < ftc)
			{
				isCollide = true;
				ftc = ftcA0;
			}
			//
			qeal ftcB0 = 1.0;
			bool collideB0 = checkBetaIsZeroCCD(eventId, *JA_2, *JA_1, *JA_0, *JB_2, *JB_1, *JB_0, *JC_2, *JC_1, *JC_0, *JD_2, *JD_1, *JD_0, *JE_2, *JE_1, *JE_0, *JF_2, *JF_1, *JF_0, ftcB0, ss, norrow);
			if (collideB0 && ftcB0 < ftc)
			{
				isCollide = true;
				ftc = ftcB0;
			}

			//
			if (ss)
			{
				qeal ftcAB1 = 1.0;
				bool collideAB1 = checkAlphaPlusBetaIsOneCCD(eventId, *JA_2, *JA_1, *JA_0, *JB_2, *JB_1, *JB_0, *JC_2, *JC_1, *JC_0, *JD_2, *JD_1, *JD_0, *JE_2, *JE_1, *JE_0, *JF_2, *JF_1, *JF_0, ftcAB1, ss, norrow);
				if (collideAB1 && ftcAB1 < ftc)
				{
					isCollide = true;
					ftc = ftcAB1;
				}
			}
			else
			{
				qeal ftcA1 = 1.0;
				bool collideA1 = checkAlphaIsOneCCD(eventId, *JA_2, *JA_1, *JA_0, *JB_2, *JB_1, *JB_0, *JC_2, *JC_1, *JC_0, *JD_2, *JD_1, *JD_0, *JE_2, *JE_1, *JE_0, *JF_2, *JF_1, *JF_0, ftcA1, ss, norrow);
				if (collideA1 && ftcA1 < ftc)
				{
					isCollide = true;
					ftc = ftcA1;
				}

				qeal ftcB1 = 1.0;
				bool collideB1 = checkBetaIsOneCCD(eventId, *JA_2, *JA_1, *JA_0, *JB_2, *JB_1, *JB_0, *JC_2, *JC_1, *JC_0, *JD_2, *JD_1, *JD_0, *JE_2, *JE_1, *JE_0, *JF_2, *JF_1, *JF_0, ftcB1, ss, norrow);
				if (collideB1 && ftcB1 < ftc)
				{
					isCollide = true;
					ftc = ftcB1;
				}
			}

			qeal ftcAB = 1.0;
			bool collideAB = checkAlphaBetaCCD(eventId, *JA_2, *JA_1, *JA_0, *JB_2, *JB_1, *JB_0, *JC_2, *JC_1, *JC_0, *JD_2, *JD_1, *JD_0, *JE_2, *JE_1, *JE_0, *JF_2, *JF_1, *JF_0, ftcAB, collideType, ss, norrow);

			if (collideAB && ftcAB < ftc)
			{
				isCollide = true;
				ftc = ftcAB;
			}

			//
			devCCD[eventId] = ftc;
		}
	}


	__device__ __forceinline__
		bool dcd(qeal * C1, qeal * V1, qeal * C2, qeal * V2, qeal * C3, qeal * V3, qeal& R1, qeal& R2, qeal& R3, bool is_ss, qeal norrow)
	{
		qeal CV1[3], CV2[3], CV3[3];
		for (int i = 0; i < 3; i++)
		{
			CV1[i] = C1[i] + V1[i];
			CV2[i] = C2[i] + V2[i];
			CV3[i] = C3[i] + V3[i];
		}

		qeal A = VectorDot(CV1, CV1) - R1 * R1;
		qeal B = 2.0 * (VectorDot(CV1, CV2) - R1 * R2);
		qeal C = VectorDot(CV2, CV2) - R2 * R2;
		qeal D = 2.0 * (VectorDot(CV1, CV3) - R1 * R3);
		qeal E = 2.0 * (VectorDot(CV2, CV3) - R2 * R3);
		qeal F = VectorDot(CV3, CV3) - R3 * R3;

		return dcd(A, B, C, D, E, F, is_ss, norrow);
	}

	__device__ __forceinline__
		bool dcd(qeal* C1, qeal* C2, qeal* C3, qeal& R1, qeal& R2, qeal& R3, bool is_ss, qeal norrow)
	{
		qeal A = VectorDot(C1, C1) - R1 * R1;
		qeal B = 2.0 * (VectorDot(C1, C2) - R1 * R2);
		qeal C = VectorDot(C2, C2) - R2 * R2;
		qeal D = 2.0 * (VectorDot(C1, C3) - R1 * R3);
		qeal E = 2.0 * (VectorDot(C2, C3) - R2 * R3);
		qeal F = VectorDot(C3, C3) - R3 * R3;

		return dcd(A, B, C, D, E, F, is_ss, norrow);
	}

	__device__ __forceinline__
		bool dcd(qeal A, qeal B, qeal C, qeal D, qeal E, qeal F, bool is_ss, qeal norrow)
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
			if (!IS_CUDA_ZERO(A))
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
			if (!IS_CUDA_ZERO(4 * A*C - B * B))
			{
				qeal x = (B*E - 2 * C*D) / (4 * A*C - B * B);
				qeal y = (B*D - 2 * A*E) / (4 * A*C - B * B);
				if (x > 0 && x < 1 && y > 0 && y < 1)
				{
					if (!is_ss)
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
			if (!is_ss)
			{
				if (!IS_CUDA_ZERO(A))
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
			if (!is_ss)
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
			if (!is_ss)
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
		min = Check_CUDA_ZERO(min);

		if (min > 0.0)
			return false;

		return true;
	}

	__device__ __forceinline__
		void computeCollisionDistance(qeal * C1, qeal * C2, qeal * C3, qeal& R1, qeal& R2, qeal& R3, bool ss, qeal& dist, qeal& alpha, qeal& beta)
	{
		qeal A = VectorDot(C1, C1) - R1 * R1;
		qeal B = 2.0 * (VectorDot(C1, C2) - R1 * R2);
		qeal C = VectorDot(C2, C2) - R2 * R2;
		qeal D = 2.0 * (VectorDot(C1, C3) - R1 * R3);
		qeal E = 2.0 * (VectorDot(C2, C3) - R2 * R3);
		qeal F = VectorDot(C3, C3) - R3 * R3;

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
			alpha = temp_alpha; beta = temp_beta;
		}

		temp_alpha = 0.0, temp_beta = 1.0;
		temp_dist = valueOfQuadircSurface2D(temp_alpha, temp_beta, A, B, C, D, E, F);
		if (dist > temp_dist)
		{
			dist = temp_dist;
			alpha = temp_alpha; beta = temp_beta;
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

		if (delta != 0.0)
		{
			temp_alpha = (B * E - 2.0 * C * D) / delta; temp_beta = (B * D - 2.0 * A * E) / delta;
						
			if (temp_alpha > 0.0 && temp_alpha < 1.0 && temp_beta> 0.0 && temp_beta < 1.0 && (ss == 1 ? temp_alpha + temp_beta < 1.0: true))
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

	__device__ __forceinline__
		bool checkEndpointAlphaBetaCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		bool collide = false;
		//set alpha = 0.0, beta = 0.0
		//F(t) <= 0
		qeal pos_00[3];
		pos_00[RANGE_FLAG] = 1;
		pos_00[RANGE_MIN_INDEX] = 0;
		pos_00[RANGE_MAX_INDEX] = 1;
		qeal neg_00[3];
		neg_00[RANGE_FLAG] = 1;
		neg_00[RANGE_MIN_INDEX] = 0;
		neg_00[RANGE_MAX_INDEX] = 1;
		qeal roots_00[2];
		int rootCount_00 = 0;
		J1_2 = JF_2; J1_1 = JF_1; J1_0 = JF_0;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_00, rootCount_00, pos_00, neg_00);

		qeal ft00 = -1;
		if (rootCount_00 > 0)
			ft00 = roots_00[0];
		else if (neg_00[RANGE_FLAG] == 1)
			ft00 = neg_00[RANGE_MIN_INDEX];

		if (ft00 > 0.0 && ft00 < 1.0)
		{
			collide = true;
			if (ft00 < ftc) ftc = ft00;
		}

		// set alpha = 0.0, beta = 1.0
		// C(t) + E(t) + F(t) <= 0
		qeal pos_01[3];
		pos_01[RANGE_FLAG] = 1;
		pos_01[RANGE_MIN_INDEX] = 0;
		pos_01[RANGE_MAX_INDEX] = 1;
		qeal neg_01[3];
		neg_01[RANGE_FLAG] = 1;
		neg_01[RANGE_MIN_INDEX] = 0;
		neg_01[RANGE_MAX_INDEX] = 1;
		qeal roots_01[2];
		int rootCount_01 = 0;
		J1_2 = JC_2 + JE_2 + JF_2;
		J1_1 = JC_1 + JE_1 + JF_1;
		J1_0 = JC_0 + JE_0 + JF_0;

		solveQuadricNEq(J1_2, J1_1, J1_0, roots_01, rootCount_01, pos_01, neg_01);

		qeal ft01 = -1;
		if (rootCount_01 > 0)
			ft01 = roots_01[0];
		else if (neg_01[RANGE_FLAG] == 1)
			ft01 = neg_01[RANGE_MIN_INDEX];
		if (ft01 > 0.0 && ft01 < 1.0)
		{
			collide = true;
			if (ft01 < ftc) ftc = ft01;
		}

		// set alpha = 1.0, beta = 0.0
		//A(t) + D(t) + F(t) <= 0
		qeal pos_10[3];
		pos_10[RANGE_FLAG] = 1;
		pos_10[RANGE_MIN_INDEX] = 0;
		pos_10[RANGE_MAX_INDEX] = 1;
		qeal neg_10[3];
		neg_10[RANGE_FLAG] = 1;
		neg_10[RANGE_MIN_INDEX] = 0;
		neg_10[RANGE_MAX_INDEX] = 1;
		qeal roots_10[2];
		int rootCount_10 = 0;
		J1_2 = JA_2 + JD_2 + JF_2;
		J1_1 = JA_1 + JD_1 + JF_1;
		J1_0 = JA_0 + JD_0 + JF_0;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_10, rootCount_10, pos_10, neg_10);

		qeal ft10 = -1;
		if (rootCount_10 > 0)
			ft10 = roots_10[0];
		else if (neg_10[RANGE_FLAG] == 1)
			ft10 = neg_10[RANGE_MIN_INDEX];
		if (ft10 >= 0.0 && ft10 <= 1.0)
		{
			collide = true;
			if (ft10 < ftc) ftc = ft10;
		}

		if (!is_ss)
		{
			// set alpha = 1.0, beta = 1.0
			//A(t) + B(t) + C(t) + D(t) + E(t) + F(t) <= 0
			qeal pos_11[3];
			pos_11[RANGE_FLAG] = 1;
			pos_11[RANGE_MIN_INDEX] = 0;
			pos_11[RANGE_MAX_INDEX] = 1;
			qeal neg_11[3];
			neg_11[RANGE_FLAG] = 1;
			neg_11[RANGE_MIN_INDEX] = 0;
			neg_11[RANGE_MAX_INDEX] = 1;
			qeal roots_11[2];
			int rootCount_11 = 0;
			J1_2 = JA_2 + JB_2 + JC_2 + JD_2 + JE_2 + JF_2;
			J1_1 = JA_1 + JB_1 + JC_1 + JD_1 + JE_1 + JF_1;
			J1_0 = JA_0 + JB_0 + JC_0 + JD_0 + JE_0 + JF_0;
			solveQuadricNEq(J1_2, J1_1, J1_0, roots_11, rootCount_11, pos_11, neg_11);
			qeal ft11 = -1;
			if (rootCount_11 > 0)
				ft11 = roots_11[0];
			else if (neg_11[RANGE_FLAG] == 1)
				ft11 = neg_11[RANGE_MIN_INDEX];
			if (ft11 > 0.0 && ft11 <= 1.0)
			{
				collide = true;
				if (ft11 < ftc) ftc = ft11;
			}
		}
		return collide;
	}

	__device__ __forceinline__
		bool checkAlphaIsZeroCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal pos_0b[3];
		pos_0b[RANGE_FLAG] = 1;
		pos_0b[RANGE_MIN_INDEX] = 0;
		pos_0b[RANGE_MAX_INDEX] = 1;
		qeal neg_0b[3];
		neg_0b[RANGE_FLAG] = 1;
		neg_0b[RANGE_MIN_INDEX] = 0;
		neg_0b[RANGE_MAX_INDEX] = 1;
		qeal roots_0b[6];
		int rootCount_0b = 0;

		// -E > 0
		J1_2 = -JE_2;
		J1_1 = -JE_1;
		J1_0 = -JE_0;

		if (IS_CUDA_ZERO(J1_2) && IS_CUDA_ZERO(J1_1) && IS_CUDA_ZERO(J1_0))
			return false;

		solveQuadricNEq(J1_2, J1_1, J1_0, roots_0b, rootCount_0b, pos_0b, neg_0b);
		if (pos_0b[RANGE_FLAG] == 0)
			return false;

		//2C + E > 0
		J1_2 = 2.0 * JC_2 + JE_2;
		J1_1 = 2.0 * JC_1 + JE_1;
		J1_0 = 2.0 * JC_0 + JE_0;
		if (IS_CUDA_ZERO(J1_2) && IS_CUDA_ZERO(J1_1) && IS_CUDA_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_0b, rootCount_0b, pos_0b, neg_0b);
		if (pos_0b[RANGE_FLAG] == 0)
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
		if (pos_0b[RANGE_FLAG] == 0)
			return false;

		qeal ft = pos_0b[RANGE_MIN_INDEX];

		if (ft > 0.0 && ft < ftc)
		{
			ftc = ft;
		}
		return true;
	}

	__device__ __forceinline__
		bool checkBetaIsZeroCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal pos_a0[3];
		pos_a0[RANGE_FLAG] = 1;
		pos_a0[RANGE_MIN_INDEX] = 0;
		pos_a0[RANGE_MAX_INDEX] = 1;
		qeal neg_a0[3];
		neg_a0[RANGE_FLAG] = 1;
		neg_a0[RANGE_MIN_INDEX] = 0;
		neg_a0[RANGE_MAX_INDEX] = 1;
		qeal roots_a0[6];
		int rootCount_a0 = 0;
		// -D > 0
		J1_2 = -JD_2;
		J1_1 = -JD_1;
		J1_0 = -JD_0;
		if (IS_CUDA_ZERO(J1_2) && IS_CUDA_ZERO(J1_1) && IS_CUDA_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_a0, rootCount_a0, pos_a0, neg_a0);
		if (pos_a0[RANGE_FLAG] == 0)
			return false;

		//2A + D > 0
		J1_2 = 2.0 * JA_2 + JD_2;
		J1_1 = 2.0 * JA_1 + JD_1;
		J1_0 = 2.0 * JA_0 + JD_0;
		if (IS_CUDA_ZERO(J1_2) && IS_CUDA_ZERO(J1_1) && IS_CUDA_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_a0, rootCount_a0, pos_a0, neg_a0);
		if (pos_a0[RANGE_FLAG] == 0)
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
		if (pos_a0[RANGE_FLAG] == 0)
			return false;

		qeal ft = pos_a0[RANGE_MIN_INDEX];
		// Numerical safe
		if (ft > 0.0 && ft < ftc)
		{
			ftc = ft;
		}
		return true;
	}

	__device__ __forceinline__
		bool checkAlphaIsOneCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal J4_2, J4_1, J4_0;
		qeal pos_1b[3];
		pos_1b[RANGE_FLAG] = 1;
		pos_1b[RANGE_MIN_INDEX] = 0;
		pos_1b[RANGE_MAX_INDEX] = 1;
		qeal neg_1b[3];
		neg_1b[RANGE_FLAG] = 1;
		neg_1b[RANGE_MIN_INDEX] = 0;
		neg_1b[RANGE_MAX_INDEX] = 1;

		qeal roots_1b[6];
		int rootCount_1b = 0;

		// -(B+ E) > 0
		J1_2 = -(JB_2 + JE_2);
		J1_1 = -(JB_1 + JE_1);
		J1_0 = -(JB_0 + JE_0);
		if (IS_CUDA_ZERO(J1_2) && IS_CUDA_ZERO(J1_1) && IS_CUDA_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_1b, rootCount_1b, pos_1b, neg_1b);
		if (pos_1b[RANGE_FLAG] == 0)
			return false;

		// 2C + B + E > 0
		J1_2 = 2.0 * JC_2 + JB_2 + JE_2;
		J1_1 = 2.0 * JC_1 + JB_1 + JE_1;
		J1_0 = 2.0 * JC_0 + JB_0 + JE_0;
		if (IS_CUDA_ZERO(J1_2) && IS_CUDA_ZERO(J1_1) && IS_CUDA_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_1b, rootCount_1b, pos_1b, neg_1b);
		if (pos_1b[RANGE_FLAG] == 0)
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
		qeal ft = pos_1b[RANGE_MIN_INDEX];
		if (pos_1b[RANGE_FLAG] == 0)
			return false;
		if (ft > 0.0 && ft < ftc)
		{
			ftc = ft;
		}
		return true;
	}

	__device__ __forceinline__
		bool checkBetaIsOneCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal J4_2, J4_1, J4_0;

		qeal pos_a1[3];
		pos_a1[RANGE_FLAG] = 1;
		pos_a1[RANGE_MIN_INDEX] = 0;
		pos_a1[RANGE_MAX_INDEX] = 1;
		qeal neg_a1[3];
		neg_a1[RANGE_FLAG] = 1;
		neg_a1[RANGE_MIN_INDEX] = 0;
		neg_a1[RANGE_MAX_INDEX] = 1;
		qeal roots_a1[6];
		int rootCount_a1 = 0;

		// -(B+ D) > 0
		J1_2 = -(JB_2 + JD_2);
		J1_1 = -(JB_1 + JD_1);
		J1_0 = -(JB_0 + JD_0);
		if (IS_CUDA_ZERO(J1_2) && IS_CUDA_ZERO(J1_1) && IS_CUDA_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_a1, rootCount_a1, pos_a1, neg_a1);
		if (pos_a1[RANGE_FLAG] == 0)
			return false;

		// 2A + B + D > 0
		J1_2 = 2.0 * JA_2 + JB_2 + JD_2;
		J1_1 = 2.0 * JA_1 + JB_1 + JD_1;
		J1_0 = 2.0 * JA_0 + JB_0 + JD_0;
		if (IS_CUDA_ZERO(J1_2) && IS_CUDA_ZERO(J1_1) && IS_CUDA_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_a1, rootCount_a1, pos_a1, neg_a1);
		if (pos_a1[RANGE_FLAG] == 0)
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
		if (pos_a1[RANGE_FLAG] == 0)
			return false;

		qeal ft = pos_a1[RANGE_MIN_INDEX];
		if (ft > 0.0 && ft < ftc)
		{
			ftc = ft;
		}

		return true;
	}

	__device__ __forceinline__
		bool checkAlphaPlusBetaIsOneCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, bool is_ss, qeal norrow)
	{
		qeal J1_2, J1_1, J1_0;
		qeal J2_2, J2_1, J2_0;
		qeal J3_2, J3_1, J3_0;
		qeal J4_2, J4_1, J4_0;

		//
		qeal pos_ab1[3];
		pos_ab1[RANGE_FLAG] = 1;
		pos_ab1[RANGE_MIN_INDEX] = 0;
		pos_ab1[RANGE_MAX_INDEX] = 1;
		qeal neg_ab1[3];
		neg_ab1[RANGE_FLAG] = 1;
		neg_ab1[RANGE_MIN_INDEX] = 0;
		neg_ab1[RANGE_MAX_INDEX] = 1;
		qeal roots_ab1[6];
		int rootCount_ab1 = 0;

		// -(B+ D - 2C - E) > 0
		J1_2 = -(JB_2 + JD_2 - 2.0 * JC_2 - JE_2);
		J1_1 = -(JB_1 + JD_1 - 2.0 * JC_1 - JE_1);
		J1_0 = -(JB_0 + JD_0 - 2.0 * JC_0 - JE_0);
		if (IS_CUDA_ZERO(J1_2) && IS_CUDA_ZERO(J1_1) && IS_CUDA_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_ab1, rootCount_ab1, pos_ab1, neg_ab1);
		if (pos_ab1[RANGE_FLAG] == 0)
			return false;

		// 2.0 *  (A + C - B) + (B + D - 2.0 * C - E) > 0
		J1_2 = (2.0 * JA_2 + JD_2) - (JB_2 + JE_2);
		J1_1 = (2.0 * JA_1 + JD_1) - (JB_1 + JE_1);
		J1_0 = (2.0* JA_0 + JD_0) - (JB_0 + JE_0);

		if (IS_CUDA_ZERO(J1_2) && IS_CUDA_ZERO(J1_1) && IS_CUDA_ZERO(J1_0))
			return false;
		solveQuadricNEq(J1_2, J1_1, J1_0, roots_ab1, rootCount_ab1, pos_ab1, neg_ab1);
		if (pos_ab1[RANGE_FLAG] == 0)
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
		if (pos_ab1[RANGE_FLAG] == 0)
			return false;
		qeal ft = pos_ab1[RANGE_MIN_INDEX];

		if (ft > 0.0 && ft < ftc)
		{
			ftc = ft;
		}
		return true;
	}

	__device__ __forceinline__
		bool checkAlphaBetaCCD(int id, qeal& JA_2, qeal& JA_1, qeal& JA_0, qeal&  JB_2, qeal&  JB_1, qeal&  JB_0, qeal&  JC_2, qeal&  JC_1, qeal& JC_0, qeal&  JD_2, qeal&  JD_1, qeal&  JD_0, qeal&  JE_2, qeal&  JE_1, qeal&  JE_0, qeal&  JF_2, qeal&  JF_1, qeal&  JF_0, qeal& ftc, int collideType, bool is_ss, qeal norrow)
	{
		bool collide = false;
		qeal pos[3];
		pos[RANGE_FLAG] = 1;
		pos[RANGE_MIN_INDEX] = 0;
		pos[RANGE_MAX_INDEX] = 1;
		qeal neg[3];
		neg[RANGE_FLAG] = 1;
		neg[RANGE_MIN_INDEX] = 0;
		neg[RANGE_MAX_INDEX] = 1;
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

		for (int i = 0; i < rootCount_delta; i++)
		{
			qeal ft = roots_delta[i];
			qeal At = JA_2 * ft * ft + JA_1 * ft + JA_0;
			qeal Bt = JB_2 * ft * ft + JB_1 * ft + JB_0;
			qeal Ct = JC_2 * ft * ft + JC_1 * ft + JC_0;
			qeal Dt = JD_2 * ft * ft + JD_1 * ft + JD_0;
			qeal Et = JE_2 * ft * ft + JE_1 * ft + JE_0;
			qeal Ft = JF_2 * ft * ft + JF_1 * ft + JF_0;
			if (dcd(At, Bt, Ct, Dt, Et, Ft, is_ss, norrow) && ft < ftc)
			{
				collide = true;
				ftc = ft;
			}
		}

		//
		//BE - 2CD
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
		if (pos[RANGE_FLAG] == 0)
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
		if (pos[RANGE_FLAG] == 0)
			return collide;

		if (!is_ss)
		{
			// (4AC - B^2) -  (BE - 2CD)
			qeal roots_ga1[6];
			int rootCount_ga1 = 0;
			
			P_4 = H_4 - L_4;
			P_3 = H_3 - L_3;
			P_2 = H_2 - L_2;
			P_1 = H_1 - L_1;
			P_0 = H_0 - L_0;
			solveQuarticNEq	(P_4, P_3, P_2, P_1, P_0, roots_ga1, rootCount_ga1, pos, neg);
			if (pos[RANGE_FLAG] == 0)
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
			if (pos[RANGE_FLAG] == 0)
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
			if (pos[RANGE_FLAG] == 0)
				return collide;
		}

		qeal ft = 1.0;
		qeal W_6 = 0, W_5 = 0, W_4 = 0, W_3 = 0, W_2 = 0, W_1 = 0, W_0 = 0;

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

		W_6 = Check_CUDA_ZERO(W_6);
		W_5 = Check_CUDA_ZERO(W_5);
		W_4 = Check_CUDA_ZERO(W_4);
		W_3 = Check_CUDA_ZERO(W_3);
		W_2 = Check_CUDA_ZERO(W_2);
		W_1 = Check_CUDA_ZERO(W_1);
		W_0 = Check_CUDA_ZERO(W_0);
		if (CUDA_ABS(W_6) < 1e-8 &&  CUDA_ABS(W_5) < 1e-8 &&  CUDA_ABS(W_4) < 1e-8 &&  CUDA_ABS(W_3) < 1e-8 &&  CUDA_ABS(W_2) < 1e-8 &&  CUDA_ABS(W_1) < 1e-8)
		{
			W_6 *= 1e8; W_5 *= 1e8; W_4 *= 1e8; W_3 *= 1e8; W_2 *= 1e8; W_1 *= 1e8; W_0 *= 1e8;
		}
		ft = SolveSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, rootCount_delta, roots_delta, pos, norrow);
		if (ft > 0.0 && ft < ftc)
		{
			// double check
			qeal A = JA_2 * ft * ft + JA_1 * ft + JA_0;
			qeal B = JB_2 * ft * ft + JB_1 * ft + JB_0;
			qeal C = JC_2 * ft * ft + JC_1 * ft + JC_0;
			qeal D = JD_2 * ft * ft + JD_1 * ft + JD_0;
			qeal E = JE_2 * ft * ft + JE_1 * ft + JE_0;
			qeal F = JF_2 * ft * ft + JF_1 * ft + JF_0;
			qeal delta = 4 * A * C - B * B;
			qeal alpha, beta;
			alpha = Check_CUDA_ZERO((B * E - 2.0 * C * D) / delta); beta = Check_CUDA_ZERO((B * D - 2.0 * A * E) / delta);
			if (delta != 0.0 && alpha >= 0 && alpha <= 1.0 && beta >= 0.0 && beta <= 1.0 && (is_ss == 1 ? alpha + beta < 1.0 : true))
			{
				ftc = ft;
				collide = true;
			}
		}
		return collide;
	}

	__device__ __forceinline__
		qeal valueOfQuadircSurface2D(qeal x, qeal y, qeal& A, qeal& B, qeal& C, qeal& D, qeal& E, qeal& F)
	{
		return A * x * x + B * x * y + C * y * y + D * x + E * y + F;
	}

	__device__ __forceinline__
		qeal computeQuadricEquation(qeal x, qeal& a, qeal& b, qeal& c)
	{
		return a * x * x + b * x + c;
	}

	__device__ __forceinline__
		qeal computeQuarticEquation(qeal x, qeal& a, qeal& b, qeal& c, qeal& d, qeal& e)
	{
		return a * x * x * x * x + b * x * x * x + c * x * x + d * x + e;
	}

	__device__ __forceinline__
		qeal computeSexticEquation(qeal x, qeal& a, qeal& b, qeal& c, qeal& d, qeal& e, qeal& f, qeal& g)
	{
		return a * x * x * x * x * x * x+ b * x * x * x * x * x+ c * x * x * x * x + d * x * x * x+ e * x * x + f * x + g;
	}

	__device__ __forceinline__
		void solveQuadricNEq(qeal& a, qeal& b, qeal& c, qeal * x, int& countRoot, qeal* posRange, qeal* negRange, qeal mini, qeal maxi)
	{
		countRoot = 0;
		if (IS_Numerical_ZERO(a) && IS_Numerical_ZERO(b))
		{
			if (c > 0.0)
				negRange[RANGE_FLAG] = 0;
			else if (c < 0.0)
				posRange[RANGE_FLAG] = 0;
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
				v = Check_CUDA_ZERO(v);
				if (v >= 0.0)
				{
					overlapRange(posRange, lx, hx);
					posInsert = true;
				}
				else
				{
					overlapRange(negRange, lx, hx);
					negInsert = true;
				}

			}
			lx = x[0];
			if (countRoot > 1)
				hx = x[1];
			else hx = maxi;
			qeal midx = lx + (hx - lx) / 2;
			qeal v = computeQuadricEquation(midx, a, b, c);
			v = computeQuadricEquation(midx, a, b, c);
			if (v >= 0.0)
			{
				overlapRange(posRange, lx, hx);
				posInsert = true;
			}
			else
			{
				overlapRange(negRange, lx, hx);
				negInsert = true;
			}
			if (!posInsert)posRange[RANGE_FLAG] = 0;
			if (!negInsert)negRange[RANGE_FLAG] = 0;
		}
		else
		{
			qeal midx = mini + (maxi - mini) / 2.0;
			qeal v = computeQuadricEquation(midx, a, b, c);
			if (v >= 0.0)
			{
				overlapRange(posRange, mini, maxi);
				negRange[RANGE_FLAG] = 0;
			}
			else
			{
				overlapRange(negRange, mini, maxi);
				posRange[RANGE_FLAG] = 0;
			}
		}
	}

	__device__ __forceinline__
		bool solveQuadricEquation(qeal& a, qeal& b, qeal& c, qeal* x, int& countRoot, bool banSame, qeal mini, qeal maxi)
	{
		countRoot = 0;
		qeal temp_x[2];
		qeal root;
		if (IS_Numerical_ZERO(a))
		{
			if (IS_Numerical_ZERO(b))
			{
				if (IS_Numerical_ZERO(c))
				{
					temp_x[0] = -1.0;
					temp_x[1] = -1.0;
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
			if (IS_Numerical_ZERO(delta))
			{
				root = -b_ / (2.0 * a_);
				temp_x[countRoot++] = root;
				if (!banSame)
					temp_x[countRoot++] = root;
			}
			else if (delta > 0)
			{
				if (IS_Numerical_ZERO(c_))
				{
					root = -1; //invalid value
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

	__device__ __forceinline__
		void solveQuarticNEq(qeal&  a, qeal& b, qeal& c, qeal& d, qeal& e, qeal * x, int& countRoot, qeal * posRange, qeal * negRange, qeal mini, qeal maxi)
	{
		countRoot = 0;
		if (IS_Numerical_ZERO(a) && IS_Numerical_ZERO(b))
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
				v = Check_CUDA_ZERO(v);
				if (v >= 0.0)
				{
					overlapRange(posRange, lx, hx);
					posInsert = true;
				}
				else
				{
					overlapRange(negRange, lx, hx);
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
			v = Check_CUDA_ZERO(v);
			if (v >= 0.0)
			{
				overlapRange(posRange, lx, hx);
				posInsert = true;
			}
			else
			{
				overlapRange(negRange, lx, hx);
				negInsert = true;
			}

			if (!posInsert)posRange[RANGE_FLAG] = 0;
			if (!negInsert)negRange[RANGE_FLAG] = 0;
		}
		else
		{
			qeal midx = mini + (maxi - mini) / 2.0;
			qeal v = computeQuarticEquation(midx, a, b, c, d, e);
			v = Check_CUDA_ZERO(v);
			if (v >= 0.0)
			{
				overlapRange(posRange, mini, maxi);
				negRange[RANGE_FLAG] = 0;
			}
			else
			{
				overlapRange(negRange, mini, maxi);
				posRange[RANGE_FLAG] = 0;
			}
		}
	}

	__device__ __forceinline__
		bool solveCubicEquation(qeal& a, qeal& b, qeal& c, qeal& d, qeal* x, int& countRoot, bool banSame, qeal mini, qeal maxi)
	{
		countRoot = 0;
		qeal temp_x[3];
		qeal root;

		if (IS_Numerical_ZERO(a))
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

			if (IS_Numerical_ZERO(A) && IS_Numerical_ZERO(B))
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
			else if (IS_Numerical_ZERO(delta))
			{
				if (!IS_Numerical_ZERO(A))
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
		qeal temp;
		for (int i = 0; i < countRoot; i++)
		{
			qeal p = temp_x[i];
			for (int j = i + 1; j < countRoot; j++)
			{
				qeal q = temp_x[j];
				if (p > q)
				{
					temp = temp_x[i];
					temp_x[i] = temp_x[j];
					temp_x[j] = temp;
					p = q;
				}
			}
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

	__device__ __forceinline__
		bool solveQuarticEquation(qeal& a, qeal& b, qeal& c, qeal& d, qeal& e, qeal* x, int& countRoot, bool banSame, qeal mini, qeal maxi)
	{
		countRoot = 0;
		qeal temp_x[3];

		if (IS_Numerical_ZERO(a))
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
			if (!solveCubicEquation(J0, J1, J2, J3, Jx, Jc, false, -QEAL_MAX, QEAL_MAX))
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
		qeal temp;
		for (int i = 0; i < countRoot; i++)
		{
			qeal p = temp_x[i];
			for (int j = i + 1; j < countRoot; j++)
			{
				qeal q = temp_x[j];
				if (p > q)
				{
					temp = temp_x[i];
					temp_x[i] = temp_x[j];
					temp_x[j] = temp;
					p = q;
				}
			}
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

	__device__ __forceinline__
		bool overlapRange(qeal* a, qeal mini, qeal maxi)
	{
		if (a[RANGE_FLAG] == 0)
			return false;
		if (a[RANGE_MAX_INDEX] < mini || a[RANGE_MIN_INDEX] > maxi)
			return false;

		qeal set[2];
		if (a[RANGE_MIN_INDEX] >= mini)
			set[0] = a[RANGE_MIN_INDEX];
		else set[0] = mini;

		if (a[RANGE_MAX_INDEX] >= maxi)
			set[1] = a[RANGE_MAX_INDEX];
		else set[1] = maxi;

		a[RANGE_MIN_INDEX] = set[0];
		a[RANGE_MAX_INDEX] = set[1];
		return true;
	}

	__device__ __forceinline__
		void genQuarticCoeffs(qeal& P_4, qeal& P_3, qeal& P_2, qeal& P_1, qeal& P_0, qeal& A_2, qeal& A_1, qeal& A_0, qeal& B_2, qeal& B_1, qeal& B_0, qeal& C_2, qeal& C_1, qeal& C_0, qeal& D_2, qeal& D_1, qeal& D_0, qeal S)
	{
		P_4 += (A_2 * B_2 - C_2 * D_2) * S;
		P_3 += (A_2 * B_1 + A_1 * B_2 - C_2 * D_1 - C_1 * D_2) * S;
		P_2 += (A_0 * B_2 + A_2 * B_0 + A_1 * B_1 - C_0 * D_2 - C_2 * D_0 - C_1 * D_1) * S;
		P_1 += (A_0 * B_1 + A_1 * B_0 - C_0 * D_1 - C_1 * D_0) * S;
		P_0 += (A_0 * B_0 - C_0 * D_0) * S;
	}

	__device__ __forceinline__
		void genSexticCoeffs(qeal& P_6, qeal& P_5, qeal& P_4, qeal& P_3, qeal& P_2, qeal& P_1, qeal& P_0, qeal& A_2, qeal& A_1, qeal& A_0, qeal& B_2, qeal& B_1, qeal& B_0, qeal& C_2, qeal& C_1, qeal& C_0, qeal S)
	{
		P_6 += A_2 * B_2 * C_2 * S;
		P_5 += (A_2 * B_1 * C_2 + A_1 * B_2 * C_2 + A_2 * B_2 * C_1) * S;
		P_4 += (A_2 * B_0 * C_2 + A_2 * B_2 * C_0 + A_0 * B_2 * C_2 + A_2 * B_1 * C_1 + A_1 * B_2*C_1 + A_1 * B_1 * C_2) * S;
		P_3 += (A_2 * B_1 * C_0 + A_1 * B_2 * C_0 + A_0 * B_2 * C_1 + A_2 * B_0 * C_1 + A_1 * B_1 * C_1 + A_0 * B_1*C_2 + A_1 * B_0*C_2) * S;
		P_2 += (A_2 * B_0 * C_0 + A_0 * B_2 * C_0 + A_0 * B_0 * C_2 + A_1 * B_1*C_0 + A_0 * B_1 * C_1 + A_1 * B_0 * C_1) * S;
		P_1 += (A_1 * B_0 * C_0 + A_0 * B_1 * C_0 + A_0 * B_0 * C_1) * S;
		P_0 += A_0 * B_0 * C_0 * S;
	}

	__device__ __forceinline__
		qeal SolveSexticEqForFTC(qeal& W_6, qeal& W_5, qeal& W_4, qeal& W_3, qeal& W_2, qeal& W_1, qeal& W_0, int& banRootsNum, qeal* banRoots, qeal* range, qeal norrow)
	{
		qeal ft = 1.0;
		int count = 0;
		qeal x[6];
		//range[RANGE_MIN_INDEX] = 0.0;
		//range[RANGE_MAX_INDEX] = 1.0;

		if (IS_Numerical_ZERO(W_6) && IS_Numerical_ZERO(W_5))
		{
			if (solveQuarticEquation(W_4, W_3, W_2, W_1, W_0, x, count, true, range[RANGE_MIN_INDEX], range[RANGE_MAX_INDEX]))
			{
				if (x[0] > 0 && x[0] < 1)
					ft = x[0];
			}
			else
			{
				qeal fx = computeSexticEquation(range[RANGE_MAX_INDEX], W_6, W_5, W_4, W_3, W_2, W_1, W_0);
				if (fx > 0) ft = range[RANGE_MAX_INDEX];
				else ft = NewtonSolverSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, range, norrow);
			}
		}
		else ft = NewtonSolverSexticEqForFTC(W_6, W_5, W_4, W_3, W_2, W_1, W_0, range, norrow);

		for (int j = 0; j < banRootsNum; j++)
		{
			if (IS_Numerical_ZERO(ft - banRoots[j]))
			{
				ft = banRoots[j];
				break;
			}
		}

		if (ft > 0 && ft < 1.0)
			return ft;

		return 1.0;
	}

	__device__ __forceinline__ qeal NewtonSolverSexticEqForFTC(qeal W_6, qeal W_5, qeal W_4, qeal W_3, qeal W_2, qeal W_1, qeal W_0, qeal* r, qeal norrow)
	{
		qeal f_rmin = computeSexticEquation(r[RANGE_MIN_INDEX], W_6, W_5, W_4, W_3, W_2, W_1, W_0);
		qeal f_rmax = computeSexticEquation(r[RANGE_MAX_INDEX], W_6, W_5, W_4, W_3, W_2, W_1, W_0);
		if (f_rmax <= 0) // intersection
			return NewtonSolverForIntersectionToi(W_6, W_5, W_4, W_3, W_2, W_1, W_0, r[RANGE_MIN_INDEX], r[RANGE_MAX_INDEX]);
		else
		{
			int maxIter = 20;
			int iter = 0;
			qeal has_neg = -1;
			qeal x = r[RANGE_MIN_INDEX];
			qeal fd;
			qeal fx = f_rmin;
			do
			{
				fd = 6.0 * W_6 * x * x * x * x  * x + 5.0 * W_5 * x * x * x  * x + 4.0 * W_4 * x * x  * x + 3.0 * W_3 * x  * x + 2.0 * W_2 * x + W_1;
				if (IS_CUDA_ZERO(fd))
					break;
				x = x -fx / fd;
				if (x < r[RANGE_MIN_INDEX] || x > r[RANGE_MAX_INDEX])
					x = r[RANGE_MAX_INDEX];
				fx = computeSexticEquation(x, W_6, W_5, W_4, W_3, W_2, W_1, W_0);
				if (fx < 0)
				{
					has_neg = x;
					break;
				}
				iter++;
			} while (iter < maxIter);

			if(has_neg > r[RANGE_MIN_INDEX] && has_neg < r[RANGE_MAX_INDEX]) // penetrate
				return NewtonSolverForIntersectionToi(W_6, W_5, W_4, W_3, W_2, W_1, W_0, r[RANGE_MIN_INDEX], has_neg);
			else return r[RANGE_MAX_INDEX]; // away
		}
	}

	__device__ __forceinline__ qeal NewtonSolverForIntersectionToi(qeal& W_6, qeal& W_5, qeal& W_4, qeal& W_3, qeal& W_2, qeal& W_1, qeal& W_0, qeal mini, qeal maxi)
	{
		int maxIter = 200;
		int iter = 0;
		qeal ftc = 1.0;
		//st. F(x) = f^2(x)
		qeal x = mini;
		qeal fx, fd, Fx, Fd;
		qeal f_rmin, f_rmax;
		fx = computeSexticEquation(x, W_6, W_5, W_4, W_3, W_2, W_1, W_0);
		Fx = fx * fx;
		do
		{
			fd = 6.0 * W_6 * x * x * x * x  * x + 5.0 * W_5 * x * x * x  * x + 4.0 * W_4 * x * x  * x + 3.0 * W_3 * x  * x + 2.0 * W_2 * x + W_1;
			Fd = 2.0 * fx * fd;
			if (IS_CUDA_ZERO(Fd))
				break;
			x = x - Fx / Fd;
			fx = computeSexticEquation(x, W_6, W_5, W_4, W_3, W_2, W_1, W_0);
			Fx = fx * fx;
			iter++;
		} while (CUDA_ABS(Fx) > 1e-10 && iter < maxIter);
		x = Check_CUDA_ZERO(x);
		
		if (x > mini && x <= maxi)
			ftc = x;
		else ftc = (mini == 0.0 ? maxi : mini);
		return ftc;
	}

	__device__ __forceinline__ qeal VectorDot(qeal * v1, qeal * v2)
	{
		return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
	}



}

