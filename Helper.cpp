/*****************************************************************************
*  LPSolver, An experimental implementation of Homogeneous and Self-Dual     *
*  Algorithm for Linear Programming.                                         *
*  Author: Tian Xie (Research Center for Management Science and Information  *
*          Analytics, Shanghai University of Finance and Economics)          *
*  Credits: (1) Fundamental implementation idea originated from COPL_LP.     *
*               (Xiong Zhang and Yinyu Ye)                                   *
*               See http://web.stanford.edu/~yyye/Col.html .                 *
*           (2) Sparse Cholesky Decomposition is supported by CHOLMOD.       *
*               (Timothy A. Davis)                                           *
******************************************************************************/

#include "LP.h"
#ifdef _MSC_VER
// NOTHING
#else
#include <sys/time.h>
#endif

void CheckError(int ExitID, char* ErrMsg)
{
	if (ExitID)
	{
		printf("ERROR: Msg = %s\n", ErrMsg);
		exit(ExitID);
	}
}

double GetTime()
{
#ifdef _MSC_VER
	return clock() / (double) CLOCKS_PER_SEC;
#else
	timeval start;
	gettimeofday(&start, NULL);
	return (double) start.tv_sec + (double) start.tv_usec * 1e-6;
#endif
}

double DotProduct(int n, double* a, double* b)
{
	double Ret = 0;
	for (int i = 0; i < n; i ++)
		Ret += a[i] * b[i];
	return Ret;
}

void SetScaledVector(int n, double alpha, double* src, double* dest) // dest = alpha * src
{
	for (int i = 0; i < n; i ++)
	{
		if (alpha == 1)
			dest[i] = src[i];
		else if (alpha == -1)
			dest[i] = -src[i];
		else if (alpha == 0)
			dest[i] = 0;
		else
			dest[i] = alpha * src[i];
	}
}

// Transpose = 0: dest = dest + Sign * A * v
// Transpose = 1: dest = dest + Sign * A^T * v
// Sign \in {1, -1}
void SetATimesVector(int Transpose, int Sign, double* v, double* dest) 
{
	for (int Col = 0; Col < n_Col; Col ++)
		for (int j = V_Matrix_Col_Head[Col]; j != -1; j = V_Matrix_Col_Next[j])
		{
			int Row = V_Matrix_Row[j];
			int Value = V_Matrix_Value[j];
			if (! Transpose) // A * v
			{
				if (Sign == 1)
					dest[Row] += Value * v[Col];
				else
					dest[Row] -= Value * v[Col];
			}
			else // A^T * v
			{
				if (Sign == 1)
					dest[Col] += Value * v[Row];
				else
					dest[Col] -= Value * v[Row];
			}
		}
}
