/*****************************************************************************
*  LPSolver, An experimental implementation of Homogeneous and Self-Dual     *
*  Algorithm for Linear Programming.                                         *
*  Author: Tian Xie (Research Center for Management Science and Information  *
*          Analytics, Shanghai University of Finance and Economics)          *
*  Credits: Fundamental implementation idea originated from COPL_LP.         *
*           (Xiong Zhang and Yinyu Ye)                                       *
*           See http://web.stanford.edu/~yyye/Col.html .                     *
******************************************************************************/

#include "LP.h"
#ifdef _MSC_VER
// NOTHING
#else
#include <sys/time.h>
#endif
#include "omp.h"

void CheckError(int ExitID, char* ErrMsg)
{
	if (ExitID)
	{
		printf("ERROR: Msg = %s, ErrorCode = %d\n", ErrMsg, ExitID);
		LinearEquation_Destruct();
#ifdef _MSC_VER
		system("pause");
#endif
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
			double Value = V_Matrix_Value[j];
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

int SPMM_FLAG[MAX_COLS][OMP_THREADS_MAX];
double SPMM_Value[MAX_COLS][OMP_THREADS_MAX];
int nnzRow[MAX_ROWS];

// Make sure that csrRow has been allocated as int[n_Row + 1]!
void ADAt_Allocate(int* nnzADAt, double** p_csrVal, int* csrRow, int** p_csrCol)
{
	memset(SPMM_FLAG, -1, sizeof(int) * n_Col * OMP_THREADS_MAX);
	#pragma omp parallel for
	for (int i = 0; i < n_Row; i ++)
	{
		int threadId = omp_get_thread_num();
		nnzRow[i] = 0;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int j = V_Matrix_Col[p];
			for (int q = p; q != -1; q = V_Matrix_Col_Next[q])
			{
				int k = V_Matrix_Row[q];
				if (SPMM_FLAG[k][threadId] != i)
				{
					SPMM_FLAG[k][threadId] = i;
					nnzRow[i] ++;
				}
			}
		}
	}
	csrRow[0] = 0;
	for (int i = 0; i < n_Row; i ++)
		csrRow[i + 1] = csrRow[i] + nnzRow[i];

	*nnzADAt = csrRow[n_Row];
	double* csrVal = (double*) malloc(sizeof(double) * csrRow[n_Row]);
	*p_csrVal = csrVal;
	int* csrCol = (int*) malloc(sizeof(int) * csrRow[n_Row]);
	*p_csrCol = csrCol;
	
	memset(SPMM_FLAG, -1, sizeof(int) * n_Col * OMP_THREADS_MAX);
	#pragma omp parallel for
	for (int i = 0; i < n_Row; i ++)
	{
		int threadId = omp_get_thread_num();
		int nnz_i = csrRow[i];
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int j = V_Matrix_Col[p];
			for (int q = p; q != -1; q = V_Matrix_Col_Next[q])
			{
				int k = V_Matrix_Row[q];
				if (SPMM_FLAG[k][threadId] != i)
				{
					SPMM_FLAG[k][threadId] = i;
					csrVal[nnz_i] = 1; // Symbolic
					csrCol[nnz_i] = k;
					nnz_i ++;
				}
			}
		}
	}
}

void ADAt_Calc(double* d, double* csrVal, int* csrRow, int* csrCol)
{
	memset(SPMM_FLAG, -1, sizeof(int) * n_Col * OMP_THREADS_MAX);
	#pragma omp parallel for
	for (int i = 0; i < n_Row; i ++)
	{
		int threadId = omp_get_thread_num();
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int j = V_Matrix_Col[p];
			for (int q = p; q != -1; q = V_Matrix_Col_Next[q])
			{
				int k = V_Matrix_Row[q];
				if (SPMM_FLAG[k][threadId] != i)
				{
					SPMM_FLAG[k][threadId] = i;
					SPMM_Value[k][threadId] = 0;
				}
				SPMM_Value[k][threadId] += V_Matrix_Value[p] * V_Matrix_Value[q] * d[j];
			}
		}
		for (int p = csrRow[i]; p < csrRow[i + 1]; p ++)
			csrVal[p] = SPMM_Value[csrCol[p]][threadId];
	}
}
