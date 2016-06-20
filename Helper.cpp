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
//#include "immintrin.h"

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
	if (alpha == 1)
	{
		for (int i = 0; i < n; i ++)
			dest[i] = src[i];
	}
	else if (alpha == -1)
	{
		for (int i = 0; i < n; i ++)
			dest[i] = -src[i];
	}
	else if (alpha == 0)
	{
		for (int i = 0; i < n; i ++)
			dest[i] = 0;
	}
	else
	{
		for (int i = 0; i < n; i ++)
			dest[i] = alpha * src[i];
	}
}

// Transpose = 0: dest = dest + Sign * A * v
// Transpose = 1: dest = dest + Sign * A^T * v
// Sign \in {1, -1}
void SetATimesVector(int Transpose, int Sign, double* v, double* dest) 
{
	if (! Transpose) // A * v
	{
		if (Sign == 1)
		{
			for (int Col = 0; Col < n_Col; Col ++)
			{
				double vCol = v[Col];
				for (int j = V_Matrix_Col_Head[Col]; j != -1; j = V_Matrix_Col_Next[j])
					dest[V_Matrix_Row[j]] += V_Matrix_Value[j] * vCol;
			}
		}
		else
		{
			for (int Col = 0; Col < n_Col; Col ++)
			{
				double vCol = v[Col];
				for (int j = V_Matrix_Col_Head[Col]; j != -1; j = V_Matrix_Col_Next[j])
					dest[V_Matrix_Row[j]] -= V_Matrix_Value[j] * vCol;
			}
		}
	}
	else // A^T * v
	{
		if (Sign == 1)
		{
			for (int Col = 0; Col < n_Col; Col ++)
			{
				double destCol = dest[Col];
				for (int j = V_Matrix_Col_Head[Col]; j != -1; j = V_Matrix_Col_Next[j])
					destCol += V_Matrix_Value[j] * v[V_Matrix_Row[j]];
				dest[Col] = destCol;
			}
		}
		else
		{
			for (int Col = 0; Col < n_Col; Col ++)
			{
				double destCol = dest[Col];
				for (int j = V_Matrix_Col_Head[Col]; j != -1; j = V_Matrix_Col_Next[j])
					destCol -= V_Matrix_Value[j] * v[V_Matrix_Row[j]];
				dest[Col] = destCol;
			}
		}
	}
}

int SPMM_FLAG[OMP_THREADS_MAX][MAX_COLS];
double SPMM_Value[OMP_THREADS_MAX][MAX_COLS];
int nnzRow[MAX_ROWS];

// Make sure that csrRow has been allocated as int[n_Row + 1]!
void ADAt_Allocate(int* nnzADAt, double** p_csrVal, int* csrRow, int** p_csrCol)
{
	for (int i = 0; i < OMP_THREADS_MAX; i ++)
		memset(SPMM_FLAG[i], -1, sizeof(int) * n_Col);
	#pragma omp parallel for
	for (int i = 0; i < n_Row; i ++)
	{
		int* p_Flag = SPMM_FLAG[omp_get_thread_num()];
		nnzRow[i] = 0;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int j = V_Matrix_Col[p];
			for (int q = p; q != -1; q = V_Matrix_Col_Next[q])
			{
				int k = V_Matrix_Row[q];
				if (p_Flag[k] != i)
				{
					p_Flag[k] = i;
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
	
	for (int i = 0; i < OMP_THREADS_MAX; i ++)
		memset(SPMM_FLAG[i], -1, sizeof(int) * n_Col);
	#pragma omp parallel for
	for (int i = 0; i < n_Row; i ++)
	{
		int* p_Flag = SPMM_FLAG[omp_get_thread_num()];
		int nnz_i = csrRow[i];
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int j = V_Matrix_Col[p];
			for (int q = p; q != -1; q = V_Matrix_Col_Next[q])
			{
				int k = V_Matrix_Row[q];
				if (p_Flag[k] != i)
				{
					p_Flag[k] = i;
					csrVal[nnz_i] = 1; // Symbolic
					csrCol[nnz_i] = k;
					nnz_i ++;
				}
			}
		}
	}
}

void ADAt_Allocate(int* nnzADAt, double** p_csrVal, int* csrRow, int** p_csrCol, int* LinkerTocsrAt, int* csrRowAt, int* csrColAt)
{
	for (int i = 0; i < OMP_THREADS_MAX; i ++)
		memset(SPMM_FLAG[i], -1, sizeof(int) * n_Col);
	#pragma omp parallel for
	for (int i = 0; i < n_Row; i ++)
	{
		int* p_Flag = SPMM_FLAG[omp_get_thread_num()];
		nnzRow[i] = 0;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int j = V_Matrix_Col[p];
			for (int q = LinkerTocsrAt[p]; q < csrRowAt[j + 1]; q ++)
			{
				int k = csrColAt[q];
				if (p_Flag[k] != i)
				{
					p_Flag[k] = i;
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
	
	for (int i = 0; i < OMP_THREADS_MAX; i ++)
		memset(SPMM_FLAG[i], -1, sizeof(int) * n_Col);
	#pragma omp parallel for
	for (int i = 0; i < n_Row; i ++)
	{
		int* p_Flag = SPMM_FLAG[omp_get_thread_num()];
		int nnz_i = csrRow[i];
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int j = V_Matrix_Col[p];
			for (int q = LinkerTocsrAt[p]; q < csrRowAt[j + 1]; q ++)
			{
				int k = csrColAt[q];
				if (p_Flag[k] != i)
				{
					p_Flag[k] = i;
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
	#pragma omp parallel for
	for (int i = 0; i < n_Row; i ++)
	{
		double* p_Value = SPMM_Value[omp_get_thread_num()];
		for (int p = csrRow[i]; p < csrRow[i + 1]; p ++)
			p_Value[csrCol[p]] = 0;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			double Vp_Dj = V_Matrix_Value[p] * d[V_Matrix_Col[p]];
			for (int q = p; q != -1; q = V_Matrix_Col_Next[q])
				p_Value[V_Matrix_Row[q]] += Vp_Dj * V_Matrix_Value[q];
		}
		for (int p = csrRow[i]; p < csrRow[i + 1]; p ++)
			csrVal[p] = p_Value[csrCol[p]];
	}
}

void ADAt_Calc(double* d, double* csrVal, int* csrRow, int* csrCol, int* LinkerTocsrAt, double* csrValAt, int* csrRowAt, int* csrColAt)
{
	#pragma omp parallel for
	for (int i = 0; i < n_Row; i ++)
	{
		double* p_Value = SPMM_Value[omp_get_thread_num()];
		for (int p = csrRow[i]; p < csrRow[i + 1]; p ++)
			p_Value[csrCol[p]] = 0;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			double Vp_Dj = V_Matrix_Value[p] * d[V_Matrix_Col[p]];
			int q_end = csrRowAt[V_Matrix_Col[p] + 1];
			int q = LinkerTocsrAt[p];
			for (; q + 3 < q_end; q += 4)
			{
				p_Value[csrColAt[q]] += Vp_Dj * csrValAt[q];
				p_Value[csrColAt[q + 1]] += Vp_Dj * csrValAt[q + 1];
				p_Value[csrColAt[q + 2]] += Vp_Dj * csrValAt[q + 2];
				p_Value[csrColAt[q + 3]] += Vp_Dj * csrValAt[q + 3];
			}
			if (q < q_end) p_Value[csrColAt[q]] += Vp_Dj * csrValAt[q];
			if (q + 1 < q_end) p_Value[csrColAt[q + 1]] += Vp_Dj * csrValAt[q + 1];
			if (q + 2 < q_end) p_Value[csrColAt[q + 2]] += Vp_Dj * csrValAt[q + 2];
		}
		for (int p = csrRow[i]; p < csrRow[i + 1]; p ++)
			csrVal[p] = p_Value[csrCol[p]];
	}
}

/*
void ADAt_Calc_FMA(double* d, double* csrVal, int* csrRow, int* csrCol, int* LinkerTocsrAt, double* csrValAt, int* csrRowAt, int* csrColAt)
{
	//#pragma omp parallel for
	for (int i = 0; i < n_Row; i ++)
	{
		double* v_Value = SPMM_Value[omp_get_thread_num()];
		ATTR_ALIGN(32) double Scalar[4], Mult[4], Res[4], Add[4];
		for (int p = csrRow[i]; p < csrRow[i + 1]; p ++)
			v_Value[csrCol[p]] = 0;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			double Vp_Dj = V_Matrix_Value[p] * d[V_Matrix_Col[p]];
			int q_end = csrRowAt[V_Matrix_Col[p] + 1];
			int q = LinkerTocsrAt[p];
			Scalar[0] = Scalar[1] = Scalar[2] = Scalar[3] = Vp_Dj;
			__m256d v_Scalar = _mm256_load_pd(Scalar);
			for (; q + 3 < q_end; q += 4)
			{
				Add[0] = v_Value[csrColAt[q]];
				Add[1] = v_Value[csrColAt[q + 1]];
				Add[2] = v_Value[csrColAt[q + 2]];
				Add[3] = v_Value[csrColAt[q + 3]];
				Mult[0] = csrValAt[q];
				Mult[1] = csrValAt[q + 1];
				Mult[2] = csrValAt[q + 2];
				Mult[3] = csrValAt[q + 3];
				__m256d v_Add = _mm256_load_pd(Add);
				__m256d v_Mult = _mm256_load_pd(Mult);
				__m256d v_Res = _mm256_fmadd_pd(v_Scalar, v_Mult, v_Add);
				_mm256_store_pd(Res, v_Res);
				v_Value[csrColAt[q]] = Res[0];
				v_Value[csrColAt[q + 1]] = Res[1];
				v_Value[csrColAt[q + 2]] = Res[2];
				v_Value[csrColAt[q + 3]] = Res[3];
			}
			if (q < q_end) v_Value[csrColAt[q]] += Vp_Dj * csrValAt[q];
			if (q + 1 < q_end) v_Value[csrColAt[q + 1]] += Vp_Dj * csrValAt[q + 1];
			if (q + 2 < q_end) v_Value[csrColAt[q + 2]] += Vp_Dj * csrValAt[q + 2];
		}
		for (int p = csrRow[i]; p < csrRow[i + 1]; p ++)
			csrVal[p] = v_Value[csrCol[p]];
	}
}
*/