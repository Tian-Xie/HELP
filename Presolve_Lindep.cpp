/******************************************************************************
 *  HELP (Homogeneous and self-dual Efficient Linear Programming)             *
 *                                                                            *
 *  An experimental implementation of Homogeneous and Self-Dual Algorithm     *
 *  for Linear Programming.                                                   *
 *                                                                            *
 *  Author: Tian Xie (Research Center for Management Science and Information  *
 *          Analytics, Shanghai University of Finance and Economics)          *
 *                                                                            *
 *  Credits: Fundamental implementation idea originated from COPL_LP.         *
 *           (Xiong Zhang and Yinyu Ye)                                       *
 *           See http://web.stanford.edu/~yyye/Col.html .                     *
 ******************************************************************************/

// The following program is mainly based on: 
// Erling D. Andersen. Finding all linearly dependent rows in large-scale linear programming, Optimization Methods and Software 6:3 (1995), 219-227.

#include "LP.h"
#include "Presolve.h"
#include "lusol/lusol.h"
#include "lusol/commonlib.h"
#include "lusol/myblas.h"

int Real_n_Row;
int max_nnz_B;
double Presolve_Vij[MAX_ROWS + 1];
int Presolve_iV[MAX_ROWS + 1];
double Presolve_RHS[MAX_ROWS + 1], Presolve_PI[MAX_ROWS + 1];
LUSOLrec* LUSOL;
int Presolve_iTemp[MAX_ROWS + MAX_COLS + 2];
double Presolve_dTemp[MAX_ROWS + MAX_COLS + 2];
int Presolve_In_Base[MAX_COLS], Presolve_Cur_Base[MAX_ROWS];

// The procedure will be activated as soon as the basic presolving procedure is finished. 
void Presolve_Linear_Dependent_Init()
{
	// The post-presolving real row count. 
	memset(Row_Element_Count, 0, sizeof(int) * n_Row);
	memset(Col_Element_Count, 0, sizeof(int) * n_Col);
	for (int j = 0; j < n_Col; j ++)
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
		{
			Row_Element_Count[V_Matrix_Row[p]] ++;
			Col_Element_Count[j] ++;
		}

	Real_n_Row = 0;
	for (int i = 0; i < n_Row; i ++)
		if (Row_Element_Count[i] > 0)
			Real_n_Row ++;

	for (int j = 0; j < n_Col; j ++)
		Presolve_In_Base[j] = 0;
	for (int i = 0; i < n_Row; i ++)
		Presolve_Cur_Base[i] = -1; // Artificial Index
	LUSOL = LUSOL_create(stdout, 0, LUSOL_PIVMOD_TPP, 0);
	LUSOL -> luparm[LUSOL_IP_USEROWL0] = LUSOL_OTHERORDER;
	LUSOL -> luparm[LUSOL_IP_SCALAR_NZA] = LUSOL_MULT_nz_a;
	LUSOL -> luparm[LUSOL_IP_KEEPLU] = TRUE;
	LUSOL_sizeto(LUSOL, MAX_ROWS, MAX_COLS, MAX_ELEMENTS * LUSOL_MULT_nz_a);
	// Initially, B = V, which is an identity matrix
	for (int i = 1; i <= n_Row; i ++)
	{
		Presolve_Vij[i] = 1;
		Presolve_iV[i] = i;
	}
	if (! LUSOL_assign(LUSOL, Presolve_iV, Presolve_iV, Presolve_Vij, n_Row, TRUE))
		CheckError(1, "LUSOL failed due to insufficient memory.");
	int LUSOL_inform = LUSOL_factorize(LUSOL);
	if (LUSOL_inform > LUSOL_INFORM_SERIOUS)
		CheckError(LUSOL_inform, LUSOL_informstr(LUSOL, LUSOL_inform));
}

void Presolve_Linear_Dependent_Destruct()
{
	LUSOL_free(LUSOL);
}

void Presolve_Linear_Dependent_Renew_LU()
{
	int Size = 1;
	for (int i = 0; i < n_Row; i ++)
		Size += (Presolve_Cur_Base[i] == -1) ? 1 : Col_Element_Count[Presolve_Cur_Base[i]];
	//printf("Renewing LU: Size = %d\n", Size);
	double* Bij = (double*) calloc(Size, sizeof(double));
	int* iB = (int*) calloc(Size, sizeof(int));
	int* jB = (int*) calloc(Size, sizeof(int));
	Size = 1;
	for (int i = 0; i < n_Row; i ++)
	{
		if (Presolve_Cur_Base[i] == -1)
		{
			Bij[Size] = 1;
			iB[Size] = jB[Size] = i + 1;
			Size ++;
			continue;
		}
		for (int p = V_Matrix_Col_Head[Presolve_Cur_Base[i]]; p != -1; p = V_Matrix_Col_Next[p])
		{
			Bij[Size] = V_Matrix_Value[p];
			iB[Size] = V_Matrix_Row[p] + 1;
			jB[Size] = i + 1;
			Size ++;
		}
	}
	if (! LUSOL_assign(LUSOL, iB, jB, Bij, Size - 1, TRUE))
		CheckError(1, "LUSOL failed due to insufficient memory.");
	int LUSOL_inform = LUSOL_factorize(LUSOL);
	if (LUSOL_inform > LUSOL_INFORM_SERIOUS)
		CheckError(LUSOL_inform, LUSOL_informstr(LUSOL, LUSOL_inform));
	free(Bij);
	free(iB);
	free(jB);
}

/* 
	Assume a set of artificial columns is appended to A such that A = [I, hatA] where hatA is the original LP constraint matrix. 
	Let the set V contain the indices of the artificial columns. 
	1. Let B = V
	2. For each i \in B \cup V
	3.     PI^T = e_i^T B^{-1}
	4.     If Exist j: |PI^T * Col(A, j)| > 0 and j \not \in B \cup V
	5.         B = B + {j} - {i}
	6.     Else
	7.         The i-th row is dependent.
	8. End For
*/

void Presolve_Linear_Dependent_Solve()
{
	Real_n_Row = 0;
	for (int i = 0; i < n_Row; i ++)
		if (Row_Element_Count[i] > 0)
			Real_n_Row ++;

	for (int j = 0; j < n_Col; j ++)
		Presolve_In_Base[j] = 0;

	memset(Presolve_dTemp, 0, sizeof(double) * n_Row);
	for (int j = 0; j < n_Col; j ++)
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
			Presolve_dTemp[V_Matrix_Row[p]] += fabs(V_Matrix_Value[p]);

	int Count_LDP = 0;
	int Count_Replace = 0;
	// PI^T A = e_i^T B^{-1} A: The i-th row of (B^{-1} A)
	for (int i = 0; i < n_Row; i ++) // The i-th artificial variable
	{
		/*
		if (i % 100 == 0)
			printf("LDP Row %d\n", i);
		*/
		// Solve B^T PI = e_i
		// As the input parameter RHS, Presolve_PI will be overwritten
		memset(Presolve_PI, 0, sizeof(double) * (n_Row + 1));
		Presolve_PI[i + 1] = 1; // 1-based e_i
		int LUSOL_inform = LUSOL_btran(LUSOL, Presolve_PI, NULL);
		if (LUSOL_inform > LUSOL_INFORM_SERIOUS)
			CheckError(LUSOL_inform, LUSOL_informstr(LUSOL, LUSOL_inform));
		
		memset(Presolve_dTemp, 0, sizeof(double) * n_Col);
		SetATimesVector(1, 1, Presolve_PI + 1, Presolve_dTemp); // Sparse dgemv
		int j;
		for (j = 0; j < n_Col; j ++)
			if ((Presolve_dTemp[j] >= TOLAPIV || Presolve_dTemp[j] <= -TOLAPIV) && ! Presolve_In_Base[j])
				break;
		if (j == n_Col) // Linearly Dependent
		{
			Count_LDP ++;
			Row_Disable[i] = 1;
			/*
			if (Row_Element_Count[i] != 0)
				printf("Non-empty Row %d (1-indexed) is linearly dependent! Row_Element_Count = %d\n", i + 1, Row_Element_Count[i]);
			*/
			for (int k = 0; k < n_Row; k ++)
				Presolve_dTemp[k + 1] = V_RHS[k];
			int LUSOL_inform = LUSOL_ftran(LUSOL, Presolve_dTemp, NULL, FALSE);
			if (LUSOL_inform > LUSOL_INFORM_SERIOUS)
				CheckError(LUSOL_inform, LUSOL_informstr(LUSOL, LUSOL_inform));
			double DOWN = 1.0;
			for (int k = 1; k <= n_Row; k ++)
				DOWN += Presolve_dTemp[k] * Presolve_dTemp[k];
			if (fabs(Presolve_dTemp[i + 1]) > TOLPRIMAL * sqrt(DOWN)) // Infeasible
			{
				printf("LinDep: Problem is Primal Infeasible!\n");
				LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
				return;
			}
		}
		else
		{
			if (Row_Element_Count[i] == 0)
				printf("Error: Empty Row %d but not linearly dependent!\n", i + 1);
			Presolve_In_Base[j] = 1;
			Presolve_Cur_Base[i] = j;
			Count_Replace ++;
			if (Count_Replace >= STEPS_MODIFICATION)
			{
				Presolve_Linear_Dependent_Renew_LU();
				Count_Replace = 0;
			}
			else
			{
				memset(Presolve_dTemp, 0, sizeof(double) * (n_Row + 1));
				for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
					Presolve_dTemp[V_Matrix_Row[p] + 1] = V_Matrix_Value[p];
				LUSOL_inform = LUSOL_replaceColumn(LUSOL, i + 1, Presolve_dTemp);
				if (LUSOL_inform > LUSOL_INFORM_SERIOUS)
					CheckError(LUSOL_inform, LUSOL_informstr(LUSOL, LUSOL_inform));
				if (LUSOL_inform == LUSOL_INFORM_RANKLOSS)
				{
					Presolve_Linear_Dependent_Renew_LU();
					Count_Replace = 0;
				}
			}
		}
	}
	printf("    %d Linearly Dependent Row(s) Detected.\n", Count_LDP);
}

void Presolve_Linear_Dependent_Main()
{
	Presolve_Linear_Dependent_Init();
	Presolve_Linear_Dependent_Solve();
	Presolve_Linear_Dependent_Destruct();
}