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

int Count_LDP;
LUSOLrec* LUSOL;
double Presolve_RHS[MAX_ROWS + 1], Presolve_PI[MAX_ROWS + 1];
int Presolve_iTemp[MAX_ROWS + MAX_COLS + 2];
double Presolve_dTemp[MAX_ROWS + MAX_COLS + 2];
int Presolve_In_Base[MAX_COLS], Presolve_Cur_Base[MAX_ROWS];
int Presolve_LD_Out[MAX_ROWS]; // 1: Empty rows or Independent rows
int n_Row_Reduced, Presolve_LD_Reduced_NewToOld[MAX_ROWS], Presolve_LD_Reduced_OldToNew[MAX_ROWS];

// The procedure will be activated as soon as the basic presolving procedure is finished. 
void Presolve_Linear_Dependent_Init()
{
	Count_LDP = 0;
	// For temporary use
	memset(Row_Element_Count, 0, sizeof(int) * n_Row);
	memset(Col_Element_Count, 0, sizeof(int) * n_Col);
	for (int j = 0; j < n_Col; j ++)
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
		{
			Row_Element_Count[V_Matrix_Row[p]] ++;
			Col_Element_Count[j] ++;
		}
	memset(Presolve_LD_Out, 0, sizeof(int) * n_Row);
	// 1. Disable all empty rows
	for (int i = 0; i < n_Row; i ++)
		if (Row_Element_Count[i] == 0)
		{
			if (fabs(V_RHS[i]) > TOLPRIMAL)
			{
				printf("    LinDep: Problem is Primal Infeasible!\n");
				LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
				return;
			}
			Row_Disable[i] = 1;
			Presolve_LD_Out[i] = 1;
			Count_LDP ++;
		}
	// 2. Rows with (implied) column singletons -> Independent
	// iTemp: Queue
	int Tail = 0;
	for (int j = 0; j < n_Col; j ++)
		if (Col_Element_Count[j] == 1)
			Presolve_iTemp[Tail ++] = j;
	for (int h = 0; h < Tail; h ++)
	{
		int j = Presolve_iTemp[h];
		int i;
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
			if (! Presolve_LD_Out[V_Matrix_Row[p]])
			{
				i = V_Matrix_Row[p];
				break;
			}
		// Row i is independent
		Presolve_LD_Out[i] = 1;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int k = V_Matrix_Col[p];
			Col_Element_Count[k] --;
			if (Col_Element_Count[k] == 1) // New column singleton found!
				Presolve_iTemp[Tail ++] = k;
		}
	}
	// 3. Arrange new model
	n_Row_Reduced = 0;
	for (int i = 0; i < n_Row; i ++)
		if (! Presolve_LD_Out[i])
		{
			Presolve_LD_Reduced_NewToOld[n_Row_Reduced] = i;
			Presolve_LD_Reduced_OldToNew[i] = n_Row_Reduced;
			n_Row_Reduced ++;
		}
		else
			Presolve_LD_Reduced_OldToNew[i] = -1;
	printf("    LU Elimination: %d (out of %d) Rows in LU Decomposition.\n", n_Row_Reduced, n_Row);
}

int Presolve_InitBase_Order[MAX_ROWS];

void Presolve_Linear_Dependent_FindInitBase()
{
	for (int j = 0; j < n_Col; j ++)
		Presolve_In_Base[j] = 0;
	for (int i = 0; i < n_Row_Reduced; i ++)
		Presolve_Cur_Base[i] = -1; // Artificial Index

	// Maintain row density, from high to low
	for (int i = 0; i < n_Row_Reduced; i ++)
		Presolve_InitBase_Order[i] = i;
	for (int i = 0; i < n_Row_Reduced; i ++)
		for (int j = 0; j < n_Row_Reduced; j ++)
			if (Row_Element_Count[Presolve_LD_Reduced_NewToOld[Presolve_InitBase_Order[i]]] > Row_Element_Count[Presolve_LD_Reduced_NewToOld[Presolve_InitBase_Order[j]]])
				swap(Presolve_InitBase_Order[i], Presolve_InitBase_Order[j]);

	int Init_Structured_Cols = 0, Init_Artificial_Cols = 0;
	// Note that Presolve_LD_Out is contaminated after the following code, use Presolve_LD_Reduced_OldToNew instead!
	memset(Presolve_LD_Out, 0, sizeof(int) * n_Row_Reduced);
	for (int r = 0; r < n_Row_Reduced; r ++)
	{
		// Find the unassigned row with highest density
		int rm = Presolve_InitBase_Order[r];
		if (Presolve_LD_Out[rm])
			continue;
		// Temporarily remove Row rm, push some linear independent rows into the initial basis
		Presolve_LD_Out[rm] = 1;
		Init_Artificial_Cols ++;
		int Tail = 0;
		for (int p = V_Matrix_Row_Head[Presolve_LD_Reduced_NewToOld[rm]]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int k = V_Matrix_Col[p];
			Col_Element_Count[k] --;
			if (Col_Element_Count[k] == 1) // New column singleton found!
				Presolve_iTemp[Tail ++] = k;
		}
		for (int h = 0; h < Tail; h ++)
		{
			int j = Presolve_iTemp[h];
			int i;
			for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
				if (Presolve_LD_Reduced_OldToNew[V_Matrix_Row[p]] != -1 && ! Presolve_LD_Out[Presolve_LD_Reduced_OldToNew[V_Matrix_Row[p]]])
				{
					i = Presolve_LD_Reduced_OldToNew[V_Matrix_Row[p]];
					break;
				}
			if (Presolve_LD_Out[i])
				continue;
			// Row i is in the initial basis for column j
			Presolve_LD_Out[i] = 1;
			Presolve_In_Base[j] = 1;
			Presolve_Cur_Base[i] = j;
			Init_Structured_Cols ++;
			for (int p = V_Matrix_Row_Head[Presolve_LD_Reduced_NewToOld[i]]; p != -1; p = V_Matrix_Row_Next[p])
			{
				int k = V_Matrix_Col[p];
				Col_Element_Count[k] --;
				if (Col_Element_Count[k] == 1) // New column singleton found!
					Presolve_iTemp[Tail ++] = k;
			}
		}
	}
	printf("    Initial Basis: %d Structured Column(s), %d Artificial Column(s)\n", Init_Structured_Cols, Init_Artificial_Cols);
}

void Presolve_Linear_Dependent_InitLU()
{
	LUSOL = LUSOL_create(stdout, 0, LUSOL_PIVMOD_TPP, 0);
	LUSOL -> luparm[LUSOL_IP_USEROWL0] = LUSOL_OTHERORDER;
	LUSOL -> luparm[LUSOL_IP_SCALAR_NZA] = LUSOL_MULT_nz_a;
	LUSOL -> luparm[LUSOL_IP_KEEPLU] = TRUE;
	LUSOL_sizeto(LUSOL, n_Row_Reduced, n_Col, n_Element * LUSOL_MULT_nz_a);

	// Calculate Col_Element_Count for reduced model
	memset(Col_Element_Count, 0, sizeof(int) * n_Col);
	for (int j = 0; j < n_Col; j ++)
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
			Col_Element_Count[j] += (Presolve_LD_Reduced_OldToNew[V_Matrix_Row[p]] != -1);
}

void Presolve_Linear_Dependent_Destruct()
{
	LUSOL_free(LUSOL);
}

void Presolve_Linear_Dependent_Renew_LU()
{
	int Size = 1;
	for (int i = 0; i < n_Row_Reduced; i ++)
		Size += (Presolve_Cur_Base[i] == -1) ? 1 : Col_Element_Count[Presolve_Cur_Base[i]];
	//printf("    Renewing LU: Size = %d\n", Size);
	double* Bij = (double*) calloc(Size, sizeof(double));
	int* iB = (int*) calloc(Size, sizeof(int));
	int* jB = (int*) calloc(Size, sizeof(int));
	Size = 1;
	for (int i = 0; i < n_Row_Reduced; i ++)
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
			if (Presolve_LD_Reduced_OldToNew[V_Matrix_Row[p]] == -1)
				continue;
			Bij[Size] = V_Matrix_Value[p];
			iB[Size] = Presolve_LD_Reduced_OldToNew[V_Matrix_Row[p]] + 1;
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
	int Count_Replace = 0;
	// PI^T A = e_i^T B^{-1} A: The i-th row of (B^{-1} A)
	for (int i = 0; i < n_Row_Reduced; i ++) // The i-th artificial variable
	{
		if (Presolve_Cur_Base[i] != -1)
			continue;
		// Solve B^T PI = e_i
		// As the input parameter RHS, Presolve_PI will be overwritten
		memset(Presolve_PI, 0, sizeof(double) * (n_Row_Reduced + 1));
		Presolve_PI[i + 1] = 1; // 1-based e_i
		int LUSOL_inform = LUSOL_btran(LUSOL, Presolve_PI, NULL);
		if (LUSOL_inform > LUSOL_INFORM_SERIOUS)
			CheckError(LUSOL_inform, LUSOL_informstr(LUSOL, LUSOL_inform));
		
		memset(Presolve_dTemp, 0, sizeof(double) * n_Col);
		SetATimesVector(1, 1, Presolve_PI + 1, Presolve_dTemp, Presolve_LD_Reduced_OldToNew); // Sparse dgemv
		int j;
		for (j = 0; j < n_Col; j ++)
			if ((Presolve_dTemp[j] >= TOLAPIV || Presolve_dTemp[j] <= -TOLAPIV) && ! Presolve_In_Base[j])
				break;
		if (j == n_Col) // Linearly Dependent
		{
			Count_LDP ++;
			Row_Disable[Presolve_LD_Reduced_NewToOld[i]] = 1;
			/*
			if (Row_Element_Count[i] != 0)
				printf("    Non-empty Row %d (1-indexed) is linearly dependent! Row_Element_Count = %d\n", i + 1, Row_Element_Count[i]);
			*/
			for (int k = 0; k < n_Row_Reduced; k ++)
				Presolve_dTemp[k + 1] = V_RHS[Presolve_LD_Reduced_NewToOld[k]];
			int LUSOL_inform = LUSOL_ftran(LUSOL, Presolve_dTemp, NULL, FALSE);
			if (LUSOL_inform > LUSOL_INFORM_SERIOUS)
				CheckError(LUSOL_inform, LUSOL_informstr(LUSOL, LUSOL_inform));
			double DOWN = 1.0;
			for (int k = 1; k <= n_Row_Reduced; k ++)
				DOWN += Presolve_dTemp[k] * Presolve_dTemp[k];
			if (fabs(Presolve_dTemp[i + 1]) > TOLPRIMAL * sqrt(DOWN)) // Infeasible
			{
				printf("    LinDep: Problem is Primal Infeasible!\n");
				LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
				return;
			}
		}
		else
		{
			if (Row_Element_Count[Presolve_LD_Reduced_NewToOld[i]] == 0)
				printf("    Error: Empty Row %d but not linearly dependent!\n", i + 1);
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
					Presolve_dTemp[Presolve_LD_Reduced_OldToNew[V_Matrix_Row[p]] + 1] = V_Matrix_Value[p];
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
	if (LP_Status != LP_STATUS_OK)
		return;
	if (n_Row_Reduced == 0)
	{
		printf("    Empty LU.\n");
		return;
	}
	Presolve_Linear_Dependent_FindInitBase();
	Presolve_Linear_Dependent_InitLU();
	Presolve_Linear_Dependent_Renew_LU();
	Presolve_Linear_Dependent_Solve();
	Presolve_Linear_Dependent_Destruct();
}