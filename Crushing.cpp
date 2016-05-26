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

int CRUSH_Main()
{
	printf("CRUSHING BEGIN\n");
	V_Cost_Intercept = 0.0;

	// 1. Standardize Variable Bounds. If x[i] is not free, adjust to 0 <= x[i] <= u[i].
	for (int Col = 0; Col < n_Col; Col ++)
	{
		V_Crushing_Times[Col] = 1;
		V_Crushing_Add[Col] = 0.0;

		if (V_LB[Col] > V_UB[Col])
			CheckError(1, "CRUSH_Main: Lower Bound > Upper Bound!");
		if (fabs(V_LB[Col]) < Input_Tolerance) // Lower Bound is already 0
			continue;
		if (V_LB[Col] >= -MaxFinite) // Lower Bound is Finite, use (x[i] - V_LB[i]) to instead (x[i])
		{
			V_Crushing_Add[Col] = V_LB[Col];
			for (int i = V_Matrix_Head[Col]; i != -1; i = V_Matrix_Next[i]) // Adjust Each Row
			{
				int Row = V_Matrix_Row[i];
				double Shift = V_Matrix_Value[i] * V_LB[Col];
				V_RHS[Row] -= Shift;
				V_RHS_r[Row] -= Shift;
			}
			V_Cost_Intercept += V_Cost[Col] * V_LB[Col]; // Adjust Objective
			if (V_UB[Col] <= MaxFinite)
				V_UB[Col] -= V_LB[Col];
			V_LB[Col] = 0.0;
		}
		else if (V_UB[Col] <= MaxFinite) // Lower Bound = -Inf, Upper Bound is Finite, use (-x[i] + V_UB[i]) to instead (x[i])
		{
			V_Crushing_Times[Col] = -1;
			V_Crushing_Add[Col] = V_LB[Col];
			for (int i = V_Matrix_Head[Col]; i != -1; i = V_Matrix_Next[i]) // Adjust Each Row
			{
				int Row = V_Matrix_Row[i];
				V_Matrix_Value[i] = -V_Matrix_Value[i];
				double Shift = V_Matrix_Value[i] * V_UB[Col];
				V_RHS[Row] += Shift;
				V_RHS_r[Row] += Shift;
			}
			V_Cost[Col] = -V_Cost[Col];
			V_Cost_Intercept += V_Cost[Col] * V_UB[Col]; // Adjust Objective
			V_LB[Col] = 0.0;
			V_UB[Col] = MaxPositive;
		}
		// Note: Ignore Free Variables
	}

	// 2. Standardize Rows to Equalities. A(i, :)' * x = b(i) for all i
	for (int Row = 0; Row < n_Row; Row ++)
	{
		if (Row_Type[Row] == 'E') // Already equality
			continue;
		// Create a slack variable and the corresponding column
		int Col = n_Col;
		V_Matrix_Head[Col] = -1;
		V_Cost[Col] = 0.0;
		V_LB[Col] = 0.0;
		V_UB[Col] = MaxPositive;
		n_Col ++;

		V_Matrix_Row[n_Element] = Row;
		if (Row_Type[Row] == 'L') // A(i, :)' * x + s = b(i), s >= 0
			V_Matrix_Value[n_Element] = 1.0;
		else if (Row_Type[Row] == 'G') // A(i, :)' * x - s = b(i), s >= 0
			V_Matrix_Value[n_Element] = -1.0;
		else if (Row_Type[Row] == 'R')
		{
			// r(i) <= A(i, :)' * x <= b(i)  <=>  A(i, :)' * x - s = r(i) && 0 <= s <= b(i) - r(i)
			V_Matrix_Value[n_Element] = -1.0;
			V_UB[Col] = V_RHS[Row] - V_RHS_r[Row];
			V_RHS[Row] = V_RHS_r[Row];
		}
		V_Matrix_Next[n_Element] = V_Matrix_Head[Col];
		V_Matrix_Head[Col] = n_Element;
		n_Element ++;
	}
	printf("    After Crushing, %d Rows, %d Columns, %d Elements.\n", n_Row, n_Col, n_Element);
	printf("CRUSHING END\n");
	return 0;
}
