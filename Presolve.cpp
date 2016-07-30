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

// This program is mainly based on: 
// Erling D. Andersen and Knud D. Andersen, Presolving in linear programming, Mathematical Programming 71 (1995) 221-245.

#include "LP.h"
#include "Presolve.h"

int Presolve_Modified;
double Row_1Norm[MAX_ROWS], Col_1Norm[MAX_COLS];
int Row_Disable[MAX_ROWS], Col_Disable[MAX_COLS];
int Row_Element_Count[MAX_ROWS], Col_Element_Count[MAX_COLS];
int Presolve_Linked_List_Head, Presolve_Linked_List_Tail, Presolve_Linked_List_Next[MAX_ROWS + MAX_COLS];

// Need some sorting in case of checking duplicate rows / cols

bool Presolve_Row_Sort_Cmp(const int& p1, const int& p2)
{
	return V_Matrix_Col[p1] < V_Matrix_Col[p2];
}

bool Presolve_Col_Sort_Cmp(const int& p1, const int& p2)
{
	return V_Matrix_Row[p1] < V_Matrix_Row[p2];
}

int Presolve_Sort_Tmp[MAX_ROWS + MAX_COLS];

void Presolve_Row_Sort(int i)
{
	if (Row_Disable[i] || Row_Element_Count[i] == 0)
		return;
	int Cnt = 0;
	for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		Presolve_Sort_Tmp[Cnt ++] = p;
	sort(Presolve_Sort_Tmp, Presolve_Sort_Tmp + Cnt, Presolve_Row_Sort_Cmp);

	V_Matrix_Row_Head[i] = Presolve_Sort_Tmp[0];
	for (int t = 0; t < Cnt; t ++)
	{
		V_Matrix_Row_Next[Presolve_Sort_Tmp[t]] = (t == Cnt - 1) ? -1 : Presolve_Sort_Tmp[t + 1];
		V_Matrix_Row_Prev[Presolve_Sort_Tmp[t]] = (t == 0) ? -1 : Presolve_Sort_Tmp[t - 1];
	}
}

void Presolve_Col_Sort(int j)
{
	if (Col_Disable[j] || Col_Element_Count[j] == 0)
		return;
	int Cnt = 0;
	for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
		Presolve_Sort_Tmp[Cnt ++] = p;
	sort(Presolve_Sort_Tmp, Presolve_Sort_Tmp + Cnt, Presolve_Col_Sort_Cmp);

	V_Matrix_Col_Head[j] = Presolve_Sort_Tmp[0];
	for (int t = 0; t < Cnt; t ++)
	{
		V_Matrix_Col_Next[Presolve_Sort_Tmp[t]] = (t == Cnt - 1) ? -1 : Presolve_Sort_Tmp[t + 1];
		V_Matrix_Col_Prev[Presolve_Sort_Tmp[t]] = (t == 0) ? -1 : Presolve_Sort_Tmp[t - 1];
	}
}

void Presolve_Init()
{
	for (int i = 0; i < n_Row; i ++)
	{
		Row_Disable[i] = Row_Element_Count[i] = 0;
		Row_1Norm[i] = 0.0;
	}
	for (int j = 0; j < n_Col; j ++)
	{
		Col_Disable[j] = Col_Element_Count[j] = 0;
		Col_1Norm[j] = 0.0;
		V_Presolve_Linear_Replace[j] = -1;
	}
	V_Presolve_Duplicate_Column_Cnt = 0;
	// Build Horizontal Links
	for (int i = 0; i < n_Row; i ++)
		V_Matrix_Row_Head[i] = -1;
	for (int j = 0; j < n_Col; j ++)
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
		{
			int i = V_Matrix_Row[p];
			V_Matrix_Col[p] = j;
			V_Matrix_Row_Next[p] = V_Matrix_Row_Head[i];
			V_Matrix_Row_Head[i] = p;
			// Count Elements
			Row_Element_Count[i] ++;
			Col_Element_Count[j] ++;
			Row_1Norm[i] += fabs(V_Matrix_Value[p]);
			Col_1Norm[j] += fabs(V_Matrix_Value[p]);
		}
	// Need some sorting, and build reverse links

	for (int i = 0; i < n_Row; i ++)
		Presolve_Row_Sort(i);
	for (int j = 0; j < n_Col; j ++)
		Presolve_Col_Sort(j);
}

void Presolve_DEBUG()
{
	int cnt1 = 0, cnt2 = 0, cnt3 = 0, cnt4 = 0;
	for (int i = 0; i < n_Row; i ++)
	{
		if (Row_Disable[i])
			continue;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
			cnt1 ++;
		cnt3 += Row_Element_Count[i];
	}
	for (int i = 0; i < n_Col; i ++)
	{
		if (Col_Disable[i])
			continue;
		for (int p = V_Matrix_Col_Head[i]; p != -1; p = V_Matrix_Col_Next[p])
			cnt2 ++;
		cnt4 += Col_Element_Count[i];
	}
	printf("cnt1 = %d, cnt2 = %d, cnt3 = %d, cnt4 = %d\n", cnt1, cnt2, cnt3, cnt4);
}

void Presolve_Delete_Row(int i)
{
	for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
	{
		if (V_Matrix_Col_Prev[p] != -1)
			V_Matrix_Col_Next[V_Matrix_Col_Prev[p]] = V_Matrix_Col_Next[p];
		else
			V_Matrix_Col_Head[V_Matrix_Col[p]] = V_Matrix_Col_Next[p];
		if (V_Matrix_Col_Next[p] != -1)
			V_Matrix_Col_Prev[V_Matrix_Col_Next[p]] = V_Matrix_Col_Prev[p];
		Col_Element_Count[V_Matrix_Col[p]] --;
		Col_1Norm[V_Matrix_Col[p]] -= fabs(V_Matrix_Value[p]);
	}
	V_Matrix_Row_Head[i] = -1;
	Row_Element_Count[i] = 0; // Delete the whole row, but not disabled. 
	Presolve_Modified ++;
}

void Presolve_Delete_Col(int i)
{
	for (int p = V_Matrix_Col_Head[i]; p != -1; p = V_Matrix_Col_Next[p])
	{
		if (V_Matrix_Row_Prev[p] != -1)
			V_Matrix_Row_Next[V_Matrix_Row_Prev[p]] = V_Matrix_Row_Next[p];
		else
			V_Matrix_Row_Head[V_Matrix_Row[p]] = V_Matrix_Row_Next[p];
		if (V_Matrix_Row_Next[p] != -1)
			V_Matrix_Row_Prev[V_Matrix_Row_Next[p]] = V_Matrix_Row_Prev[p];
		Row_Element_Count[V_Matrix_Row[p]] --;
		Row_1Norm[V_Matrix_Row[p]] -= fabs(V_Matrix_Value[p]);
	}
	V_Matrix_Col_Head[i] = -1;
	Col_Element_Count[i] = 0; // Delete the whole column, but not disabled. 
	Presolve_Modified ++;
}

// Fix x_j = Val. Essentially set x_j = 0 and then linear transformation.
void Presolve_Fix_Variable(int j, double Val)
{
	for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
	{
		int i = V_Matrix_Row[p];
		V_RHS[i] -= V_Matrix_Value[p] * Val;
	}
	V_Cost_Intercept += V_Cost[j] * Val;

	V_Crushing_Add[j] += V_Crushing_Times[j] * Val;
	V_Crushing_Times[j] = 0;
	V_LB[j] = V_UB[j] = 0;
	Presolve_Delete_Col(j);
	Col_Disable[j] = 1;

	Presolve_Modified ++;
}

void Presolve_Set_Variable_LB(int j, double lb)
{
	if (lb < -MaxFinite || lb <= V_LB[j])
		return;
	// Set x_j >= lb, like crushing
	// If x_j is currently free, set (x_j - lb) >= 0; otherwise x_j >= 0, also set (x_j - lb) >= 0
	V_Crushing_Add[j] = lb * V_Crushing_Times[j];
	for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
	{
		int i = V_Matrix_Row[p];
		V_RHS[i] -= V_Matrix_Value[i] * lb;
	}
	V_Cost_Intercept += V_Cost[j] * lb;
	if (V_UB[j] <= MaxFinite)
		V_UB[j] -= V_LB[j];
	V_LB[j] = 0.0;

	Presolve_Modified ++;
}

int Presolve_Simple_Col_Check()
{
	for (int j = 0; j < n_Col; j ++)
	{
		if (Col_Disable[j])
			continue;
		// (iii) An infeasible variable
		// Exist j: l_j > u_j, trivially infeasible.
		if (V_LB[j] > V_UB[j])
		{
			LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
			return 0;
		}
		// (iv) A fixed variable
		// Exist j: l_j = u_j, x_j can be substituted out of the problem.
		if (V_UB[j] - V_LB[j] < Variable_Tolerance)
		{
			Presolve_Fix_Variable(j, V_LB[j]);
			continue;
		}
		// (ii) An empty column
		// Exist j: A_{ij} = 0, forall i
		// a. If c_j = 0, arbitrary x_j
		// b. If c_j < 0 and x_j have no upper bound: Unbounded
		// c. If c_j < 0 and x_j have upper bound: Set x_j = u_j
		// d. If c_j > 0 and x_j have no lower bound: Unbounded
		// e. If c_j > 0 and x_j have lower bound: Set x_j = l_j
		if (Col_Element_Count[j] == 0)
		{
			if (fabs(V_Cost[j]) < Variable_Tolerance)
				Presolve_Fix_Variable(j, 0); // After crushing, either 0 <= x_j <= u_j (can be +Inf) or x_j free
			else if (V_Cost[j] < 0)
			{
				if (V_UB[j] > MaxFinite)
				{
					LP_Status = LP_STATUS_PRIMAL_UNBOUNDED;
					return 0;
				}
				else
					Presolve_Fix_Variable(j, V_UB[j]);
			}
			else if (V_Cost[j] > 0)
			{
				if (V_LB[j] < -MaxFinite)
				{
					LP_Status = LP_STATUS_PRIMAL_UNBOUNDED;
					return 0;
				}
				else
					Presolve_Fix_Variable(j, V_LB[j]);
			}
		}
	}
	return 1;
}

int Presolve_Null_Row()
{
	// (i) An empty row
	// Exist i: a_{ij} = 0, forall j
	// Either redundant or infeasible.
	for (int i = 0; i < n_Row; i ++)
	{
		if (Row_Disable[i])
			continue;
		
		if (Row_Element_Count[i] == 0 || Row_1Norm[i] < Input_Tolerance)
		{
			if (fabs(V_RHS[i]) > Input_Tolerance)
			{
				LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
				return 0;
			}
			Presolve_Delete_Row(i);
			Row_Disable[i] = 1;
		}
	}

	return 1;
}

int Presolve_Singleton_Row()
{
	// (v) A singleton row
	// Exist (i, k): a_{ij} = 0, forall j != k, a_{ik} != 0
	// x_k = b_i / a_{ik}
	Presolve_Linked_List_Head = Presolve_Linked_List_Tail = -1;
	for (int i = 0; i < n_Row; i ++)
	{
		if (Row_Disable[i])
			continue;
		if (Row_Element_Count[i] == 1)
		{
			if (Presolve_Linked_List_Head == -1)
				Presolve_Linked_List_Head = i;
			else
				Presolve_Linked_List_Next[Presolve_Linked_List_Tail] = i;
			Presolve_Linked_List_Tail = i;
		}
	}
	if (Presolve_Linked_List_Tail == -1)
		return 1;
	for (int i = Presolve_Linked_List_Head; ; i = Presolve_Linked_List_Next[i])
	{
		int p = V_Matrix_Row_Head[i];
		int j = V_Matrix_Col[p];
		double a_ij = V_Matrix_Value[p];
		double b_i = V_RHS[i];
		double x_j = b_i / a_ij;
		if (x_j < V_LB[j] || x_j > V_UB[j])
		{
			LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
			return 0;
		}
		for (p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
		{
			int i = V_Matrix_Row[p];
			if (Row_Element_Count[i] == 2) // Eliminate a doubleton row
			{
				Presolve_Linked_List_Next[Presolve_Linked_List_Tail] = i;
				Presolve_Linked_List_Tail = i;
			}
		}
		Presolve_Fix_Variable(j, x_j);
		Presolve_Delete_Row(i);
		Row_Disable[i] = 1;
		if (i == Presolve_Linked_List_Tail)
			break;
	}
	return 1;
}

double Row_UB[MAX_ROWS], Row_LB[MAX_ROWS];

int Presolve_Forcing_Row()
{
	// Sect 3.3
	// g_i = sum_{j \in P_i} a_{ij} l_j + sum_{j \in M_i} a_{ij} u_j
	// h_i = sum_{j \in M_i} a_{ij} l_j + sum_{j \in P_i} a_{ij} u_j
	// where P_i = {j: a_{ij} > 0}, M_i = {j: a_{ij} < 0}
	for (int i = 0; i < n_Row; i ++)
	{
		if (Row_Disable[i])
			continue;
		Row_UB[i] = Row_LB[i] = 0;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int j = V_Matrix_Col[p];
			double a_ij = V_Matrix_Value[p];

			if (Row_UB[i] <= MaxFinite)
			{
				if (a_ij > 0)
				{
					if (V_UB[j] > MaxFinite)
						Row_UB[i] = MaxPositive;
					else
						Row_UB[i] += a_ij * V_UB[j];
				}
				else
				{
					if (V_LB[j] < -MaxFinite)
						Row_UB[i] = MaxPositive;
					else
						Row_UB[i] += a_ij * V_LB[j];
				}
			}

			if (Row_LB[i] >= -MaxFinite)
			{
				if (a_ij > 0)
				{
					if (V_LB[j] < -MaxFinite)
						Row_LB[i] = -MaxPositive;
					else
						Row_LB[i] += a_ij * V_LB[j];
				}
				else
				{
					if (V_UB[j] > MaxFinite)
						Row_LB[i] = -MaxPositive;
					else
						Row_LB[i] += a_ij * V_UB[j];
				}
			}
		}
	}

	for (int i = 0; i < n_Row; i ++)
	{
		if (Row_Disable[i])
			continue;
		// (ix) An infeasible constraint: 
		// Exist i: h_i < b_i or b_i < g_i
		if (V_RHS[i] < Row_LB[i] || V_RHS[i] > Row_UB[i])
		{
			LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
			return 0;
		}
		if (V_RHS[i] - Row_LB[i] >= Variable_Tolerance && Row_UB[i] - V_RHS[i] >= Variable_Tolerance)
			continue;
		// (x) A forcing constraint:
		// Exist i: g_i = b_i or h_i = b_i
		// Fix all variables in the i-th constraint.
		for (int p = V_Matrix_Row_Head[i]; p != -1; )
		{
			int j = V_Matrix_Col[p];
			double a_ij = V_Matrix_Value[p];
			int next_p = V_Matrix_Row_Next[p];
			if ((a_ij > 0 && Row_UB[i] - V_RHS[i] < Variable_Tolerance) || 
				(a_ij < 0 && V_RHS[i] - Row_LB[i] < Variable_Tolerance))
				Presolve_Fix_Variable(j, V_UB[i]);
			else
				Presolve_Fix_Variable(j, V_LB[i]);
			p = next_p;
		}
	}
	return 1;
}

int Presolve_Dominated_Row()
{
	for (int j = 0; j < n_Col; j ++)
	{
		if (Col_Disable[j])
			continue;
		if (Col_Element_Count[j] != 1)
			continue;
		int p = V_Matrix_Col_Head[j];
		int i = V_Matrix_Row[p];
		double a_ij = V_Matrix_Value[p];
		if (fabs(a_ij) < Input_Tolerance)
		{
			LP_Status = LP_STATUS_OVERFLOW;
			return 0;
		}
		double min_j = V_RHS[i]; // min_j: lower bound of x_j by the bounds of other variables in row i
		double max_j = V_RHS[i]; // max_j: upper bound of x_j by the bounds of other variables in row i
		for (p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int k = V_Matrix_Col[p];
			if (j == k)
				continue;
			double a_ik = V_Matrix_Value[p];
			if (min_j >= -MaxFinite)
			{
				if (a_ik > 0)
				{
					if (V_UB[k] > MaxFinite)
						min_j = -MaxPositive;
					else
						min_j -= a_ik * V_UB[k];
				}
				else
				{
					if (V_LB[k] < -MaxFinite)
						min_j = -MaxPositive;
					else
						min_j -= a_ik * V_LB[k];
				}
			}
			if (max_j <= MaxFinite)
			{
				if (a_ik > 0)
				{
					if (V_LB[k] < -MaxFinite)
						max_j = MaxPositive;
					else
						max_j -= a_ik * V_LB[k];
				}
				else
				{
					if (V_UB[k] > MaxFinite)
						max_j = MaxPositive;
					else
						max_j -= a_ik * V_UB[k];
				}
			}
		}

		if (min_j >= -MaxFinite)
			min_j /= a_ij;
		else if (a_ij < 0)
			min_j = -min_j;

		if (max_j <= MaxFinite)
			max_j /= a_ij;
		else if (a_ij < 0)
			max_j = -max_j;

		if (min_j > max_j)
			swap(min_j, max_j);

		// (viii) An implied free column singleton: 
		// Exist j, k: (a_{ij} = 0, forall i != k; a_{kj} != 0), l_j <= min_j <= max_j <= u_j.
		// 
		// [min_j, max_j] is the implied bound for x_j
		// Case 1: If l_j <= min_j && max_j <= u_j, 
		//         the original bound x_j \in [l_j, u_j] is equivalent to x_j free, 
		//         the free variable x_j can be replaced by (b_i - sum{k != j} a_{ik} x_k) / a_{ij}
		// Case 2: Otherwise, use [min_j, max_j] to tighten the original bound. 
		if (V_LB[j] <= min_j && max_j <= V_UB[j])
		{
			// c_j x_j = (c_j / a_{ij}) * (b_i - sum{k != j} a_{ik} x_k)
			double c_j_DIV_a_ij = V_Cost[j] / a_ij;
			for (p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
			{
				int k = V_Matrix_Col[p];
				if (k != j)
					V_Cost[k] -= c_j_DIV_a_ij * V_Matrix_Value[p];
			}
			V_Cost_Intercept += c_j_DIV_a_ij * V_RHS[i];

			V_Presolve_Linear_Replace[j] = V_Matrix_Row_Head[i]; // Finally recover x_j by a linear relation
			Presolve_Delete_Row(i);
			Row_Disable[i] = 1;
		}
		else
		{
			V_UB[j] = min(max_j, V_UB[j]);
			if (V_LB[j] < min_j)
				Presolve_Set_Variable_LB(j, min_j);
			if (V_LB[j] > V_UB[j])
			{
				LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
				return 0;
			}
		}
	}
	return 1;
}

int Presolve_Doubleton_Row_Singleton_Col()
{
	for (int j = 0; j < n_Col; j ++)
	{
		if (Col_Disable[j])
			continue;
		if (Col_Element_Count[j] != 1)
			continue;
		int p = V_Matrix_Col_Head[j];
		int i = V_Matrix_Row[p];
		if (Row_Element_Count[i] != 2)
			continue;
		double a_ij = V_Matrix_Value[p];
		int k;
		for (p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			k = V_Matrix_Col[p];
			if (k != j)
				break;
		}
		double a_ik = V_Matrix_Value[p];
		if (fabs(a_ij) < Input_Tolerance || fabs(a_ik) < Input_Tolerance)
		{
			LP_Status = LP_STATUS_OVERFLOW;
			return 0;
		}
		// (vii) A doubleton equation combined with a column singleton:
		// Exist i, j, k: a_ij * x_j + a_ik * x_k = b_i, j != k, a_ij != 0, a_ik != 0
		// x_j is a singleton column.
		
		// Another implied bound for x_k: [min_k, max_k]
		double min_k, max_k;
		if (a_ij * a_ik > 0)
		{
			max_k = (V_LB[j] < -MaxFinite) ? MaxPositive : ((V_RHS[i] - a_ij * V_LB[j]) / a_ik);
			min_k = (V_UB[j] > MaxFinite) ? -MaxPositive : ((V_RHS[i] - a_ij * V_UB[j]) / a_ik);
		}
		else
		{
			min_k = (V_LB[j] < -MaxFinite) ? -MaxPositive : ((V_RHS[i] - a_ij * V_LB[j]) / a_ik);
			max_k = (V_UB[j] > MaxFinite) ? MaxPositive : ((V_RHS[i] - a_ij * V_UB[j]) / a_ik);
		}
		printf("j = %d, max_k = %e, min_k = %e\n", j, max_k, min_k);

		// Update the original bound
		V_UB[k] = min(max_k, V_UB[k]);
		if (V_LB[k] < min_k)
			Presolve_Set_Variable_LB(k, min_k);
		if (V_LB[k] > V_UB[k])
		{
			LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
			return 0;
		}
		// x_j = (b_i / a_ij) - (a_ik / a_ij) x_k, which can be substituted out of the problem.
		// We can show that the original bound of x_j is implied free after the bound of x_k is updated!
		double c_j_DIV_a_ij = V_Cost[j] / a_ij;
		V_Cost[k] -= c_j_DIV_a_ij * a_ik;
		V_Cost_Intercept += c_j_DIV_a_ij * V_RHS[i];

		V_Presolve_Linear_Replace[j] = V_Matrix_Row_Head[i]; // Finally recover x_j by a linear relation
		Presolve_Delete_Row(i);
		Row_Disable[i] = 1;
	}
	return 1;
}

double Presolve_Dual_UB[MAX_ROWS], Presolve_Dual_LB[MAX_COLS];

int Presolve_Dominated_Col()
{
	// The optimality condition contains: A^T y^* + z^* = c
	
	// Set original bounds for dual variables: y free
	for (int i = 0; i < n_Row; i ++)
	{
		Presolve_Dual_LB[i] = -MaxPositive;
		Presolve_Dual_UB[i] = MaxPositive;
	}
	// Use all column singletons to calc implied bounds for dual variables
	for (int j = 0; j < n_Col; j ++)
	{
		if (Col_Disable[j])
			continue;
		if (Col_Element_Count[j] != 1)
			continue;
		int p = V_Matrix_Col_Head[j];
		int i = V_Matrix_Row[p];
		double a_ij = V_Matrix_Value[p];

		double c_j_DIV_a_ij = V_Cost[j] / a_ij;
		// Table 2 in Andersen and Andersen (1995)
		if (V_UB[j] > MaxFinite)
		{
			// 1. x_j free, has been treated! 
			// 2. x_j has LB and no UB, a_ij > 0 => y_i^* <= c_j / a_ij
			// 3. x_j has LB and no UB, a_ij < 0 => y_i^* >= c_j / a_ij
			if (a_ij > 0)
				Presolve_Dual_UB[i] = min(Presolve_Dual_UB[i], c_j_DIV_a_ij);
			else
				Presolve_Dual_LB[i] = max(Presolve_Dual_LB[i], c_j_DIV_a_ij);
		}
		// 4. 5. x_j has UB and no LB, has been transformed to 2. 3. ! 
	}
	for (int j = 0; j < n_Col; j ++)
	{
		if (Col_Disable[j])
			continue;
		// e_j = sum_{i \in P_j} a_ij dl_i + sum_{i \in M_j} a_ij du_j
		// d_j = sum_{i \in P_j} a_ij du_i + sum_{i \in M_j} a_ij dl_j
		// where P_j = {i: a_ij > 0}, M_j = {i: a_ij < 0}
		// The property is e_j <= sum_i a_ij y_i^* <= d_j
		double e_j = 0, d_j = 0;
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
		{
			int i = V_Matrix_Row[p];
			double a_ij = V_Matrix_Value[p];
			if (a_ij > 0)
			{
				if (Presolve_Dual_LB[i] < -MaxFinite)
					e_j = -MaxPositive;
				else
					e_j += a_ij * Presolve_Dual_LB[i];

				if (Presolve_Dual_UB[i] > MaxFinite)
					d_j = MaxPositive;
				else
					d_j += a_ij * Presolve_Dual_UB[i];
			}
			else
			{
				if (Presolve_Dual_LB[i] < -MaxFinite)
					d_j = MaxPositive;
				else
					d_j += a_ij * Presolve_Dual_LB[i];

				if (Presolve_Dual_UB[i] > MaxFinite)
					e_j = -MaxPositive;
				else
					e_j += a_ij * Presolve_Dual_UB[i];
			}
		}
		double c_j_SUB_d_j = (d_j > MaxFinite) ? -MaxPositive : (V_Cost[j] - d_j);
		double c_j_SUB_e_j = (e_j < -MaxFinite) ? MaxPositive : (V_Cost[j] - e_j);
		
		// Table 1 in Andersen and Andersen (1995)
		// 1. l_j > -Inf, u_j = +Inf, z_j^* > 0, x^*_j = l_j
		// 2. l_j = -Inf, u_j < +Inf, z_j^* < 0, x^*_j = u_j
		// 3. l_j = -Inf, z_j^* > 0, Unbounded
		// 4. u_j = +Inf, z_j^* < 0, Unbounded
		// 5. u_j < l_j, Infeasible
		
		// By optimality, A^T y^* + z^* = c
		// c_j - d_j <= z_j^* <= c_j - e_j
		
		// (xii) A weakly dominated column
		// S = {j: a_ij is a column singleton for some i}
		// Exist j not in S: c_j - d_j = 0 and l_j > -Inf
		//                   => x^*_j = l_j
		if (fabs(c_j_SUB_d_j) < Variable_Tolerance)
		{
			// Not a column singleton
			if (Col_Element_Count[j] != 1 && V_LB[j] >= -MaxFinite)
				Presolve_Fix_Variable(j, V_LB[j]);
		}
		// Exist j not in S: c_j - e_j = 0 and u_j < +Inf 
		//                   => x^*_j = u_j
		else if (fabs(c_j_SUB_e_j) < Variable_Tolerance)
		{
			// Not a column singleton
			if (Col_Element_Count[j] != 1 && V_UB[j] <= MaxFinite)
				Presolve_Fix_Variable(j, V_UB[j]);
		}
		// (xi) A dominated column
		// Exist j: c_j - d_j > 0 or c_j - e_j < 0
		// c_j - d_j > 0 => z_j^* > 0 => { x^*_j = l_j, l_j > -Inf
		//                               { Unbounded,   l_j = -Inf
		else if (c_j_SUB_d_j > 0)
		{
			if (V_LB[j] < -MaxFinite)
			{
				LP_Status = LP_STATUS_PRIMAL_UNBOUNDED;
				return 0;
			}
			Presolve_Fix_Variable(j, V_LB[j]);
		}
		// c_j - e_j < 0 => z_j^* < 0 => { x^*_j = u_j, u_j < +Inf
		//                               { Unbounded,   u_j = +Inf
		else if (c_j_SUB_e_j < 0)
		{
			if (V_UB[j] > MaxFinite)
			{
				LP_Status = LP_STATUS_PRIMAL_UNBOUNDED;
				return 0;
			}
			Presolve_Fix_Variable(j, V_UB[j]);
		}
		// e_j and d_j are used to compute new bounds of y_j^*
		else
		{
			if (V_UB[j] > MaxFinite)
			{
				for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
				{
					int i = V_Matrix_Row[p];
					double a_ij = V_Matrix_Value[p];
					if (fabs(a_ij) < Input_Tolerance)
					{
						LP_Status = LP_STATUS_OVERFLOW;
						return 0;
					}
					if (a_ij > 0)
					{
						if (V_LB[j] >= -MaxFinite)
						{
							double max_yi = (c_j_SUB_e_j > MaxFinite) ? MaxPositive : (c_j_SUB_e_j / a_ij + Presolve_Dual_LB[i]);
							Presolve_Dual_UB[i] = min(Presolve_Dual_UB[i], max_yi);
						}
						else
						{
							double min_yi = (c_j_SUB_d_j < -MaxFinite) ? -MaxPositive : (c_j_SUB_d_j / a_ij + Presolve_Dual_UB[i]);
							Presolve_Dual_LB[i] = max(Presolve_Dual_LB[i], min_yi);
						}
					}
					else
					{
						if (V_LB[j] >= -MaxFinite)
						{
							double min_yi = (c_j_SUB_e_j > MaxFinite) ? -MaxPositive : (c_j_SUB_e_j / a_ij + Presolve_Dual_UB[i]);
							Presolve_Dual_LB[i] = max(Presolve_Dual_LB[i], min_yi);
						}
						else
						{
							double max_yi = (c_j_SUB_d_j < -MaxFinite) ? MaxPositive : (c_j_SUB_d_j / a_ij + Presolve_Dual_UB[i]);
							Presolve_Dual_UB[i] = min(Presolve_Dual_UB[i], max_yi);
						}
					}
				}
			}
		}
	}
	return 1;
}

int Presolve_Row_Cnt_NonSingleton[MAX_ROWS], Presolve_Row_Cnt_Singleton[MAX_ROWS];
int Presolve_Row_Singleton_Pos[MAX_ROWS][2]; // Store the pointer, not col ID! 
double Presolve_Row_First_NonSingleton_Val[MAX_ROWS]; // Store the value! 
int Presolve_Row_Same_NonSingleton_Head[MAX_COLS], Presolve_Row_Same_NonSingleton_Next[MAX_ROWS];

int Presolve_Duplicate_Row()
{
	for (int j = 0; j < n_Col; j ++)
		Presolve_Row_Same_NonSingleton_Head[j] = -1;
	for (int i = 0; i < n_Row; i ++)
	{
		if (Row_Disable[i])
			continue;
		Presolve_Row_Cnt_NonSingleton[i] = Presolve_Row_Cnt_Singleton[i] = 0;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			int j = V_Matrix_Col[p];
			if (Col_Element_Count[j] == 1) // Column singleton
			{
				if (Presolve_Row_Cnt_Singleton[i] <= 1) // Store the first two
					Presolve_Row_Singleton_Pos[i][Presolve_Row_Cnt_Singleton[i]] = p;
				Presolve_Row_Cnt_Singleton[i] ++;
			}
			else
			{
				if (Presolve_Row_Cnt_NonSingleton[i] == 0)
					Presolve_Row_First_NonSingleton_Val[i] = V_Matrix_Value[p];
				Presolve_Row_Cnt_NonSingleton[i] ++;
			}
		}
		if (Presolve_Row_Cnt_Singleton[i] <= 2) // Only process rows with no more than 2 column singletons
		{
			Presolve_Row_Same_NonSingleton_Next[i] = Presolve_Row_Same_NonSingleton_Head[Presolve_Row_Cnt_NonSingleton[i]];
			Presolve_Row_Same_NonSingleton_Head[Presolve_Row_Cnt_NonSingleton[i]] = i;
		}
	}
	for (int ns = 1; ns < n_Col; ns ++) // Ignore no nonsingletons
	{
		for (int i = Presolve_Row_Same_NonSingleton_Head[ns]; i != -1; i = Presolve_Row_Same_NonSingleton_Next[i])
			for (int k = Presolve_Row_Same_NonSingleton_Head[ns]; k != -1; k = Presolve_Row_Same_NonSingleton_Next[k])
			{
				if (i == k || Row_Disable[i] || Row_Disable[k])
					continue;
				if (Presolve_Row_Cnt_Singleton[i] > Presolve_Row_Cnt_Singleton[k]) // Symmetric
					continue;
				if (Presolve_Row_Cnt_Singleton[i] + Presolve_Row_Cnt_Singleton[k] > 2) // Only process rows with no more than 2 column singletons (added up)
					continue;
				// Check whether Row i and Row k is proportional (except column singletons)
				int pi = V_Matrix_Row_Head[i];
				int pk = V_Matrix_Row_Head[k];
				// Suppose Row_i * v = Row_k (except column singletons)
				// Intend to clear Row k
				double v = Presolve_Row_First_NonSingleton_Val[k] / Presolve_Row_First_NonSingleton_Val[i];
				while (true)
				{
					while (pi != -1 && Col_Element_Count[V_Matrix_Col[pi]] == 1)
						pi = V_Matrix_Row_Next[pi];
					while (pk != -1 && Col_Element_Count[V_Matrix_Col[pk]] == 1)
						pk = V_Matrix_Row_Next[pk];
					if (pi == -1 || pk == -1) // At least one row finish
						break;
					if (V_Matrix_Col[pi] != V_Matrix_Col[pk])
						break;
					if (fabs(V_Matrix_Value[pi] * v - V_Matrix_Col[pk]) > Input_Tolerance)
						break;
					pi = V_Matrix_Row_Next[pi];
					pk = V_Matrix_Row_Next[pk];
				}
				if (pi != -1 || pk != -1) // Not proportional
					continue;
				// Case 1: (0, 0), Either Row k can be completely cleared or infeasible
				if (Presolve_Row_Cnt_Singleton[i] == 0 && Presolve_Row_Cnt_Singleton[k] == 0)
				{
					if (fabs(V_RHS[i] * v - V_RHS[k]) > Input_Tolerance)
					{
						LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
						return 0;
					}
					Presolve_Delete_Row(k);
					Row_Disable[k] = 1;
				}
				// Case 2: (0, 1), Either fix x_j = (b_k - v * b_i) / a_kj or infeasible
				else if (Presolve_Row_Cnt_Singleton[i] == 0 && Presolve_Row_Cnt_Singleton[k] == 1)
				{
					int j = V_Matrix_Col[Presolve_Row_Singleton_Pos[k][0]];
					double a_kj = V_Matrix_Value[Presolve_Row_Singleton_Pos[k][0]];
					double x_j = (V_RHS[k] - v * V_RHS[i]) / a_kj;
					if (x_j < V_LB[j] || x_j > V_UB[j])
					{
						LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
						return 0;
					}
					Presolve_Fix_Variable(j, x_j);
					Presolve_Delete_Row(k);
					Row_Disable[k] = 1;
				}
				// Case 3: (0, 2) and Case 4: (1, 1)
				else if (Presolve_Row_Cnt_Singleton[i] + Presolve_Row_Cnt_Singleton[k] == 2)
				{
					int j, w;
					double a_j, a_w;
					// a_j x_j - a_w x_w = b_k - v * b_i
					if (Presolve_Row_Cnt_Singleton[i] == 0 && Presolve_Row_Cnt_Singleton[i] == 2)
					{
						// Column singletons: a_kj and a_kw
						// After elimination, a_kj x_j - (-a_kw) x_w = b_k - v * b_i
						j = V_Matrix_Col[Presolve_Row_Singleton_Pos[k][0]];
						w = V_Matrix_Col[Presolve_Row_Singleton_Pos[k][1]];
						a_j = V_Matrix_Value[Presolve_Row_Singleton_Pos[k][0]];
						a_w = -V_Matrix_Value[Presolve_Row_Singleton_Pos[k][0]];
					}
					else // (1, 1)
					{
						// Column singletons: a_ij and a_kw
						// After elimination, a_kw x_w - v * a_ij x_j = b_k - v * b_i
						j = V_Matrix_Col[Presolve_Row_Singleton_Pos[i][0]];
						w = V_Matrix_Col[Presolve_Row_Singleton_Pos[k][0]];
						a_j = v * V_Matrix_Value[Presolve_Row_Singleton_Pos[i][0]];
						a_w = V_Matrix_Value[Presolve_Row_Singleton_Pos[k][0]];
					}
					// x_j = (a_w x_w - (b_k - v * b_i)) / a_j;
					double b_k_MINUS_vb_i = V_RHS[k] - v * V_RHS[i];
					double l_j, u_j;
					if (a_j * a_w > 0)
					{
						u_j = (V_UB[w] > MaxFinite) ? MaxPositive : ((a_w * V_UB[w] - b_k_MINUS_vb_i) / a_j);
						l_j = (V_LB[w] < -MaxFinite) ? -MaxPositive : ((a_w * V_LB[w] - b_k_MINUS_vb_i) / a_j);
					}
					else
					{
						l_j = (V_UB[w] > MaxFinite) ? -MaxPositive : ((a_w * V_UB[w] - b_k_MINUS_vb_i) / a_j);
						u_j = (V_LB[w] < -MaxFinite) ? MaxPositive : ((a_w * V_LB[w] - b_k_MINUS_vb_i) / a_j);
					}
					l_j = max(l_j, V_LB[j]);
					u_j = min(u_j, V_UB[j]);
					if (l_j > u_j)
					{
						LP_Status = LP_STATUS_PRIMAL_INFEASIBLE;
						return 0;
					}
					// After bound of x_j is modified to [l_j, u_j], x_w can be moved out like (vii)
					V_UB[j] = u_j;
					if (V_LB[j] < l_j)
						Presolve_Set_Variable_LB(j, l_j);
					// x_w = (b_k - v * b_i) / a_w + (a_j / a_w) x_j
					double c_w_DIV_a_w = V_Cost[w] / a_w;
					V_Cost[j] += c_w_DIV_a_w * a_j;
					V_Cost_Intercept += c_w_DIV_a_w * b_k_MINUS_vb_i;
					V_Presolve_Linear_Replace[w] = V_Matrix_Row_Head[k]; // Finally recover x_w by a linear relation
					Presolve_Delete_Row(k);
					Row_Disable[k] = 1;
				}
			}
	}
	return 1;
}

int Presolve_Col_Same_Head[MAX_ROWS], Presolve_Col_Same_Next[MAX_COLS];

int Presolve_Duplicate_Col()
{
	// Duplicate Cols: Exist j, k: a_ij = v * a_ik, forall i, j != k
	for (int i = 0; i < n_Row; i ++)
		Presolve_Col_Same_Head[i] = -1;
	for (int j = 0; j < n_Col; j ++)
	{
		if (Col_Disable[j])
			continue;
		Presolve_Col_Same_Next[j] = Presolve_Col_Same_Head[Col_Element_Count[j]];
		Presolve_Col_Same_Head[Col_Element_Count[j]] = j;
	}
	for (int cnt = 1; cnt < n_Row; cnt ++)
	{
		for (int j = Presolve_Col_Same_Head[cnt]; j != -1; j = Presolve_Col_Same_Next[j])
			for (int k = Presolve_Col_Same_Head[cnt]; k != -1; k = Presolve_Col_Same_Next[k])
			{
				if (j == k || Col_Disable[j] || Col_Disable[k])
					continue;
				// Check whether Col j and Col k is proportional
				int pj = V_Matrix_Col_Head[j];
				int pk = V_Matrix_Col_Head[k];
				// Suppose Col_j = v * Col_k
				// Intend to clear Row k
				double v = V_Matrix_Value[pk] / V_Matrix_Value[pj];
				while (pj != -1 && pk != -1)
				{
					if (V_Matrix_Col[pj] != V_Matrix_Col[pk])
						break;
					if (fabs(V_Matrix_Value[pj] - V_Matrix_Col[pk] * v) > Input_Tolerance)
						break;
					pj = V_Matrix_Col_Next[pj];
					pk = V_Matrix_Col_Next[pk];
				}
				if (pj != -1 || pk != -1) // Not proportional
					continue;
				double c_j_SUB_vc_k = V_Cost[j] - V_Cost[k] * v;
				if (fabs(c_j_SUB_vc_k) < Input_Tolerance)
				{
					// Because it is somehow too sophisticated for this case in Postsolve, we don't implement this part. 
					// (xvii) Replacing two duplicate columns by one
					// c_j - v c_k = 0
					// Original ...  c_j x_j +  c_k x_k ...
					//          ... a_.j x_j + a_.k x_k ... = b
					//          l_j <= x_j <= u_j, l_k <= x_k <= u_k
					// Define x_k = -v x_j + x'_k, 
					// Transformed ...  c_k x'_k ...
					//             ... a_.k x'_k ... = b
					//             l_j <= x_j <= u_j, l_k <= x'_k - v x_j <= u_k
					// The bound of x'_k is determined by the sign of v
					// x_k is replaced by x'_k and x_j is removed
					/*
					double d_k, e_k;
					// l_k + v x_j <= x'_k <= u_k + v x_j
					// d_k <= x'_k <= e_k
					if (v > 0)
					{
						d_k = (V_LB[k] < -MaxFinite || V_LB[j] < -MaxFinite) ? -MaxPositive : (V_LB[k] + v * V_LB[j]);
						e_k = (V_UB[k] > MaxFinite || V_UB[j] > MaxFinite) ? MaxPositive : (V_UB[k] + v * V_UB[j]);
					}
					else
					{
						d_k = (V_LB[k] < -MaxFinite || V_UB[j] > MaxFinite) ? -MaxPositive : (V_LB[k] + v * V_UB[j]);
						e_k = (V_UB[k] > MaxFinite || V_LB[j] < -MaxFinite) ? MaxPositive : (V_UB[k] + v * V_LB[j]);
					}
					*/
				}
				// (xvi) Fixing a duplicate column
				// c_j - v c_k != 0
				// Table 3:
				// u_k = +Inf, z_k^* >= 0, v >= 0, c_j - v c_k > 0, z_j^* > 0
				// l_k = -Inf, z_k^* <= 0, v <= 0, c_j - v c_k > 0, z_j^* > 0
				// u_k = +Inf, z_k^* >= 0, v <= 0, c_j - v c_k < 0, z_j^* < 0
				// l_k = -Inf, z_k^* <= 0, v >= 0, c_j - v c_k < 0, z_j^* < 0
				else
				{
					// u_k = +Inf, z_k^* >= 0, v >= 0, c_j - v c_k > 0 => z_j^* > 0
					// l_k = -Inf, z_k^* <= 0, v <= 0, c_j - v c_k > 0 => z_j^* > 0
					if (c_j_SUB_vc_k > 0 && ((v > 0 && V_UB[k] > MaxFinite) || (v < 0 && V_LB[k] < -MaxFinite)))
					{
						// z_j^* > 0 => { x_j^* = l_j, l_j > -Inf
						//              { Unbounded,   l_j = -Inf
						if (V_LB[j] < -MaxFinite)
						{
							LP_Status = LP_STATUS_PRIMAL_UNBOUNDED;
							return 0;
						}
						Presolve_Fix_Variable(j, V_LB[j]);
					}
					// u_k = +Inf, z_k^* >= 0, v <= 0, c_j - v c_k < 0 => z_j^* < 0
					// l_k = -Inf, z_k^* <= 0, v >= 0, c_j - v c_k < 0 => z_j^* < 0
					else if (c_j_SUB_vc_k < 0 && ((v < 0 && V_UB[k] > MaxFinite) || (v > 0 && V_LB[k] < -MaxFinite)))
					{
						// z_j^* < 0 => { x_j^* = u_j, u_j < +Inf
						//              { Unbounded,   u_j = +Inf
						if (V_UB[j] > MaxFinite)
						{
							LP_Status = LP_STATUS_PRIMAL_UNBOUNDED;
							return 0;
						}
						Presolve_Fix_Variable(j, V_UB[j]);
					}
				}
			}
	}
	return 1;
}

int Temp[MAX_COLS + MAX_ROWS];

void Presolve_FinalizeModel()
{
	// Preserve the Original Model and Finalize the Reduced Model
	// First backup some necessary information
	n_Row_ORIG = n_Row;
	n_Col_ORIG = n_Col;

	// I. Columns Reorder
	// Counting n_LB, n_UB, n_FR
	int n_RealCol = 0;
	n_LB = n_UB = n_FR = 0;
	for (int j = 0; j < n_Col; j ++)
	{
		if (Col_Disable[j])
			continue;
		n_RealCol ++;
		if (V_LB[j] == 0.0) // With LB
		{
			n_LB ++;
			if (V_UB[j] <= MaxFinite) // With LB and UB
				n_UB ++;
		}
		else
			n_FR ++;
	}
	// Rearranged as: (a) x[0 ~ (n_UB - 1)]: With LB and UB; 
	//                (b) x[n_UB ~ (n_LB - 1)]: With LB only; 
	//                (c) x[n_LB ~ (n_LB + n_FR - 1)]: Free; 
	//                (d) x[(n_LB + n_FR) ~ n_Col]: Deleted.
	int Ptr1 = 0, Ptr2 = n_UB, Ptr3 = n_LB, Ptr4 = n_LB + n_FR;
	for (int j = 0; j < n_Col; j ++)
	{
		if (Col_Disable[j])
			Col_NewToOld[Ptr4 ++] = j;
		else if (V_LB[j] == 0.0 && V_UB[j] <= MaxFinite) // With LB and UB
			Col_NewToOld[Ptr1 ++] = j;
		else if (V_LB[j] == 0.0) // With LB only
			Col_NewToOld[Ptr2 ++] = j;
		else
			Col_NewToOld[Ptr3 ++] = j;
	}
	// To call original index, just use x[Col_OldToNew[i]] instead of x[i]. 
	for (int j = 0; j < n_Col; j ++)
		Col_OldToNew[Col_NewToOld[j]] = j;
	// Perform Reordering
	memset(Temp, 0, sizeof(int) * n_Col);
	for (int i = 0; i < n_Col; i ++)
	{
		int j = i;
		int Next = Col_NewToOld[j];
		if (j == Next || Temp[i])
			continue;
		Temp[j] = 1;
		while (true)
		{
			swap(V_Cost[j], V_Cost[Next]);
			swap(V_LB[j], V_LB[Next]);
			swap(V_UB[j], V_UB[Next]);
			swap(V_Matrix_Col_Head[j], V_Matrix_Col_Head[Next]);
			swap(V_Crushing_Times[j], V_Crushing_Times[Next]);
			swap(V_Crushing_Add[j], V_Crushing_Add[Next]);
			swap(Col_Element_Count[j], Col_Element_Count[Next]);
			swap(Col_Disable[j], Col_Disable[Next]);
			Temp[Next] = 1;
			j = Next;
			Next = Col_NewToOld[j];
			if (Temp[Next])
				break;
		}
	}

	// II. Rows Reorder
	int n_RealRow = 0;
	for (int i = 0; i < n_Row; i ++)
	{
		if (Row_Disable[i])
			continue;
		n_RealRow ++;
	}
	Ptr1 = 0;
	Ptr2 = n_RealRow;
	for (int i = 0; i < n_Row; i ++)
		if (Row_Disable[i])
			Row_NewToOld[Ptr2 ++] = i;
		else
			Row_NewToOld[Ptr1 ++] = i;
	
	for (int i = 0; i < n_Row; i ++)
		Row_OldToNew[Row_NewToOld[i]] = i;
	memset(Temp, 0, sizeof(int) * n_Row);
	for (int i = 0; i < n_Row; i ++)
	{
		int j = i;
		int Next = Row_NewToOld[j];
		if (j == Next || Temp[i])
			continue;
		Temp[j] = 1;
		while (true)
		{
			swap(V_RHS[j], V_RHS[Next]);
			swap(V_Matrix_Row_Head[j], V_Matrix_Row_Head[Next]);
			swap(Row_Element_Count[j], Row_Element_Count[Next]);
			swap(Row_Disable[j], Row_Disable[Next]);
			Temp[Next] = 1;
			j = Next;
			Next = Row_NewToOld[j];
			if (Temp[Next])
				break;
		}
	}

	// III. Update V_Matrix_Row and V_Matrix_Col
	for (int i = 0; i < n_Element; i ++)
	{
		V_Matrix_Row[i] = Row_OldToNew[V_Matrix_Row[i]];
		V_Matrix_Col[i] = Col_OldToNew[V_Matrix_Col[i]];
	}

	// Note that this reorder does not affect every sorted columns

	// IV. Update n_Row and n_Col
	n_Row = n_RealRow;
	n_Col = n_RealCol;
	int nnz = 0;
	for (int j = 0; j < n_Col; j ++)
		nnz += Col_Element_Count[j];
	printf("After Presolving, n_Row = %d, n_Col = %d, n_Element = %d\n", n_Row, n_Col, nnz);
}

void Presolve_Recount_Delete_Disabled()
{
	for (int i = 0; i < n_Row; i ++)
		Row_Element_Count[i] = 0;
	for (int j = 0; j < n_Col; j ++)
		Col_Element_Count[j] = 0;
	for (int j = 0; j < n_Col; j ++)
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
		{
			int i = V_Matrix_Row[p];
			Row_Element_Count[i] ++;
			Col_Element_Count[j] ++;
		}
	for (int i = 0; i < n_Row; i ++)
		if (Row_Disable[i])
			Presolve_Delete_Row(i);
}

int Presolve_Main()
{
	Presolve_Init();

	if (PRESOLVE_LINDEP) // Check Linear Dependent Rows
	{
		Presolve_Linear_Dependent_Main();
		if (LP_Status != LP_STATUS_OK)
			return 0;
		Presolve_Recount_Delete_Disabled();
	}

	int Loop_Count = 0;
	do
	{
		Loop_Count ++;
		Presolve_Modified = 0;
		if (PRESOLVE_LEVEL >= 1)
		{
			Presolve_Simple_Col_Check();
			if (LP_Status != LP_STATUS_OK)
				return 0;
			Presolve_Null_Row();
			if (LP_Status != LP_STATUS_OK)
				return 0;
			Presolve_Singleton_Row();
			if (LP_Status != LP_STATUS_OK)
				return 0;
			Presolve_Forcing_Row();
			if (LP_Status != LP_STATUS_OK)
				return 0;
			Presolve_Dominated_Row();
			if (LP_Status != LP_STATUS_OK)
				return 0;
			Presolve_Doubleton_Row_Singleton_Col();
			if (LP_Status != LP_STATUS_OK)
				return 0;
		}
		if (PRESOLVE_LEVEL >= 2)
		{
			Presolve_Dominated_Col();
			if (LP_Status != LP_STATUS_OK)
				return 0;
		}
		if (PRESOLVE_LEVEL >= 3)
		{
			Presolve_Duplicate_Row();
			if (LP_Status != LP_STATUS_OK)
				return 0;
		}
		if (PRESOLVE_LEVEL >= 4)
		{
			Presolve_Duplicate_Col();
			if (LP_Status != LP_STATUS_OK)
				return 0;
		}
	}
	while (Presolve_Modified && Loop_Count < PRESOLVE_LOOP);

	Presolve_Null_Row();
	if (LP_Status != LP_STATUS_OK)
		return 0;
	
	Presolve_FinalizeModel();
	return 0;
}
