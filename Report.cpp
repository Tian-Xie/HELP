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

#include "LP.h"

void Report_Solution()
{
	if (LP_Status == LP_STATUS_OVERFLOW)
	{
		printf("Problem Status: Overflow.\n");
		return;
	}
	if (LP_Status == LP_STATUS_PRIMAL_INFEASIBLE)
	{
		printf("Problem Status: Primal Infeasible.\n");
		return;
	}
	if (LP_Status == LP_STATUS_PRIMAL_UNBOUNDED)
	{
		printf("Problem Status: Primal Unbounded.\n");
		return;
	}
	if (HSD_Status == HSD_STATUS_STALLED)
	{
		printf("Problem Status: Numerical Difficulty.\n");
		return;
	}
	if (HSD_Status == HSD_STATUS_ITER_LIMIT)
	{
		printf("Problem Status: Iteration Limit Exceeded.\n");
		return;
	}
	if (HSD_Status == HSD_STATUS_UNKNOWN_ERROR)
	{
		printf("Problem Status: Unknown Error in HSD Algorithm.\n");
		return;
	}
	
	// TODO: Infeasibility from HSD Results

	// A Post-solve Procedure, 
	printf("Problem Status: Problem Solved.\n");



	/*
extern double V_Cost_Intercept; // After crushing, objective may have nonzero intercept
extern int V_Crushing_Times[MAX_COLS];
extern double V_Crushing_Add[MAX_COLS]; // Output (x[i] * V_Crushing_Times[i] + V_Crushing_Add[i])

// Presolve
extern int V_Presolve_Linear_Replace[MAX_COLS];
// { x_k = x'_k - v x_j
// { l_j <= x_j <= u_j
// { l_k <= x'_k - v x_j <= u_k
// Recover from newest to oldest!
extern int V_Presolve_Duplicate_Column_Cnt;
extern double V_Presolve_Duplicate_Column_v[MAX_COLS];
extern int V_Presolve_Duplicate_Column_k[MAX_COLS], V_Presolve_Duplicate_Column_j[MAX_COLS];
	*/

}