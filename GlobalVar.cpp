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

char Filename[MAX_FILENAME];
char Enabled_RHS[] = "";
char Enabled_RANGES[] = "";
char Enabled_BOUNDS[] = "";

// Problem Related
char Problem_Name[MAX_PROBLEM_NAME];
int n_Row, n_Col, n_Element;
char Row_Type[MAX_ROWS];
double V_Cost[MAX_COLS]; // c, Cost Row
double V_RHS[MAX_ROWS]; // b, RHS
double V_RHS_r[MAX_ROWS]; // RANGES, if Row_Type[i] == 'R', then [V_RHS_r[i], V_RHS[i]]
double V_LB[MAX_COLS], V_UB[MAX_COLS]; // For Variable x[j], V_LB[j] <= x[j] <= V_UB[j];
// A, Column Majored Matrix, Linked List
int V_Matrix_Col_Head[MAX_COLS], V_Matrix_Col_Next[MAX_ELEMENTS], V_Matrix_Col_Prev[MAX_ELEMENTS];
int V_Matrix_Row_Head[MAX_COLS], V_Matrix_Row_Next[MAX_ELEMENTS], V_Matrix_Row_Prev[MAX_ELEMENTS];
int V_Matrix_Row[MAX_ELEMENTS], V_Matrix_Col[MAX_ELEMENTS];
double V_Matrix_Value[MAX_ELEMENTS];

// Crushing
double V_Cost_Intercept; // After crushing, objective may have nonzero intercept
int V_Crushing_Times[MAX_COLS];
double V_Crushing_Add[MAX_COLS]; // Output (x[i] * V_Crushing_Times[i] + V_Crushing_Add[i])

// Presolve
int V_Presolve_Linear_Replace[MAX_COLS];
// { x_k = x'_k - v x_j
// { l_j <= x_j <= u_j
// { l_k <= x'_k - v x_j <= u_k
int V_Presolve_Duplicate_Column_Cnt;
double V_Presolve_Duplicate_Column_v[MAX_COLS];
int V_Presolve_Duplicate_Column_k[MAX_COLS], V_Presolve_Duplicate_Column_j[MAX_COLS];
// Final Model
int n_LB, n_UB, n_FR; 
// Rearranged as: (a) x[0 ~ (n_UB - 1)]: With LB and UB; 
//                (b) x[n_UB ~ (n_LB - 1)]: With LB only; 
//                (c) x[n_LB ~ (n_LB + n_FR - 1)]: Free.
int RecoverOrder[MAX_COLS], TransOrder[MAX_COLS];

// Problem Status
int LP_Status; // 0 - ok, 1 - infeasible, 2 - unbounded 
