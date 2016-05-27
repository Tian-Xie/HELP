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

#ifndef _LP_H

//#define DEBUG_TRACK 
//#define PRINT_DEBUG
//#define PRINT_TIME

#define _LP_H
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <map>
#include <vector>

using namespace std;

// GlobalVar.cpp
const int MAX_FILENAME = 100;
extern char Filename[MAX_FILENAME];
extern char Enabled_RHS[];
extern char Enabled_RANGES[];
extern char Enabled_BOUNDS[];

// Problem Related
const int MAX_ROWS = 10000;
const int MAX_COLS = 10000;
const int MAX_ELEMENTS = 10000000;
const int MAX_PROBLEM_NAME = 100;
const double Input_Tolerance = 1e-8;
const double Variable_Tolerance = 1e-12;
const double MaxPositive = 1e+30;
const double MaxFinite = 0.99995 * MaxPositive;
const double Var_Lower_Bound = 0;
const double Var_Upper_Bound = MaxPositive;

// Homogeneous Algorithm Parameter
const double STEPSIZE_GAMMA = 0.99995;
const double STEPSIZE_BETA = 1e-8;
const int Max_Iterations = 10000;
const double Mu_Tolerance = 1e-12;
const double Infeasibility_Tolerance = 1e-8;
const double Gap_Tolerance = 1e-8;
const double Primal_Infeasibility_Tolerance = 1e-8;
const double Dual_Infeasibility_Tolerance = 1e-8;

extern char Problem_Name[MAX_PROBLEM_NAME];
extern int n_Row, n_Col, n_Element;
extern char Row_Type[MAX_ROWS];
extern double V_Cost[MAX_COLS]; // c, Cost Row
extern double V_RHS[MAX_ROWS]; // b, RHS
extern double V_RHS_r[MAX_ROWS]; // RANGES, if Row_Type[i] == 'R', then [V_RHS_r[i], V_RHS[i]]
extern double V_LB[MAX_COLS], V_UB[MAX_COLS]; // For Variable x[j], V_LB[j] <= x[j] <= V_UB[j];
// A, Column Majored Matrix, Linked List
extern long V_Matrix_Head[MAX_COLS], V_Matrix_Next[MAX_ELEMENTS], V_Matrix_Row[MAX_ELEMENTS];
extern double V_Matrix_Value[MAX_ELEMENTS];

// Crushing
extern double V_Cost_Intercept; // After crushing, objective may have nonzero intercept
extern int V_Crushing_Times[MAX_COLS];
extern double V_Crushing_Add[MAX_COLS]; // Output (x[i] * V_Crushing_Times[i] + V_Crushing_Add[i])

// Presolve
extern int n_LB, n_UB, n_FR; 
// Rearranged as: (a) x[0 ~ (n_UB - 1)]: With LB and UB; 
//                (b) x[n_UB ~ (n_LB - 1)]: With LB only; 
//                (c) x[n_LB ~ (n_LB + n_FR - 1)]: Free.
extern int RecoverOrder[MAX_COLS], TransOrder[MAX_COLS];

// Crushing.cpp
int CRUSH_Main();

// Helper.cpp
void CheckError(int ExitID, char* ErrMsg);
clock_t GetTime();
double DotProduct(int n, double* a, double* b);
void SetScaledVector(int n, double alpha, double* src, double* dest); // dest = alpha * src
int CHOLMOD_Construct();
int CHOLMOD_Destruct();
void SetATimesVector(int Transpose, double alpha, double beta, double* v, double* dest);
void RenewCholesky(double* d);
int SolveLinearEquation(double* d, double* b_1, double* b_2, double* x_1, double* x_2);

// HSDSolver.cpp
const int HSD_STATUS_OK = 0;
const int HSD_STATUS_STALLED = 1;
const int HSD_STATUS_ITER_LIMIT = 2;
const int HSD_STATUS_UNKNOWN_ERROR = 3;
int HSD_Init();
int HSD_Main();

// Init.cpp
int Prog_Init();

// MPSRead.cpp
int MPS_ReadFile();

// Presolve.cpp
int Presolve_Main();

#endif
