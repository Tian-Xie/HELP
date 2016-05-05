#ifndef _LP_H

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
const double MaxPositive = 1e+30;
const double MaxFinite = 0.99995 * MaxPositive;
const double Var_Lower_Bound = 0;
const double Var_Upper_Bound = MaxPositive;

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

// Crushing.cpp
int CRUSH_Main();

// Helper.cpp
void CheckError(int ExitID, char* ErrMsg);
clock_t GetTime();

// Init.cpp
int Prog_Init();

// MPSRead.cpp
int MPS_ReadFile();

#endif

