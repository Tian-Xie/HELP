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

const int MAX_FILENAME = 100;
char Filename[MAX_FILENAME];
char Enabled_RHS[] = "";
char Enabled_RANGES[] = "";
char Enabled_BOUNDS[] = "";

// Problem Related
const int MAX_PROBLEM_NAME = 100;
const double Input_Tolerance = 1e-8;
const double MaxPositive = 1e+30;
const double Var_Lower_Bound = 0;
const double Var_Upper_Bound = MaxPositive;

char Problem_Name[MAX_PROBLEM_NAME];
int n_Row, n_Col, n_Element;
vector <char> Row_Type;
vector <double> V_Cost; // c, Cost Row
vector < map <int, double> > V_Matrix; // A, Column Majored Matrix
vector <double> V_RHS; // b, RHS
vector <double> V_RHS_r; // RANGES, if Row_Type[i] == 'R', then [V_RHS_r[i], V_RHS[i]]
vector <double> V_LB, V_UB; // For Variable x[j], V_LB[j] <= x[j] <= V_UB[j];

// Helper.cpp
void CheckError(int ExitID, char* ErrMsg);

// Init.cpp
int Prog_Init();

// MPSRead.cpp
int MPS_ReadFile();

#endif

