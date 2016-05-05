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
long V_Matrix_Head[MAX_COLS], V_Matrix_Next[MAX_ELEMENTS], V_Matrix_Row[MAX_ELEMENTS];
double V_Matrix_Value[MAX_ELEMENTS];

// Crushing
double V_Cost_Intercept; // After crushing, objective may have nonzero intercept
int V_Crushing_Times[MAX_COLS];
double V_Crushing_Add[MAX_COLS]; // Output (x[i] * V_Crushing_Times[i] + V_Crushing_Add[i])
