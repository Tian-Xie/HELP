#include "LP.h"

char Filename[MAX_FILENAME];
char Enabled_RHS[] = "";
char Enabled_RANGES[] = "";
char Enabled_BOUNDS[] = "";

// Problem Related
char Problem_Name[MAX_PROBLEM_NAME];
int n_Row, n_Col, n_Element;
vector <char> Row_Type;
vector <double> V_Cost; // c, Cost Row
vector < map <int, double> > V_Matrix; // A, Column Majored Matrix
vector <double> V_RHS; // b, RHS
vector <double> V_RHS_r; // RANGES, if Row_Type[i] == 'R', then [V_RHS_r[i], V_RHS[i]]
vector <double> V_LB, V_UB; // For Variable x[j], V_LB[j] <= x[j] <= V_UB[j];
