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
#include "MKL_Util.h"
#include "PCG.h"

int nnzA; 
double *csrValAt; int csrRowPtrAt[MAX_COLS + 1]; int *csrColIndAt;
double *csrValA; int csrRowPtrA[MAX_ROWS + 1]; int *csrColIndA;

double tmp_row[MAX_ROWS], tmp_col[MAX_COLS];
double X1[MAX_COLS], X2[MAX_ROWS];

int LinearEquation_Construct()
{
	//mkl_set_num_threads(4);
	// Count A size
	nnzA = 0;
	for (int j = 0; j < n_Col; j ++)
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
			nnzA ++;
	csrValAt = (double*) malloc(sizeof(double) * nnzA);
	csrColIndAt = (int*) malloc(sizeof(double) * nnzA);
	csrValA = (double*) malloc(sizeof(double) * nnzA);
	csrColIndA = (int*) malloc(sizeof(double) * nnzA);

	// Construct matrix At in CSR format. 
	// After Presolve_Init(), we can assume A is sorted.
	nnzA = 0;
	for (int j = 0; j < n_Col; j ++)
	{
		csrRowPtrAt[j] = nnzA;
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
		{
			csrValAt[nnzA] = V_Matrix_Value[p];
			csrColIndAt[nnzA] = V_Matrix_Row[p]; // One-based
			nnzA ++;
		}
	}
	csrRowPtrAt[n_Col] = nnzA;

	// Construct matrix A in CSR format. 
	nnzA = 0;
	for (int i = 0; i < n_Row; i ++)
	{
		csrRowPtrA[i] = nnzA;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			csrValA[nnzA] = V_Matrix_Value[p];
			csrColIndA[nnzA] = V_Matrix_Col[p]; // One-based
			nnzA ++;
		}
	}
	csrRowPtrA[n_Row] = nnzA;
	return 0;
}

int LinearEquation_Destruct()
{
	if (csrValAt) { free(csrValAt); csrValAt = 0; }
	if (csrColIndAt) { free(csrColIndAt); csrColIndAt = 0; }
	return 0;
}

// Renew Preconditioner
void RenewLinearEquation(double* d) // d should be inversed
{
}

double dinv[MAX_COLS], CG_WorkVar[MAX_COLS + MAX_ROWS * 4];

// Solve [ D  A^T ][x_1] == [b_1]
//       [ A   0  ][x_2]    [b_2]
// A: m * n; D: Diagonal Matrix, n * n; x_1: n * 1; x_2: m * 1; b_1: n * 1; b_2: m * 1
// Solution: (A D^(-1) A^T) x_2 = (A D^(-1) b_1 - b_2)
//           x_1 = D^(-1) (b_1 - A^T x_2)
// Need to call RenewLinearEquation(d) if ADA^T is not factorized yet!

int SolveLinearEquation(double* d, double* b_1, double* b_2, double* x_1, double* x_2)
{
	double lambda = 1e-10;
	double delta = 1e-10;
#ifdef PRINT_DEBUG
	printf("b_1 = zeros(%d, 1);\n", n_Col);
	for (int i = 0; i < n_Col; i ++)
		printf("b_1(%d) = %.6lf;\n", i + 1, b_1[i]);
	printf("b_2 = zeros(%d, 1);\n", n_Row);
	for (int i = 0; i < n_Row; i ++)
		printf("b_2(%d) = %.6lf;\n", i + 1, b_2[i]);
#endif
	for (int i = 0; i < n_Row; i ++)
		tmp_row[i] = -b_2[i];
	// D^(-1) b_1
	for (int i = 0; i < n_Col; i ++)
		tmp_col[i] = b_1[i] / d[i];
	// A * (D^(-1) b_1) - b_2
	char MKL_matdescra[6] = {0};
	MKL_matdescra[0] = 'G';
	MKL_matdescra[3] = 'C';
	mkl_dcsrmv(&CHAR_T, &n_Col, &n_Row, &DOUBLE_ONE, MKL_matdescra, csrValAt, csrColIndAt, csrRowPtrAt, csrRowPtrAt + 1, tmp_col, &DOUBLE_ONE, tmp_row);

	// Solve (A D^(-1) A^T) x_2 = (A D^(-1) b_1 - b_2)
#ifdef PRINT_TIME
	double Tm = GetTime();
#endif
#ifdef DEBUG_TRACK
printf("%%Before Solve\n");
#endif

/*
	FILE* out = fopen("debug.txt", "w");
	fprintf(out, "n_Row = %d;\n", n_Row);
	fprintf(out, "n_Col = %d;\n", n_Col);
	fprintf(out, "A = zeros(n_Row, n_Col);\n");
	for (int j = 0; j < n_Col; j ++)
		for (int p = V_Matrix_Col_Head[j]; p != -1; p = V_Matrix_Col_Next[p])
			fprintf(out, "A(%d, %d) = %lf;\n", V_Matrix_Row[p] + 1, j + 1, V_Matrix_Value[p]);
	fprintf(out, "Dinv = zeros(n_Col, 1);\n");
	for (int j = 0; j < n_Col; j ++)
		fprintf(out, "Dinv(%d) = %lf;\n", j + 1, dinv[j]);
	fprintf(out, "Dinv = diag(Dinv);\n");
	fprintf(out, "rhs = zeros(n_Row, 1);\n");
	for (int j = 0; j < n_Row; j ++)
		fprintf(out, "rhs(%d) = %lf;\n", j + 1, tmp_row[j]);
	
	fprintf(out, "csrRowPtrAt = [\n");
	for (int i = 0; i < n_Col; i ++)
		fprintf(out, "%d, ", csrRowPtrAt[i]);
	fprintf(out, "]\n");
	fprintf(out, "csrColIndAt = [\n");
	for (int i = 0; i < nnzA; i ++)
		fprintf(out, "%d, ", csrColIndAt[i]);
	fprintf(out, "]\n");
	fprintf(out, "csrValAt = [\n");
	for (int i = 0; i < nnzA; i ++)
		fprintf(out, "%lf, ", csrValAt[i]);
	fprintf(out, "]\n");

	fclose(out);
*/

	double CG_gamma = 0;
	double CG_delta = 0;
	ConjugateGradient(CG_gamma, CG_delta, n_Row, n_Col, csrValAt, csrColIndAt, csrRowPtrAt, csrValA, csrColIndA, csrRowPtrA, d, tmp_row, X2, 
		CG_WorkVar, CG_WorkVar + n_Col, CG_WorkVar + n_Col + n_Row, CG_WorkVar + n_Col + n_Row * 2, CG_WorkVar + n_Col + n_Row * 3);

#ifdef DEBUG_TRACK
printf("%%After Solve\n");
#endif
#ifdef PRINT_TIME
	printf("%%Solve: %.2lf s\n", GetTime() - Tm);
#endif
	for (int i = 0; i < n_Row; i ++)
		x_2[i] = X2[i];
#ifdef PRINT_DEBUG
	for (int i = 0; i < n_Row; i ++)
		printf("x_2(%d) = %lf\n", i + 1, x_2[i]);
#endif
	// b_1 - A^T x_2
	for (int i = 0; i < n_Col; i ++)
		tmp_col[i] = b_1[i];
	mkl_dcsrmv(&CHAR_N, &n_Col, &n_Row, &DOUBLE_NEGONE, MKL_matdescra, csrValAt, csrColIndAt, csrRowPtrAt, csrRowPtrAt + 1, x_2, &DOUBLE_ONE, tmp_col);
	// Let x_1 = D^(-1) (b_1 - A^T x_2)
	for (int i = 0; i < n_Col; i ++)
		x_1[i] = tmp_col[i] / d[i];
#ifdef PRINT_DEBUG
	for (int i = 0; i < n_Col; i ++)
		printf("x_1(%d) = %lf\n", i + 1, x_1[i]);
#endif
	return 0;
}


// Regularized KKT system solves
// [ D + lambda I    A^T  ][x_1] == [b_1]
// [      A        delta I][x_2]    [b_2]
// Solution: (A (D + lambda I)^(-1) A^T - delta I) x_2 = A (D + lambda I)^(-1) b_1 - b_2
//           x_1 = (D + lambda I)^(-1) (b_1 - A^T x_2)