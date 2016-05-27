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
#include "UseCholmod.h"

void CheckError(int ExitID, char* ErrMsg)
{
	if (ExitID)
	{
		printf("ERROR: Msg = %s\n", ErrMsg);
		exit(ExitID);
	}
}

clock_t GetTime()
{
	return clock();
}

double DotProduct(int n, double* a, double* b)
{
	double Ret = 0;
	for (int i = 0; i < n; i ++)
		Ret += a[i] * b[i];
	return Ret;
}

void SetScaledVector(int n, double alpha, double* src, double* dest) // dest = alpha * src
{
	for (int i = 0; i < n; i ++)
	{
		if (alpha == 1)
			dest[i] = src[i];
		else if (alpha == -1)
			dest[i] = -src[i];
		else if (alpha == 0)
			dest[i] = 0;
		else
			dest[i] = alpha * src[i];
	}
}

cholmod_sparse* CHOL_A;
cholmod_common CHOL_Com;
cholmod_dense* CHOL_Vector_Row;
cholmod_dense* CHOL_Vector_Col;
cholmod_dense* CHOL_Vector_Row_Reorder;
cholmod_dense* CHOL_Vector_Col_Reorder;
cholmod_dense* CHOL_TIMES_Row;
cholmod_dense* CHOL_TIMES_Col;
cholmod_factor* CHOL_Fac;

int CHOLMOD_Construct()
{
printf("SuiteSparse_long: %d\n", sizeof(SuiteSparse_long));
printf("long: %d\n", sizeof(long));
	cholmod_l_start(&CHOL_Com); // Start CHOLMOD
	// Initially open all options
	CHOL_Com.nmethods = 9;
	CHOLMOD_Setting(&CHOL_Com); // In Setting_CPU.cpp / Setting_CPU.cpp

	cholmod_triplet* Triplet = cholmod_l_allocate_triplet(n_Row, n_Col, n_Element, 0, CHOLMOD_REAL, &CHOL_Com); // stype == 0, asymmetric
	int k = 0;
	Triplet -> nnz = n_Element;
#ifdef PRINT_DEBUG
	printf("A = zeros(%d, %d);\n", n_Row, n_Col);
#endif
	for (int Col = 0; Col < n_Col; Col ++)
		for (int j = V_Matrix_Head[Col]; j != -1; j = V_Matrix_Next[j])
		{
			((long*) Triplet -> i)[k] = V_Matrix_Row[j];
			((long*) Triplet -> j)[k] = Col;
			((double*) Triplet -> x)[k] = V_Matrix_Value[j];
#ifdef PRINT_DEBUG
			printf("A(%d, %d) = %f;\n", ((long*) Triplet -> i)[k] + 1, ((long*) Triplet -> j)[k] + 1, ((double*) Triplet -> x)[k]);
#endif
			k ++;
		}
	CHOL_A = cholmod_l_triplet_to_sparse(Triplet, n_Element, &CHOL_Com);
	cholmod_l_free_triplet(&Triplet, &CHOL_Com);
	cholmod_l_sort(CHOL_A, &CHOL_Com);

	CHOL_Vector_Row = cholmod_l_allocate_dense(n_Row, 1, n_Row, CHOLMOD_REAL, &CHOL_Com);
	CHOL_Vector_Col = cholmod_l_allocate_dense(n_Col, 1, n_Col, CHOLMOD_REAL, &CHOL_Com);
	CHOL_Vector_Row_Reorder = cholmod_l_allocate_dense(n_Row, 1, n_Row, CHOLMOD_REAL, &CHOL_Com);
	CHOL_Vector_Col_Reorder = cholmod_l_allocate_dense(n_Col, 1, n_Col, CHOLMOD_REAL, &CHOL_Com);
	CHOL_TIMES_Row = cholmod_l_allocate_dense(n_Row, 1, n_Row, CHOLMOD_REAL, &CHOL_Com);
	CHOL_TIMES_Col = cholmod_l_allocate_dense(n_Col, 1, n_Col, CHOLMOD_REAL, &CHOL_Com);
	CHOL_Fac = NULL;
	return 0;
}

int CHOLMOD_Destruct()
{
	cholmod_l_free_sparse(&CHOL_A, &CHOL_Com);
	cholmod_l_free_dense(&CHOL_Vector_Row, &CHOL_Com);
	cholmod_l_free_dense(&CHOL_Vector_Col, &CHOL_Com);
	cholmod_l_free_dense(&CHOL_Vector_Row_Reorder, &CHOL_Com);
	cholmod_l_free_dense(&CHOL_Vector_Col_Reorder, &CHOL_Com);
	cholmod_l_free_dense(&CHOL_TIMES_Row, &CHOL_Com);
	cholmod_l_free_dense(&CHOL_TIMES_Col, &CHOL_Com);
	if (CHOL_Fac)
	{
		cholmod_l_free_factor(&CHOL_Fac, &CHOL_Com);
		CHOL_Fac = NULL;
	}
	cholmod_l_finish(&CHOL_Com); // Finish CHOLMOD
	return 0;
}

// Transpose = 0: dest = alpha * A * v + beta * dest
// Transpose = 1: dest = alpha * A^T * v + beta * dest
void SetATimesVector(int Transpose, double alpha, double beta, double* v, double* dest)
{
	double SDMULT_ALPHA[] = {alpha, 0};
	double SDMULT_BETA[] = {beta, 0};
	if (Transpose)
	{
		for (int i = 0; i < n_Row; i ++)
			((double*) CHOL_TIMES_Row -> x)[i] = v[i];
		for (int i = 0; i < n_Col; i ++)
			((double*) CHOL_TIMES_Col -> x)[i] = dest[i];
		cholmod_l_sdmult(CHOL_A, Transpose, SDMULT_ALPHA, SDMULT_BETA, CHOL_TIMES_Row, CHOL_TIMES_Col, &CHOL_Com);
		for (int i = 0; i < n_Col; i ++)
			dest[i] = ((double*) CHOL_TIMES_Col -> x)[i];
	}
	else
	{
		for (int i = 0; i < n_Col; i ++)
			((double*) CHOL_TIMES_Col -> x)[i] = v[i];
		for (int i = 0; i < n_Row; i ++)
			((double*) CHOL_TIMES_Row -> x)[i] = dest[i];
		cholmod_l_sdmult(CHOL_A, Transpose, SDMULT_ALPHA, SDMULT_BETA, CHOL_TIMES_Col, CHOL_TIMES_Row, &CHOL_Com);
		for (int i = 0; i < n_Row; i ++)
			dest[i] = ((double*) CHOL_TIMES_Row -> x)[i];
	}
}

// Renew CHOL_Fac: ADA^T
void RenewCholesky(double* d)
{
#ifdef DEBUG_TRACK
printf("In RenewCholesky\n");
#endif
	if (CHOL_Fac)
	{
		cholmod_l_free_factor(&CHOL_Fac, &CHOL_Com);
		CHOL_Fac = NULL;
	}
	// D^(-1/2) A
	for (int i = 0; i < n_Col; i ++)
		((double*) CHOL_Vector_Col -> x)[i] = 1.0 / sqrt(d[i]);
	cholmod_sparse* DsqrtinvA = cholmod_l_copy_sparse(CHOL_A, &CHOL_Com);
	cholmod_l_scale(CHOL_Vector_Col, CHOLMOD_COL, DsqrtinvA, &CHOL_Com);

#ifdef PRINT_TIME
	clock_t Tm;
	Tm = GetTime();
#endif
#ifdef DEBUG_TRACK
printf("RenewCholesky: Before Analyze\n");
#endif
    CHOL_Fac = cholmod_l_analyze(DsqrtinvA, &CHOL_Com);
#ifdef DEBUG_TRACK
printf("RenewCholesky: After Analyze\n");
#endif
#ifdef PRINT_TIME
	printf("Analyze: %d ms\n", GetTime() - Tm);
#endif
	if (CHOL_Com.nmethods == 9) // Choose a best method
	{
		CHOL_Com.nmethods = 1;
		CHOL_Com.method[0] = CHOL_Com.method[CHOL_Com.selected];
#ifdef DEBUG_TRACK
printf("RenewCholesky: Best Method = %d\n", CHOL_Com.selected);
#endif
	}
#ifdef PRINT_TIME
	Tm = GetTime();
#endif
#ifdef DEBUG_TRACK
printf("RenewCholesky: Before Factorize\n");
#endif
	cholmod_l_factorize(DsqrtinvA, CHOL_Fac, &CHOL_Com);
#ifdef DEBUG_TRACK
printf("RenewCholesky: After Factorize\n");
#endif
#ifdef PRINT_TIME
	printf("Factorize: %d ms\n", GetTime() - Tm);
#endif
	cholmod_l_free_sparse(&DsqrtinvA, &CHOL_Com);
#ifdef DEBUG_TRACK
printf("Out RenewCholesky\n");
#endif
}

// Solve [ D  A^T ][x_1] == [b_1]
//       [ A   0  ][x_2]    [b_2]
// A: m * n; D: Diagonal Matrix, n * n; x_1: n * 1; x_2: m * 1; b_1: n * 1; b_2: m * 1
// Solution: (A D^(-1) A^T) x_2 = (A D^(-1) b_1 - b_2)
//           x_1 = D^(-1) (b_1 - A^T x_2)
// Need to call RenewCholesky(d) if ADA^T is not factorized yet!
int SolveLinearEquation(double* d, double* b_1, double* b_2, double* x_1, double* x_2)
{
#ifdef DEBUG_TRACK
printf("In SolveLinearEquation\n");
#endif

	double SDMULT_POSITIVE[] = {1, 0};
	double SDMULT_NEGATIVE[] = {-1, 0};
	
#ifdef PRINT_DEBUG
	printf("b_1 = zeros(%d, 1);\n", n_Col);
	for (int i = 0; i < n_Col; i ++)
		printf("b_1(%d) = %.6lf;\n", i + 1, b_1[i]);
	printf("b_2 = zeros(%d, 1);\n", n_Row);
	for (int i = 0; i < n_Row; i ++)
		printf("b_2(%d) = %.6lf;\n", i + 1, b_2[i]);
#endif

	for (int i = 0; i < n_Row; i ++)
		((double*) CHOL_Vector_Row -> x)[i] = -b_2[i];
	// D^(-1) b_1
	for (int i = 0; i < n_Col; i ++)
		((double*) CHOL_Vector_Col -> x)[i] = b_1[i] / d[i];
	// A * (D^(-1) b_1) - b_2
	cholmod_l_sdmult(CHOL_A, 0, SDMULT_POSITIVE, SDMULT_POSITIVE, CHOL_Vector_Col, CHOL_Vector_Row, &CHOL_Com);
	for (int i = 0; i < n_Row; i ++)
	{
		((double*) CHOL_Vector_Row_Reorder -> x)[i] = ((double*) CHOL_Vector_Row -> x)[((long *) CHOL_Fac -> Perm) [i]];
#ifdef PRINT_DEBUG
		printf("%d\n", ((int *) CHOL_Fac -> Perm) [i]);
#endif
	}

	// Solve (A D^(-1) A^T) x_2 = (A D^(-1) b_1 - b_2)
#ifdef PRINT_TIME
	clock_t Tm = GetTime();
#endif
	cholmod_dense* X2;
	cholmod_dense* Tmp;
	/*if (CHOL_Fac -> is_ll)
	{
		Tmp = cholmod_l_solve(CHOLMOD_L, CHOL_Fac, CHOL_Vector_Row_Reorder, &CHOL_Com);
		X2 = cholmod_l_solve(CHOLMOD_Lt, CHOL_Fac, Tmp, &CHOL_Com);
		cholmod_l_free_dense(&Tmp, &CHOL_Com);
	}
	else*/
#ifdef DEBUG_TRACK
printf("SolveLinearEquation: Before Solve\n");
#endif
		X2 = cholmod_l_solve(CHOLMOD_LDLt, CHOL_Fac, CHOL_Vector_Row_Reorder, &CHOL_Com);
#ifdef DEBUG_TRACK
printf("SolveLinearEquation: After Solve\n");
#endif

#ifdef PRINT_TIME
	printf("Solve: %d ms\n", GetTime() - Tm);
#endif
	for (int i = 0; i < n_Row; i ++)
		x_2[((long *) CHOL_Fac -> Perm) [i]] = ((double*) X2 -> x)[i];
	for (int i = 0; i < n_Row; i ++)
		((double*) X2 -> x)[i] = x_2[i];
#ifdef PRINT_DEBUG
	for (int i = 0; i < n_Row; i ++)
		printf("x_2(%d) = %lf\n", i + 1, x_2[i]);
#endif

	// b_1 - A^T x_2
	for (int i = 0; i < n_Col; i ++)
		((double*) CHOL_Vector_Col -> x)[i] = b_1[i];
	cholmod_l_sdmult(CHOL_A, 1, SDMULT_NEGATIVE, SDMULT_POSITIVE, X2, CHOL_Vector_Col, &CHOL_Com);
	// Let x_1 = D^(-1) (b_1 - A^T x_2)
	for (int i = 0; i < n_Col; i ++)
		x_1[i] = ((double*) CHOL_Vector_Col -> x)[i] / d[i];
	cholmod_l_free_dense(&X2, &CHOL_Com);
	
#ifdef DEBUG_TRACK
printf("Out SolveLinearEquation\n");
#endif
	return 0;
}
