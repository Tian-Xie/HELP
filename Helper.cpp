#include "LP.h"
// Using CHOLMOD library, by Timothy A. Davis
#include "cholmod.h"

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

void SetATimesVector(int Transpose, int Sign, double* v, double* dest) // dest = dest + Sign * A * v, or dest = dest + Sign * A^T * v, Sign \in {1, -1}
{
	for (int Col = 0; Col < n_Col; Col ++)
	{
		for (int j = V_Matrix_Head[Col]; j != -1; j = V_Matrix_Next[j])
		{
			int Row = V_Matrix_Row[j];
			int Value = V_Matrix_Value[j];
			if (! Transpose) // A * v
			{
				if (Sign == 1)
					dest[Row] += Value * v[Col];
				else
					dest[Row] -= Value * v[Col];
			}
			else // A^T * v
			{
				if (Sign == 1)
					dest[Col] += Value * v[Row];
				else
					dest[Col] -= Value * v[Row];
			}
		}
	}
}

cholmod_sparse* CHOL_A;
cholmod_common CHOL_Com;
cholmod_dense* CHOL_Vector_Row;
cholmod_dense* CHOL_Vector_Col;

int CHOLMOD_Construct()
{
	cholmod_start(&CHOL_Com); // Start CHOLMOD

	cholmod_triplet* Triplet = cholmod_allocate_triplet(n_Row, n_Col, n_Element, 0, CHOLMOD_REAL, &CHOL_Com); // stype == 0, asymmetric
	int k = 0;
	Triplet -> nnz = n_Element;
	printf("A = zeros(%d, %d);\n", n_Row, n_Col);
	for (int Col = 0; Col < n_Col; Col ++)
		for (int j = V_Matrix_Head[Col]; j != -1; j = V_Matrix_Next[j])
		{
			((int*) Triplet -> i)[k] = V_Matrix_Row[j];
			((int*) Triplet -> j)[k] = Col;
			((double*) Triplet -> x)[k] = V_Matrix_Value[j];
			printf("A(%d, %d) = %f;\n", ((int*) Triplet -> i)[k] + 1, ((int*) Triplet -> j)[k] + 1, ((double*) Triplet -> x)[k]);
			k ++;
		}
	CHOL_A = cholmod_triplet_to_sparse(Triplet, n_Element, &CHOL_Com);
	cholmod_free_triplet(&Triplet, &CHOL_Com);
	cholmod_sort(CHOL_A, &CHOL_Com);

	CHOL_Vector_Row = cholmod_allocate_dense(n_Row, 1, n_Row, CHOLMOD_REAL, &CHOL_Com);
	CHOL_Vector_Col = cholmod_allocate_dense(n_Col, 1, n_Col, CHOLMOD_REAL, &CHOL_Com);
	return 0;
}

int CHOLMOD_Destruct()
{
	cholmod_free_sparse(&CHOL_A, &CHOL_Com);
	cholmod_free_dense(&CHOL_Vector_Row, &CHOL_Com);
	cholmod_free_dense(&CHOL_Vector_Col, &CHOL_Com);
	cholmod_finish(&CHOL_Com); // Finish CHOLMOD
	return 0;
}

// Solve [ D  A^T ][x_1] == [b_1]
//       [ A   0  ][x_2]    [b_2]
// A: m * n; D: Diagonal Matrix, n * n; x_1: n * 1; x_2: m * 1; b_1: n * 1; b_2: m * 1
// Solution: (A D^(-1) A^T) x_2 = (A D^(-1) b_1 - b_2)
//           x_1 = D^(-1) (b_1 - A^T x_2)
int SolveLinearEquation(double* d, double* b_1, double* b_2, double* x_1, double* x_2)
{
	double SDMULT_POSITIVE[] = {1, 0};
	double SDMULT_NEGATIVE[] = {-1, 0};
	
	// D^(-1/2) A
	for (int i = 0; i < n_Col; i ++)
		((double*) CHOL_Vector_Col -> x)[i] = 1.0 / sqrt(d[i]);
	cholmod_sparse* DsqrtinvA = cholmod_copy_sparse(CHOL_A, &CHOL_Com);
	cholmod_scale(CHOL_Vector_Col, CHOLMOD_COL, DsqrtinvA, &CHOL_Com);
	
	for (int i = 0; i < n_Row; i ++)
		((double*) CHOL_Vector_Row -> x)[i] = -b_2[i];
	// D^(-1) b_1
	for (int i = 0; i < n_Col; i ++)
		((double*) CHOL_Vector_Col -> x)[i] = b_1[i] / d[i];
	// A * (D^(-1) b_1) - b_2
	cholmod_sdmult(CHOL_A, 0, SDMULT_POSITIVE, SDMULT_POSITIVE, CHOL_Vector_Col, CHOL_Vector_Row, &CHOL_Com);

	// Solve (A D^(-1) A^T) x_2 = (A D^(-1) b_1 - b_2)
	cholmod_factor* Factor = cholmod_analyze(DsqrtinvA, &CHOL_Com);	
	cholmod_factorize(DsqrtinvA, Factor, &CHOL_Com);
	cholmod_free_sparse(&DsqrtinvA, &CHOL_Com);
	cholmod_dense* X2 = cholmod_solve(CHOLMOD_LDLt, Factor, CHOL_Vector_Row, &CHOL_Com);
	cholmod_free_factor(&Factor, &CHOL_Com);
	for (int i = 0; i < n_Col; i ++)
		x_2[i] = ((double*) X2 -> x)[i];

	// b_1 - A^T x_2
	for (int i = 0; i < n_Col; i ++)
		((double*) CHOL_Vector_Col -> x)[i] = b_1[i];
	cholmod_sdmult(CHOL_A, 1, SDMULT_NEGATIVE, SDMULT_POSITIVE, X2, CHOL_Vector_Col, &CHOL_Com);
	// Let x_1 = D^(-1) (b_1 - A^T x_2)
	for (int i = 0; i < n_Col; i ++)
		x_1[i] = ((double*) CHOL_Vector_Col -> x)[i] / d[i];
	cholmod_free_dense(&X2, &CHOL_Com);
	return 0;
}
