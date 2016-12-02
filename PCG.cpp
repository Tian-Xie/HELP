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

#include "MKL_Util.h"
#include "PCG.h"
#include <cstdio>

int PCG_PERM[MAX_ROWS];

double DIAG_S[MAX_ROWS];
double TEMP_GR[MAX_ROWS];
double PARTIAL_CHOL_L[MAX_ROWS][MAX_PART_CHOL_COLS]; // TODO: Exchange two dimensions to achieve a better performance!
double PARTIAL_CHOL_U[MAX_ROWS][MAX_PART_CHOL_COLS];
double PARTIAL_CHOL_D[MAX_ROWS];
double TEMP_SCATTER[MAX_COLS];

int csrRowPtrA_BEG[MAX_ROWS], csrRowPtrA_END[MAX_ROWS];

double L11[MAX_PART_CHOL_COLS][MAX_PART_CHOL_COLS];
double L11Inv[MAX_PART_CHOL_COLS][MAX_PART_CHOL_COLS];
double csrValPARTIAL_CHOL_INV1[MAX_ROWS * MAX_PART_CHOL_COLS + MAX_ROWS];
int csrColIndPARTIAL_CHOL_INV1[MAX_ROWS * MAX_PART_CHOL_COLS + MAX_ROWS], csrRowPtrPARTIAL_CHOL_INV1[MAX_ROWS + 1];
double csrValPARTIAL_CHOL_INV2[MAX_PART_CHOL_COLS * MAX_PART_CHOL_COLS + MAX_ROWS];
int csrColIndPARTIAL_CHOL_INV2[MAX_PART_CHOL_COLS * MAX_PART_CHOL_COLS + MAX_ROWS], csrRowPtrPARTIAL_CHOL_INV2[MAX_ROWS + 1];

// Renew Partial Cholesky of G_R = A (D + gamma I)^(-1) A^T + delta I
void RenewPartialCholesky(int n_Row, int n_Col, double* csrValA, int* csrColIndA, int* csrRowPtrA, double* d, double gamma, double delta)
{
	// Step 1: All Diagonal Elements => DIAG_S
	for (int i = 0; i < n_Row; i ++)
	{
		DIAG_S[i] = delta;
		for (int j = csrRowPtrA[i]; j < csrRowPtrA[i + 1]; j ++)
			DIAG_S[i] += csrValA[j] * csrValA[j] / (d[csrColIndA[j]] + gamma);
	}

	// Step 2: Partial Cholesky Decomposition
	for (int i = 0; i < n_Row; i ++)
		PCG_PERM[i] = i;
	for (int k = 0; k < PART_CHOL_COLS && k < n_Row; k ++)
	{
		// Choose pivoting Column col_max with maximum DIAG_S[j]
		int col_max_ind = k;
		for (int i = k + 1; i < n_Row; i ++)
			if (DIAG_S[PCG_PERM[i]] > DIAG_S[PCG_PERM[col_max_ind]])
				col_max_ind = i;
		// Bubble the corresponding column up
		swap(PCG_PERM[k], PCG_PERM[col_max_ind]);
		int col_max = PCG_PERM[k];
		// Now, deal with Column col_max

		// Calculate the col_max-th column of GR => TEMP_GR
		// (1) Scatter A(col_max, :) => TEMP_SCATTER
		for (int j = csrRowPtrA[col_max]; j < csrRowPtrA[col_max + 1]; j ++)
			TEMP_SCATTER[csrColIndA[j]] = csrValA[j] / (d[csrColIndA[j]] + gamma);
		// (2) Calculate Inner Products
		for (int i = k; i < n_Row; i ++)
		{
			int ii = PCG_PERM[i];
			// Calculate GR(ii, col_max) => TEMP_GR(ii)
			TEMP_GR[ii] = 0;
			for (int j = csrRowPtrA[ii]; j < csrRowPtrA[ii + 1]; j ++)
				TEMP_GR[ii] += csrValA[j] * TEMP_SCATTER[csrColIndA[j]];
		}
		TEMP_GR[col_max] += delta;
		// (3) Recover TEMP_SCATTER
		for (int j = csrRowPtrA[col_max]; j < csrRowPtrA[col_max + 1]; j ++)
			TEMP_SCATTER[csrColIndA[j]] = 0;

		// Calculate Cholesky Column according to:
		// d_k = a_{kk} - sum_{r = 1, ..., (k - 1)} u_{kr} l_{kr}
		// u_{jk} = a_{jk} - sum_{r = 1, ..., (k - 1)} u_{jr} l_{kr}
		// l_{jk} = u_{jk} / d_k
		PARTIAL_CHOL_D[k] = TEMP_GR[col_max];
		for (int j = k + 1; j < n_Row; j ++)
			PARTIAL_CHOL_U[PCG_PERM[j]][k] = TEMP_GR[PCG_PERM[j]];
		for (int r = 0; r < k; r ++)
		{
			double tmp = PARTIAL_CHOL_L[col_max][r];
			PARTIAL_CHOL_D[k] -= PARTIAL_CHOL_U[col_max][r] * tmp;
			for (int j = k + 1; j < n_Row; j ++)
				PARTIAL_CHOL_U[PCG_PERM[j]][k] -= PARTIAL_CHOL_U[PCG_PERM[j]][r] * tmp;
		}
		PARTIAL_CHOL_L[col_max][k] = 1;
		for (int j = k + 1; j < n_Row; j ++)
			PARTIAL_CHOL_L[PCG_PERM[j]][k] = PARTIAL_CHOL_U[PCG_PERM[j]][k] / PARTIAL_CHOL_D[k];

		// Maintain the diagonal of Schur Complement
		double tmp_d = PARTIAL_CHOL_D[k];
		for (int j = k + 1; j < n_Row; j ++)
			DIAG_S[PCG_PERM[j]] -= tmp_d * PARTIAL_CHOL_L[PCG_PERM[j]][k] * PARTIAL_CHOL_L[PCG_PERM[j]][k];
	}

	// Step 3: Reorder A's rows
	for (int i = 0; i < n_Row; i ++)
	{
		csrRowPtrA_BEG[i] = csrRowPtrA[PCG_PERM[i]];
		csrRowPtrA_END[i] = csrRowPtrA[PCG_PERM[i] + 1];
	}

	// Step 4: Recover a Partial Cholesky Preconditioner
	// M = [ L_{11}   ][ D_L      ][ L_{11}^T  L_{21}^T ]
	//     [ L_{21} I ][      D_S ][               I    ]
	// [ L_{11}   ]^{-1} = [    I      ][ L_{11}^{-1}   ]
	// [ L_{21} I ]        [ -L_{21} I ][             I ]
	// Step 4.1: Gaussian Elimination => L_{11}^{-1}
	// Note that L11[i][j] = PARTIAL_CHOL_L[PCG_PERM[i]][j];
	for (int i = 0; i < PART_CHOL_COLS; i ++)
		for (int j = 0; j < PART_CHOL_COLS; j ++)
		{
			L11[i][j] = PARTIAL_CHOL_L[PCG_PERM[i]][j];
			L11Inv[i][j] = (double) (i == j);
		}
	// Note that the diagonal of L_{11} is all 1!
	for (int j = 0; j < PART_CHOL_COLS; j ++)
		for (int i = j + 1; i < PART_CHOL_COLS; i ++)
		{
			// Row[i] -= Row[j] * L11[i][j]
			double factor = PARTIAL_CHOL_L[PCG_PERM[i]][j];
			for (int k = 0; k <= j; k ++)
				L11Inv[i][k] -= L11Inv[j][k] * factor;
		}
	// Step 4.2: Explicitly Formulate Two Matrices
	int nnzINV1 = 0;
	for (int i = 0; i < PART_CHOL_COLS; i ++)
	{
		csrRowPtrPARTIAL_CHOL_INV1[i] = nnzINV1;
		csrColIndPARTIAL_CHOL_INV1[nnzINV1] = i; // Zero-based
		csrValPARTIAL_CHOL_INV1[nnzINV1] = 1;
		nnzINV1 ++;
	}
	for (int i = PART_CHOL_COLS; i < n_Row; i ++)
	{
		csrRowPtrPARTIAL_CHOL_INV1[i] = nnzINV1;
		for (int j = 0; j < PART_CHOL_COLS; j ++)
		{
			csrColIndPARTIAL_CHOL_INV1[nnzINV1] = j; // Zero-based
			csrValPARTIAL_CHOL_INV1[nnzINV1] = -PARTIAL_CHOL_L[PCG_PERM[i]][j];
			nnzINV1 ++;
		}
		csrColIndPARTIAL_CHOL_INV1[nnzINV1] = i; // Zero-based
		csrValPARTIAL_CHOL_INV1[nnzINV1] = 1;
		nnzINV1 ++;
	}
	csrRowPtrPARTIAL_CHOL_INV1[n_Row] = nnzINV1;
	
	int nnzINV2 = 0;
	for (int i = 0; i < PART_CHOL_COLS; i ++)
	{
		csrRowPtrPARTIAL_CHOL_INV2[i] = nnzINV2;
		for (int j = 0; j <= i; j ++)
		{
			csrColIndPARTIAL_CHOL_INV2[nnzINV2] = j; // Zero-based
			csrValPARTIAL_CHOL_INV2[nnzINV2] = L11Inv[i][j];
			nnzINV2 ++;
		}
	}
	for (int i = PART_CHOL_COLS; i < n_Row; i ++)
	{
		csrRowPtrPARTIAL_CHOL_INV2[i] = nnzINV2;
		csrColIndPARTIAL_CHOL_INV2[nnzINV2] = i; // Zero-based
		csrValPARTIAL_CHOL_INV2[nnzINV2] = 1;
		nnzINV2 ++;
	}
	csrRowPtrPARTIAL_CHOL_INV2[n_Row] = nnzINV2;
}

// Calculate ret = M^(-1) * r
// Note that ``n_Row'' and ``n_Col'' are A's dimension size!
double PREC_TMP1[MAX_ROWS], PREC_TMP2[MAX_ROWS];


void Preconditioner(int n_Row, int n_Col, double* r, double* Ret)
{
	char MKL_matdescra[6] = {0};
	MKL_matdescra[0] = 'G';
	MKL_matdescra[3] = 'C';
	
	// (LDL^T)^{-1} r = Inv2^T * Inv1^T * D^{-1} * Inv1 * Inv2 * r
	mkl_dcsrmv(&CHAR_N, &n_Row, &n_Row, &DOUBLE_ONE, MKL_matdescra, csrValPARTIAL_CHOL_INV2, csrColIndPARTIAL_CHOL_INV2, csrRowPtrPARTIAL_CHOL_INV2, csrRowPtrPARTIAL_CHOL_INV2 + 1, r, &DOUBLE_ZERO, PREC_TMP1);
	mkl_dcsrmv(&CHAR_N, &n_Row, &n_Row, &DOUBLE_ONE, MKL_matdescra, csrValPARTIAL_CHOL_INV1, csrColIndPARTIAL_CHOL_INV1, csrRowPtrPARTIAL_CHOL_INV1, csrRowPtrPARTIAL_CHOL_INV1 + 1, PREC_TMP1, &DOUBLE_ZERO, PREC_TMP2);
	for (int i = 0; i < n_Row; i ++)
		PREC_TMP2[i] /= DIAG_S[PCG_PERM[i]];
	mkl_dcsrmv(&CHAR_T, &n_Row, &n_Row, &DOUBLE_ONE, MKL_matdescra, csrValPARTIAL_CHOL_INV1, csrColIndPARTIAL_CHOL_INV1, csrRowPtrPARTIAL_CHOL_INV1, csrRowPtrPARTIAL_CHOL_INV1 + 1, PREC_TMP2, &DOUBLE_ZERO, PREC_TMP1);
	mkl_dcsrmv(&CHAR_T, &n_Row, &n_Row, &DOUBLE_ONE, MKL_matdescra, csrValPARTIAL_CHOL_INV2, csrColIndPARTIAL_CHOL_INV2, csrRowPtrPARTIAL_CHOL_INV2, csrRowPtrPARTIAL_CHOL_INV2 + 1, PREC_TMP1, &DOUBLE_ZERO, Ret);
}

/*
void Preconditioner(int n_Row, int n_Col, double* r, double* Ret)
{
	// Ret = r
	cblas_dcopy(n_Row, r, 1, Ret, 1); // y = x
}
*/

/*
double PCG_TMP_HELPER_IN[MAX_ROWS], PCG_TMP_HELPER_OUT[MAX_ROWS];

void PCGDebugHelper(int n_Row, int n_Col, double* csrValA, int* csrColIndA, double* d, double gamma, double delta)
{
	// GR = A (D + gamma I)^(-1) A^T + delta I
	FILE* O = fopen("PCGDebug.txt", "w");
	fprintf(O, "n_Row = %d\n", n_Row);
	fprintf(O, "n_Col = %d\n", n_Col);
	fprintf(O, "A = zeros(n_Row, n_Col);\n"); // Reordered
	for (int i = 0; i < n_Row; i ++)
		for (int j = csrRowPtrA_BEG[i]; j < csrRowPtrA_END[i]; j ++)
			fprintf(O, "A(%d, %d) = %lf;\n", i + 1, csrColIndA[j] + 1, csrValA[j]);
	fprintf(O, "d_diag = [");
	for (int i = 0; i < n_Col - 1; i ++)
		fprintf(O, "%lf, ", d[i]);
	fprintf(O, "%lf];\n", d[n_Col - 1]);
	fprintf(O, "gamma = %lf;\n", gamma);
	fprintf(O, "delta = %lf;\n", delta);
	fprintf(O, "GR = A * inv(diag(d_diag) + gamma * eye(n_Col)) * A' + delta * eye(n_Row);\n");

	// Partial Cholesky
	fprintf(O, "LDL_L = eye(n_Row);\n");
	for (int j = 0; j < PART_CHOL_COLS; j ++)
	{
		fprintf(O, "LDL_L(%d:%d, %d) = [", j + 1, n_Row, j + 1);
		for (int i = j; i < n_Row - 1; i ++)
			fprintf(O, "%lf, ", PARTIAL_CHOL_L[PCG_PERM[i]][j]);
		fprintf(O, "%lf];\n", PARTIAL_CHOL_L[PCG_PERM[n_Row - 1]][j]);
	}
	fprintf(O, "LDL_D = diag([");
	for (int i = 0; i < n_Row - 1; i ++)
		fprintf(O, "%lf, ", DIAG_S[PCG_PERM[i]]);
	fprintf(O, "%lf]);\n", DIAG_S[PCG_PERM[n_Row - 1]]);

	// Inverse
	fprintf(O, "Inv1 = zeros(n_Row, n_Row);\n");
	for (int i = 0; i < n_Row; i ++)
		for (int j = csrRowPtrPARTIAL_CHOL_INV1[i]; j < csrRowPtrPARTIAL_CHOL_INV1[i + 1]; j ++)
			fprintf(O, "Inv1(%d, %d) = %lf;\n", i + 1, csrColIndPARTIAL_CHOL_INV1[j] + 1, csrValPARTIAL_CHOL_INV1[j]);
	fprintf(O, "Inv2 = zeros(n_Row, n_Row);\n");
	for (int i = 0; i < n_Row; i ++)
		for (int j = csrRowPtrPARTIAL_CHOL_INV2[i]; j < csrRowPtrPARTIAL_CHOL_INV2[i + 1]; j ++)
			fprintf(O, "Inv2(%d, %d) = %lf;\n", i + 1, csrColIndPARTIAL_CHOL_INV2[j] + 1, csrValPARTIAL_CHOL_INV2[j]);
	fprintf(O, "Inv = Inv1 * Inv2;\n");

	// Try to apply a preconditioner
	fprintf(O, "tmpV = (1 : n_Row)';\n");
	for (int i = 0; i < n_Row; i ++)
		PCG_TMP_HELPER_IN[i] = i + 1.0;
	Preconditioner(n_Row, n_Col, PCG_TMP_HELPER_IN, PCG_TMP_HELPER_OUT);
	fprintf(O, "trueAns = Inv2' * Inv1' * inv(LDL_D) * Inv1 * Inv2 * tmpV;\n");
	fprintf(O, "testAns = zeros(n_Row, 1);\n");
	for (int i = 0; i < n_Row; i ++)
		fprintf(O, "testAns(%d) = %lf;\n", i + 1, PCG_TMP_HELPER_OUT[i]);
	fclose(O);
}
*/

/*
Solving (A (D + gamma I)^(-1) A^T + delta I) x = b
A: Sparse Matrix (CSR), n_Row * n_Col
D: Diagonal, n_Col * n_Col
b: Dense Vector, n_Row * 1
x: Dense Vector, n_Row * 1
M: Preconditioner

Preconditioned Conjugate Gradient (PCG)

Initialization: (k = 0)
	x_0 (initial guess)
	r_0 = (b - delta * x_0) - (A * ((D + gamma I)^(-1) * (A^T x_0)))
	z_0 = M^(-1) r_0
	p_0 = z_0

For k >= 0, while ||r_k|| / ||r_0|| > Epsilon
	q_k = delta * p_k + A * ((D + gamma I)^(-1) * (A^T * p_k))
	alpha_k = (z_k^T r_k) / (p_k^T q_k)
	x_{k + 1} = x_k + alpha_k * p_k
	r_{k + 1} = r_k - alpha_k * q_k
	z_{k + 1} = M^(-1) r_{k + 1}
	beta_k = (z_{k + 1}^T r_{k + 1}) / (z_k^T r_k)
	p_{k + 1} = z_{k + 1} + beta_k p_k
End

Here, M^(-1) is the preconditioner of (A (D + gamma I)^(-1) A^T + delta I)
*/

double PCG_TMP_ROW[MAX_ROWS];

void ConjugateGradient(double gamma, double delta, int n_Row, int n_Col, double* csrValAt, int* csrColIndAt, int* csrRowPtrAt, double* csrValA, int* csrColIndA, int* csrRowPtrA, 
					   double* d, double* b, double* x, double* tmp_col, double* r, double* z, double* p, double* q)
{
	for (int i = 0; i < n_Row; i ++)
		PCG_PERM[i] = i;
	RenewPartialCholesky(n_Row, n_Col, csrValA, csrColIndA, csrRowPtrA, d, gamma, delta); // Gamma and Delta need to be adjusted!
	//PCGDebugHelper(n_Row, n_Col, csrValA, csrColIndA, d, gamma, delta);
	
	char MKL_matdescra[6] = {0};
	MKL_matdescra[0] = 'G';
	MKL_matdescra[3] = 'C';

	// x_0: Initial Guess
	for (int i = 0; i < n_Row; i ++)
		x[i] = 0;
	
	// r_0 = b - (A (D + gamma I)^(-1) A^T + delta I) x
	//     = (b - delta x) - A * ((D + gamma I)^(-1) * (A^T * x))

	// Step 1: tmp_col = A^T * x
	//mkl_dcsrmv(&CHAR_T, &n_Row, &n_Col, &DOUBLE_ONE, MKL_matdescra, csrValA, csrColIndA, csrRowPtrA_BEG, csrRowPtrA_END, x, &DOUBLE_ZERO, tmp_col);
	CSRMV_T(n_Row, n_Col, DOUBLE_ONE, csrValA, csrColIndA, csrRowPtrA_BEG, csrRowPtrA_END, x, DOUBLE_ZERO, tmp_col);
	// Step 2: tmp_col = (D + gamma I)^(-1) * tmp_col
	for (int i = 0; i < n_Col; i ++)
		tmp_col[i] /= d[i] + gamma;
	// Step 3: r_0 = b - delta x
	// A row reorder => b reorder
	for (int i = 0; i < n_Row; i ++)
		r[i] = b[PCG_PERM[i]] - delta * x[i];
	// Step 4: r_0 = r_0 - A * tmp_col
	//mkl_dcsrmv(&CHAR_N, &n_Row, &n_Col, &DOUBLE_NEGONE, MKL_matdescra, csrValA, csrColIndA, csrRowPtrA_BEG, csrRowPtrA_END, tmp_col, &DOUBLE_ONE, r);
	CSRMV_N(n_Row, n_Col, DOUBLE_NEGONE, csrValA, csrColIndA, csrRowPtrA_BEG, csrRowPtrA_END, tmp_col, DOUBLE_ONE, r);

	// z_0 = M^(-1) r_0
	Preconditioner(n_Row, n_Col, r, z);

	// p_0 = z_0
	cblas_dcopy(n_Row, z, 1, p, 1); // y = x

	double terminate = cblas_dnrm2(n_Row, r, 1) * PREC_THRESHOLD;

	int loop = 0;
	while (1)
	{
		loop ++;
		double rnorm = cblas_dnrm2(n_Row, r, 1);
		printf("%d, %e\n", loop, rnorm);
		if (rnorm < terminate)
			break;

		// q_k = delta * p_k + A * ((D + gamma I)^(-1) * (A^T * p_k))
		// Step 1: tmp_col = A^T * p_k
		//mkl_dcsrmv(&CHAR_T, &n_Row, &n_Col, &DOUBLE_ONE, MKL_matdescra, csrValA, csrColIndA, csrRowPtrA_BEG, csrRowPtrA_END, p, &DOUBLE_ZERO, tmp_col);
		CSRMV_T(n_Row, n_Col, DOUBLE_ONE, csrValA, csrColIndA, csrRowPtrA_BEG, csrRowPtrA_END, p, DOUBLE_ZERO, tmp_col);
		// Step 2: tmp_col = (D + gamma I)^(-1) * tmp_col
		for (int i = 0; i < n_Col; i ++)
			tmp_col[i] /= d[i] + gamma;
		// Step 3: q_k = delta * p_k
		for (int i = 0; i < n_Row; i ++)
			q[i] = delta * p[i];
		// Step 4: q = q + A * tmp_col
		//mkl_dcsrmv(&CHAR_N, &n_Row, &n_Col, &DOUBLE_ONE, MKL_matdescra, csrValA, csrColIndA, csrRowPtrA_BEG, csrRowPtrA_END, tmp_col, &DOUBLE_ONE, q);
		CSRMV_N(n_Row, n_Col, DOUBLE_ONE, csrValA, csrColIndA, csrRowPtrA_BEG, csrRowPtrA_END, tmp_col, DOUBLE_ONE, q);

		// alpha_k = (z_k^T r_k) / (p_k^T q_k)
		double alpha_up = cblas_ddot(n_Row, z, 1, r, 1);
		double alpha_down = cblas_ddot(n_Row, p, 1, q, 1);
		double alpha = alpha_up / alpha_down;

		// x_{k + 1} = x_k + alpha_k * p_k
		cblas_daxpy(n_Row, alpha, p, 1, x, 1); // y = alpha * x + y

		// TEST
		//for (int i = 0; i < n_Row; i ++) PCG_TMP_ROW[i] = r[i];

		// r_{k + 1} = r_k - alpha_k * q_k
		cblas_daxpy(n_Row, -alpha, q, 1, r, 1); // y = alpha * x + y

		// z_{k + 1} = M^(-1) r_{k + 1}
		Preconditioner(n_Row, n_Col, r, z);

		// TEST
		//printf("conj = %e\n", cblas_ddot(n_Row, PCG_TMP_ROW, 1, z, 1));

		// beta_k = (z_{k + 1}^T r_{k + 1}) / (z_k^T r_k)
		double beta_up = cblas_ddot(n_Row, z, 1, r, 1);
		double beta = beta_up / alpha_up;
		
		// p_{k + 1} = z_{k + 1} + beta_k p_k
		// Step 1: p = beta * p
		cblas_dscal(n_Row, beta, p, 1);
		// Step 2: p = p + z
		cblas_daxpy(n_Row, 1, z, 1, p, 1); // y = alpha * x + y
	}

	// Reorder x
	for (int i = 0; i < n_Row; i ++)
		PCG_TMP_ROW[PCG_PERM[i]] = x[i];
	cblas_dcopy(n_Row, PCG_TMP_ROW, 1, x, 1);
}






// M = [1 0 0 0 0; 2 1 0 0 0; 3 5 1 0 0; 2 4 3 1 0; 3 3 0 2 1];
// A = M * diag([4 1 5 3 2]) * M';
void UnitTest()
{
	int n_Row = 5;
	int n_Col = 5;
	int csrRowPtrA[] = {0, 1, 3, 6, 10, 15};
	int csrColIndA[] = {1, 1, 2, 1, 2, 3, 1, 2, 3, 4, 1, 2, 3, 4, 5};
	double csrValA[] = {1, 2, 1, 3, 5, 1, 2, 4, 3, 1, 3, 3, 0, 2, 1};
	double d[] = {0.25 * 0.25, 1, 0.2 * 0.2, 1.0 / 9.0, 0.25};
	double gamma = 0;
	double delta = 0;
	RenewPartialCholesky(n_Row, n_Col, csrValA, csrColIndA, csrRowPtrA, d, gamma, delta);

	/*

M = [1 0 0 0 0; 2 1 0 0 0; 3 5 1 0 0; 2 4 3 1 0; 3 3 0 2 1];
A = M * diag([4 1 5 3 2]) * diag([4 1 5 3 2]) * M';
G = [4 2 3 1 5];
A(G, G)

>> L = [1.0000, 0, 0, 0, 0; 0.2166, 1, 0, 0, 0; 0.6083, 0, 1, 0, 0; 0.1019, 0, 0, 1, 0; 0.4013, 0, 0, 0, 1]

L =

    1.0000         0         0         0         0
    0.2166    1.0000         0         0         0
    0.6083         0    1.0000         0         0
    0.1019         0         0    1.0000         0
    0.4013         0         0         0    1.0000

inv(L) * A(G, G) * inv(L')

ans =

  314.0000   -0.0124   -0.0062    0.0034   -0.0082
   -0.0124   50.2739   59.6369   25.0701   71.7134
   -0.0062   59.6369   77.8185   28.5350   82.3567
    0.0034   25.0701   28.5350   12.7389   35.1592
   -0.0082   71.7134   82.3567   35.1592  142.4395

M = [1 0 0 0 0; 2 1 0 0 0; 3 5 1 0 0; 2 4 3 1 0; 3 3 0 2 1];
A = M * diag([4 1 5 3 2]) * diag([4 1 5 3 2]) * M';
G = [4 5 3 1 2];
pa = A(G, G)

l = chol(pa)'

partl = [l(:, 1) / l(1, 1), l(:, 2) / l(2, 2), zeros(5, 3)];
partl(3, 3) = 1;
partl(4, 4) = 1;
partl(5, 5) = 1;

inv(partl) * A(G, G) * inv(partl')

	*/
}