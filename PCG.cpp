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

// Calculate ret = M^(-1) * r
// Note that ``n_Row'' and ``n_Col'' are A's dimension size!
void Preconditioner(int n_Row, int n_Col, double* csrValAt, int* csrColIndAt, int* csrRowPtrAt, double* d, double* r, double* Ret)
{
	// Ret = r
	cblas_dcopy(n_Row, r, 1, Ret, 1); // y = x
	return;
}

/*
Solving (A D^(-1) A^T) x = b
A: Sparse Matrix (CSR), n_Row * n_Col
D: Diagonal, n_Col * n_Col
b: Dense Vector, n_Row * 1
x: Dense Vector, n_Row * 1

Preconditioned Conjugate Gradient (PCG)

Initialization: (k = 0)
	x_0 (initial guess)
	r_0 = b - (A * (D * (A^T x_0)))
	z_0 = M^(-1) r_0
	p_0 = z_0

For k >= 0, while ||r_k|| / ||r_0|| > Epsilon
	q_k = A * (D * (A^T * p_k))
	alpha_k = (z_k^T r_k) / (p_k^T q_k)
	x_{k + 1} = x_k + alpha_k * p_k
	r_{k + 1} = r_k - alpha_k * q_k
	z_{k + 1} = M^(-1) r_{k + 1}
	beta_k = (z_{k + 1}^T r_{k + 1}) / (z_k^T r_k)
	p_{k + 1} = r_{k + 1} + beta_k p_k
End

Here, M^(-1) is the preconditioner of (A D A^T)
*/

void ConjugateGradient(int n_Row, int n_Col, double* csrValAt, int* csrColIndAt, int* csrRowPtrAt, double* d, double* b, double* x, 
					   double* tmp_col, double* r, double* z, double* p, double* q)
{
	char MKL_matdescra[6] = {0};
	MKL_matdescra[0] = 'G';
	MKL_matdescra[3] = 'F';

	// x_0: Initial Guess
	for (int i = 0; i < n_Row; i ++)
		x[i] = 0;
	
	// r_0 = b - (A * (D * (A^T x)))
	// Step 1: tmp_col = A^T * x
	mkl_dcsrmv(&CHAR_N, &n_Col, &n_Row, &DOUBLE_ONE, MKL_matdescra, csrValAt, csrColIndAt, csrRowPtrAt, csrRowPtrAt + 1, x, &DOUBLE_ZERO, tmp_col);
	// Step 2: tmp_col = D * tmp_col
	for (int i = 0; i < n_Col; i ++)
		tmp_col[i] *= d[i];
	// Step 3: r_0 = b
	cblas_dcopy(n_Row, b, 1, r, 1); // y = x
	// Step 4: r_0 = r_0 - A * tmp_col
	mkl_dcsrmv(&CHAR_T, &n_Col, &n_Row, &DOUBLE_NEGONE, MKL_matdescra, csrValAt, csrColIndAt, csrRowPtrAt, csrRowPtrAt + 1, tmp_col, &DOUBLE_ONE, r);

	// z_0 = M^(-1) r_0
	Preconditioner(n_Row, n_Col, csrValAt, csrColIndAt, csrRowPtrAt, d, r, z);

	// p_0 = z_0
	cblas_dcopy(n_Row, z, 1, p, 1); // y = x

	mkl_dcsrmv(&CHAR_N, &n_Col, &n_Row, &DOUBLE_ONE, MKL_matdescra, csrValAt, csrColIndAt, csrRowPtrAt, csrRowPtrAt + 1, z, &DOUBLE_ZERO, tmp_col);

	double terminate = cblas_dnrm2(n_Row, r, 1) * prec_threshold;

	int loop = 0;
	while (1)
	{
		loop ++;
		double rnorm = cblas_dnrm2(n_Row, r, 1);
		printf("%d, %e\n", loop, rnorm);
		if (rnorm < terminate)
			break;

		// q_k = A * (D * (A^T * p_k))
		// Step 1: tmp_col = A^T * p_k
		mkl_dcsrmv(&CHAR_N, &n_Col, &n_Row, &DOUBLE_ONE, MKL_matdescra, csrValAt, csrColIndAt, csrRowPtrAt, csrRowPtrAt + 1, p, &DOUBLE_ZERO, tmp_col); // ??????????
		// Step 2: tmp_col = D * tmp_col
		for (int i = 0; i < n_Col; i ++)
			tmp_col[i] *= d[i];
		// Step 3: q = A * tmp_col
		mkl_dcsrmv(&CHAR_T, &n_Col, &n_Row, &DOUBLE_ONE, MKL_matdescra, csrValAt, csrColIndAt, csrRowPtrAt, csrRowPtrAt + 1, tmp_col, &DOUBLE_ZERO, q);

		// alpha_k = (z_k^T r_k) / (p_k^T q_k)
		double alpha_up = cblas_ddot(n_Row, z, 1, r, 1);
		double alpha_down = cblas_ddot(n_Row, p, 1, q, 1);
		double alpha = alpha_up / alpha_down;

		// x_{k + 1} = x_k + alpha_k * p_k
		cblas_daxpy(n_Row, alpha, p, 1, x, 1); // y = alpha * x + y

		// r_{k + 1} = r_k - alpha_k * q_k
		cblas_daxpy(n_Row, -alpha, q, 1, r, 1); // y = alpha * x + y

		// z_{k + 1} = M^(-1) r_{k + 1}
		Preconditioner(n_Row, n_Col, csrValAt, csrColIndAt, csrRowPtrAt, d, r, z);

		// beta_k = (z_{k + 1}^T r_{k + 1}) / (z_k^T r_k)
		double beta_up = cblas_ddot(n_Row, z, 1, r, 1);
		double beta = beta_up / alpha_up;
		
		// p_{k + 1} = r_{k + 1} + beta_k p_k
		// Step 1: p = beta * p
		cblas_dscal(n_Row, beta, p, 1);
		// Step 2: p = p + r
		cblas_daxpy(n_Row, 1, r, 1, p, 1); // y = alpha * x + y
	}
}