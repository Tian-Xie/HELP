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
#include "cusolver_common.h"
#include <cuda_runtime.h>
#include "Helper_CUDA.h"
#include "cusparse.h"
#include "cusolverSp_LOWLEVEL_PREVIEW.h"

const double ONE = 1.0;
const double NEGONE = -1.0;

cusolverSpHandle_t cusolverSpHandle = NULL;
cusparseHandle_t cusparseHandle = NULL;
cusparseMatDescr_t descrM = NULL;
int nnzA; double csrValA[MAX_ELEMENTS]; int csrRowPtrA[MAX_ROWS], csrColIndA[MAX_ELEMENTS];
int Perm[MAX_ROWS];
int csrIndD[MAX_COLS]; double dInv[MAX_COLS];

double *d_csrValA; int *d_csrRowPtrA, *d_csrColIndA;
double *d_csrValD; int *d_csrIndD;
double *d_csrValAD; int *d_csrRowPtrAD, *d_csrColIndAD;
int max_nnzADAt, nnzADAt; double *d_csrValADAt; int *d_csrRowPtrADAt, *d_csrColIndADAt;

size_t pBufferSize;
void *pBuffer;

double tmp_row[MAX_ROWS], tmp_col[MAX_COLS];
double *dev_row, *dev_col;
double X1[MAX_COLS], X2[MAX_ROWS];
double *d_X1, *d_X2;
csrcholInfo_t csrcholInfo;

int LinearEquation_Construct()
{
	// Construct matrix A in cuSPARSE. 
	// After Presolve_Init(), we can assume A is sorted.
	nnzA = 0;
	for (int i = 0; i < n_Row; i ++)
	{
		csrRowPtrA[i] = nnzA;
		for (int p = V_Matrix_Row_Head[i]; p != -1; p = V_Matrix_Row_Next[p])
		{
			csrValA[nnzA] = V_Matrix_Value[p];
			csrColIndA[nnzA] = V_Matrix_Col[p];
			nnzA ++;
		}
	}
	csrRowPtrA[n_Row] = nnzA;
	// Construct indexs for D. Note that the csrRowPtrD and csrColIndD is identical. 
	for (int j = 0; j <= n_Col; j ++)
		csrIndD[j] = j;

	checkCudaErrors(cusparseCreate(&cusparseHandle));
	checkCudaErrors(cusparseCreateMatDescr(&descrM));
	checkCudaErrors(cusparseSetMatType(descrM, CUSPARSE_MATRIX_TYPE_GENERAL));
	checkCudaErrors(cusparseSetMatIndexBase(descrM, CUSPARSE_INDEX_BASE_ZERO));

	checkCudaErrors(cudaMalloc((void **) &d_csrRowPtrA, sizeof(int) * (n_Row + 1)));
	checkCudaErrors(cudaMalloc((void **) &d_csrColIndA, sizeof(int) * nnzA));
	checkCudaErrors(cudaMalloc((void **) &d_csrValA, sizeof(double) * nnzA));
	checkCudaErrors(cudaMalloc((void **) &d_csrIndD, sizeof(int) * (n_Col + 1)));
	checkCudaErrors(cudaMalloc((void **) &d_csrValD, sizeof(double) * n_Col));
	checkCudaErrors(cudaMalloc((void **) &d_csrRowPtrAD, sizeof(int) * (n_Row + 1)));
	checkCudaErrors(cudaMalloc((void **) &d_csrColIndAD, sizeof(int) * nnzA));
	checkCudaErrors(cudaMalloc((void **) &d_csrValAD, sizeof(double) * nnzA));
	checkCudaErrors(cudaMalloc((void **) &d_csrRowPtrADAt, sizeof(int) * (n_Row + 1)));
	checkCudaErrors(cudaMalloc((void **) &dev_row, sizeof(double) * n_Row));
	checkCudaErrors(cudaMalloc((void **) &dev_col, sizeof(double) * n_Col));
	checkCudaErrors(cudaMalloc((void **) &d_X1, sizeof(double) * n_Col));
	checkCudaErrors(cudaMalloc((void **) &d_X2, sizeof(double) * n_Row));

	checkCudaErrors(cusparseSetPointerMode(cusparseHandle, CUSPARSE_POINTER_MODE_HOST));

	checkCudaErrors(cudaMemcpy(d_csrRowPtrA, csrRowPtrA, sizeof(int) * (n_Row + 1), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_csrColIndA, csrColIndA, sizeof(int) * nnzA, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_csrValA, csrValA, sizeof(double) * nnzA, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_csrIndD, csrIndD, sizeof(int) * (n_Col + 1), cudaMemcpyHostToDevice));

	// Symbolic Calculate Maximized Elements of ADAt 
	int* nnzADAt_Ptr = &max_nnzADAt;
	checkCudaErrors(cusparseXcsrgemmNnz(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_TRANSPOSE, 
		n_Row, n_Row, n_Col, 
		descrM, nnzA, d_csrRowPtrA, d_csrColIndA,
		descrM, nnzA, d_csrRowPtrA, d_csrColIndA,
		descrM, d_csrRowPtrADAt, nnzADAt_Ptr));
	if (nnzADAt_Ptr)
		max_nnzADAt = *nnzADAt_Ptr;
	else
	{
		int tmp1, tmp2;
		checkCudaErrors(cudaMemcpy(&tmp1, d_csrRowPtrADAt + n_Row, sizeof(int), cudaMemcpyDeviceToHost));
		checkCudaErrors(cudaMemcpy(&tmp2, d_csrRowPtrADAt, sizeof(int), cudaMemcpyDeviceToHost));
		max_nnzADAt = tmp1 - tmp2;
	}
	checkCudaErrors(cudaMalloc((void **) &d_csrColIndADAt, sizeof(int) * max_nnzADAt));
	checkCudaErrors(cudaMalloc((void **) &d_csrValADAt, sizeof(double) * max_nnzADAt));
	// Test: D = I
	checkCudaErrors(cusparseDcsrgemm(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_TRANSPOSE, 
		n_Row, n_Row, n_Col, 
		descrM, nnzA, d_csrValA, d_csrRowPtrA, d_csrColIndA,
		descrM, nnzA, d_csrValA, d_csrRowPtrA, d_csrColIndA,
		descrM, d_csrValADAt, d_csrRowPtrADAt, d_csrColIndADAt));
	int tmp1, tmp2;
	checkCudaErrors(cudaMemcpy(&tmp1, d_csrRowPtrADAt + n_Row, sizeof(int), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(&tmp2, d_csrRowPtrADAt, sizeof(int), cudaMemcpyDeviceToHost));
	nnzADAt = tmp1 - tmp2;
	
	// Hereafter start cuSOLVER
	checkCudaErrors(cusolverSpCreate(&cusolverSpHandle));
	checkCudaErrors(cusolverSpCreateCsrcholInfo(&csrcholInfo));

	// Symbolic Analysis: Approx Minimum Degree Permutation, only CPU version
	int* csrRowPtrADAt = (int*) malloc(sizeof(int) * (n_Row + 1));
	int* csrColIndADAt = (int*) malloc(sizeof(int) * nnzADAt);
	cudaMemcpy(csrRowPtrADAt, d_csrRowPtrADAt, sizeof(int) * (n_Row + 1), cudaMemcpyDeviceToHost);
	cudaMemcpy(csrColIndADAt, d_csrColIndADAt, sizeof(int) * nnzADAt, cudaMemcpyDeviceToHost);
	// Always, DO NOT FORGET Perm!
#ifdef PRINT_TIME
	double Tm = GetTime();
#endif
#ifdef DEBUG_TRACK
	printf("Before Symbolic AMD\n");
#endif
	checkCudaErrors(cusolverSpXcsrsymamdHost(cusolverSpHandle, n_Row, nnzADAt, descrM, csrRowPtrADAt, csrColIndADAt, Perm));
#ifdef DEBUG_TRACK
	printf("After Symbolic AMD\n");
#endif
#ifdef PRINT_TIME
	printf("AMD: %.2lf s\n", GetTime() - Tm);
#endif
	free(csrRowPtrADAt);
	free(csrColIndADAt);
	for (int i = 0; i < n_Row; i ++) Perm[i] = i;
	// Reorder A rows!
	nnzA = 0;
	for (int i = 0; i < n_Row; i ++)
	{
		csrRowPtrA[i] = nnzA;
		for (int p = V_Matrix_Row_Head[Perm[i]]; p != -1; p = V_Matrix_Row_Next[p])
		{
			csrValA[nnzA] = V_Matrix_Value[p];
			csrColIndA[nnzA] = V_Matrix_Col[p];
			nnzA ++;
		}
	}
	csrRowPtrA[n_Row] = nnzA;
	// Then copy new A into GPU device
	checkCudaErrors(cudaMemcpy(d_csrRowPtrA, csrRowPtrA, sizeof(int) * (n_Row + 1), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_csrRowPtrAD, d_csrRowPtrA, sizeof(int) * (n_Row + 1), cudaMemcpyDeviceToDevice));
	checkCudaErrors(cudaMemcpy(d_csrColIndA, csrColIndA, sizeof(int) * nnzA, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_csrValA, csrValA, sizeof(double) * nnzA, cudaMemcpyHostToDevice));

	pBufferSize = 0;
	pBuffer = NULL;
	return 0;
}

int LinearEquation_Destruct()
{
	if (pBuffer) { pBufferSize = 0; checkCudaErrors(cudaFree(pBuffer)); pBuffer = NULL; }
	if (d_csrValA) { checkCudaErrors(cudaFree(d_csrValA)); d_csrValA = NULL; }
	if (d_csrRowPtrA) { checkCudaErrors(cudaFree(d_csrRowPtrA)); d_csrRowPtrA = NULL; }
	if (d_csrColIndA) { checkCudaErrors(cudaFree(d_csrColIndA)); d_csrColIndA = NULL; }
	if (d_csrValD) { checkCudaErrors(cudaFree(d_csrValD)); d_csrValD = NULL; }
	if (d_csrIndD) { checkCudaErrors(cudaFree(d_csrIndD)); d_csrIndD = NULL; }
	if (d_csrValAD) { checkCudaErrors(cudaFree(d_csrValAD)); d_csrValAD = NULL; }
	if (d_csrRowPtrAD) { checkCudaErrors(cudaFree(d_csrRowPtrAD)); d_csrRowPtrAD = NULL; }
	if (d_csrColIndAD) { checkCudaErrors(cudaFree(d_csrColIndAD)); d_csrColIndAD = NULL; }
	if (d_csrValADAt) { checkCudaErrors(cudaFree(d_csrValADAt)); d_csrValADAt = NULL; }
	if (d_csrRowPtrADAt) { checkCudaErrors(cudaFree(d_csrRowPtrADAt)); d_csrRowPtrADAt = NULL; }
	if (d_csrColIndADAt) { checkCudaErrors(cudaFree(d_csrColIndADAt)); d_csrColIndADAt = NULL; }
	if (dev_row) { checkCudaErrors(cudaFree(dev_row)); dev_row = NULL; }
	if (dev_col) { checkCudaErrors(cudaFree(dev_col)); dev_col = NULL; }
	if (d_X1) { checkCudaErrors(cudaFree(d_X1)); d_X1 = NULL; }
	if (d_X2) { checkCudaErrors(cudaFree(d_X2)); d_X2 = NULL; }
	if (csrcholInfo) { checkCudaErrors(cusolverSpDestroyCsrcholInfo(csrcholInfo)); csrcholInfo = NULL; }
	if (cusolverSpHandle) { checkCudaErrors(cusolverSpDestroy(cusolverSpHandle)); cusolverSpHandle = NULL; }
	if (cusparseHandle) { checkCudaErrors(cusparseDestroy(cusparseHandle)); cusparseHandle = NULL; }
	cudaDeviceReset();
	return 0;
}

// Renew ADA^T and factorize
void RenewLinearEquation(double* d)
{
	// Copy d to device
	for (int i = 0; i < n_Col; i ++)
		dInv[i] = 1.0 / d[i];
#ifdef PRINT_DEBUG
printf("arrd = zeros(%d, 1);\n", n_Col); for (int i = 0; i < n_Col; i ++) printf("arrd(%d) = %.10lf;\n", i + 1, dInv[i]); printf("D = diag(arrd);\n");
#endif
	checkCudaErrors(cudaMemcpy(d_csrValD, dInv, sizeof(double) * n_Col, cudaMemcpyHostToDevice));
	// Calc ADA^T
	// Step 1: AD. nnzAD should equal to nnzA. 
	checkCudaErrors(cusparseDcsrgemm(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_NON_TRANSPOSE, 
		n_Row, n_Col, n_Col, 
		descrM, nnzA, d_csrValA, d_csrRowPtrA, d_csrColIndA, 
		descrM, n_Col, d_csrValD, d_csrIndD, d_csrIndD, 
		descrM, d_csrValAD, d_csrRowPtrAD, d_csrColIndAD));
	// Step 2: ADAt. 
	checkCudaErrors(cusparseDcsrgemm(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, CUSPARSE_OPERATION_TRANSPOSE, 
		n_Row, n_Row, n_Col, 
		descrM, nnzA, d_csrValAD, d_csrRowPtrAD, d_csrColIndAD,
		descrM, nnzA, d_csrValA, d_csrRowPtrA, d_csrColIndA,
		descrM, d_csrValADAt, d_csrRowPtrADAt, d_csrColIndADAt));
	int tmp1, tmp2;
	checkCudaErrors(cudaMemcpy(&tmp1, d_csrRowPtrADAt + n_Row, sizeof(int), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(&tmp2, d_csrRowPtrADAt, sizeof(int), cudaMemcpyDeviceToHost));
	nnzADAt = tmp1 - tmp2;

	checkCudaErrors(cusolverSpXcsrcholAnalysis(cusolverSpHandle, n_Row, nnzADAt, descrM, d_csrRowPtrADAt, d_csrColIndADAt, csrcholInfo));

	size_t internalDataInBytes, workspaceInBytes;
	cusolverSpDcsrcholBufferInfo(cusolverSpHandle, n_Row, nnzADAt, descrM, d_csrValADAt, d_csrRowPtrADAt, d_csrColIndADAt, csrcholInfo, &internalDataInBytes, &workspaceInBytes);
	if (workspaceInBytes > pBufferSize)
	{
		if (pBuffer) { checkCudaErrors(cudaFree(pBuffer)); pBuffer = NULL; }
		pBufferSize = workspaceInBytes;
		checkCudaErrors(cudaMalloc((void **) &pBuffer, pBufferSize));
	}
#ifdef PRINT_TIME
	double Tm = GetTime();
#endif
#ifdef DEBUG_TRACK
	printf("%%Before Factorize\n");
#endif
	cusolverSpDcsrcholFactor(cusolverSpHandle, n_Row, nnzADAt, descrM, d_csrValADAt, d_csrRowPtrADAt, d_csrColIndADAt, csrcholInfo, pBuffer);
#ifdef DEBUG_TRACK
	printf("%%After Factorize\n");
#endif
#ifdef PRINT_TIME
	printf("%%Factorize: %.2lf s\n", GetTime() - Tm);
#endif
}

// Solve [ D  A^T ][x_1] == [b_1]
//       [ A   0  ][x_2]    [b_2]
// A: m * n; D: Diagonal Matrix, n * n; x_1: n * 1; x_2: m * 1; b_1: n * 1; b_2: m * 1
// Solution: (A D^(-1) A^T) x_2 = (A D^(-1) b_1 - b_2)
//           x_1 = D^(-1) (b_1 - A^T x_2)
// Need to call RenewLinearEquation(d) if ADA^T is not factorized yet!
int SolveLinearEquation(double* d, double* b_1, double* b_2, double* x_1, double* x_2)
{
#ifdef PRINT_DEBUG
/*
	printf("b_1 = zeros(%d, 1);\n", n_Col);
	for (int i = 0; i < n_Col; i ++)
		printf("b_1(%d) = %.6lf;\n", i + 1, b_1[i]);
	printf("b_2 = zeros(%d, 1);\n", n_Row);
	for (int i = 0; i < n_Row; i ++)
		printf("b_2(%d) = %.6lf;\n", i + 1, b_2[i]);
*/
#endif
	for (int i = 0; i < n_Row; i ++)
		tmp_row[i] = -b_2[Perm[i]];
	// D^(-1) b_1
	for (int i = 0; i < n_Col; i ++)
		tmp_col[i] = b_1[i] / d[i];
	checkCudaErrors(cudaMemcpy(dev_row, tmp_row, sizeof(double) * n_Row, cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_col, tmp_col, sizeof(double) * n_Col, cudaMemcpyHostToDevice));
	// A * (D^(-1) b_1) - b_2
	checkCudaErrors(cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, n_Row, n_Col, nnzA, &ONE, descrM, d_csrValA, d_csrRowPtrA, d_csrColIndA, dev_col, &ONE, dev_row));

	// Solve (A D^(-1) A^T) x_2 = (A D^(-1) b_1 - b_2)
#ifdef PRINT_DEBUG
checkCudaErrors(cudaMemcpy(tmp_row, dev_row, sizeof(double) * n_Row, cudaMemcpyDeviceToHost));
printf("RHS = zeros(%d, 1);\n", n_Row); for (int i = 0; i < n_Row; i ++) printf("RHS(%d) = %.10lf;\n", i + 1, tmp_row[i]);
printf("inv(A * D * A') * RHS\n");
#endif

#ifdef PRINT_TIME
	double Tm = GetTime();
#endif
#ifdef DEBUG_TRACK
printf("%%Before Solve\n");
#endif
	cusolverSpDcsrcholSolve(cusolverSpHandle, n_Row, dev_row, d_X2, csrcholInfo, pBuffer);
#ifdef DEBUG_TRACK
printf("%%After Solve\n");
#endif
#ifdef PRINT_TIME
	printf("%%Solve: %.2lf s\n", GetTime() - Tm);
#endif
	checkCudaErrors(cudaMemcpy(X2, d_X2, sizeof(double) * n_Row, cudaMemcpyDeviceToHost));
	for (int i = 0; i < n_Row; i ++)
		x_2[Perm[i]] = X2[i];
#ifdef PRINT_DEBUG
	for (int i = 0; i < n_Row; i ++)
		printf("x_2(%d) = %lf\n", i + 1, x_2[i]);
#endif
	// b_1 - A^T x_2
	checkCudaErrors(cudaMemcpy(dev_col, b_1, sizeof(double) * n_Col, cudaMemcpyHostToDevice));
	checkCudaErrors(cusparseDcsrmv(cusparseHandle, CUSPARSE_OPERATION_TRANSPOSE, n_Row, n_Col, nnzA, &NEGONE, descrM, d_csrValA, d_csrRowPtrA, d_csrColIndA, d_X2, &ONE, dev_col));
	checkCudaErrors(cudaMemcpy(tmp_col, dev_col, sizeof(double) * n_Col, cudaMemcpyDeviceToHost));
	// Let x_1 = D^(-1) (b_1 - A^T x_2)
	for (int i = 0; i < n_Col; i ++)
		x_1[i] = tmp_col[i] / d[i];
#ifdef PRINT_DEBUG
	for (int i = 0; i < n_Col; i ++)
		printf("x_1(%d) = %lf\n", i + 1, x_1[i]);
#endif
	return 0;
}
