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

#ifndef _PCG_H

#include "LP.h"

const double PREC_THRESHOLD = 1e-10;
const int PART_CHOL_COLS = 100;
const int MAX_PART_CHOL_COLS = PART_CHOL_COLS;

extern double PARTIAL_GR[MAX_ROWS][MAX_PART_CHOL_COLS];
extern double PARTIAL_CHOL_L[MAX_ROWS][MAX_PART_CHOL_COLS];
extern double PARTIAL_CHOL_D[MAX_ROWS];
extern int PCG_PERM[MAX_ROWS];

void UnitTest();

void Preconditioner(int n_Row, int n_Col, double* r, double* Ret);

void ConjugateGradient(double gamma, double delta, int n_Row, int n_Col, double* csrValAt, int* csrColIndAt, int* csrRowPtrAt, double* csrValA, int* csrColIndA, int* csrRowPtrA, 
					   double* d, double* b, double* x, double* tmp_col, double* r, double* z, double* p, double* q);

#endif