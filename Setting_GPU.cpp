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

#include "UseCholmod.h"

void CHOLMOD_Setting(cholmod_common* CHOL_Com)
{
	//CHOL_Com -> supernodal = CHOLMOD_SIMPLICIAL;
	//CHOL_Com -> supernodal = CHOLMOD_AUTO;
	CHOL_Com -> supernodal = CHOLMOD_SUPERNODAL;
	CHOL_Com -> useGPU = 1;
	CHOL_Com -> maxGpuMemFraction = 0.8;

	// Use AMD
	CHOL_Com -> nmethods = 1;
	CHOL_Com -> method[0] = CHOL_Com -> method[1];
}
