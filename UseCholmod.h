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

#ifndef _USE_CHOLMOD_H

#define _USE_CHOLMOD_H

// Using CHOLMOD library, by Timothy A. Davis
#include "cholmod.h"

// Setting_CPU.cpp or Setting_GPU.cpp 
void CHOLMOD_Setting(cholmod_common* CHOL_Com);

#endif
