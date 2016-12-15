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

// Rearranged as: (a) x[0 ~ (n_UB - 1)]: With LB and UB; 
//                (b) x[n_UB ~ (n_LB - 1)]: With LB only; 
//                (c) x[n_LB ~ (n_LB + n_FR - 1)]: Free.
double HSD_x[MAX_COLS], HSD_y[MAX_ROWS], HSD_z[MAX_COLS]; // For variables with LB
double HSD_xu[MAX_COLS], HSD_yu[MAX_ROWS], HSD_zu[MAX_COLS]; // For variables with UB
double HSD_xf[MAX_COLS]; // For free variables
double HSD_tau, HSD_kappa;

double HSD_rho;
double HSD_inf0; // inf_0 = tau_0 / kappa_0

// Problem Status Related
double HSD_sum_r_norm; // ||r_p||_2^2 + ||r_d||_2^2 + ||r_g||_2^2
double HSD_primal_obj, HSD_dual_obj; // P. Obj = c^T * x / tau, D. Obj = b^T * y / tau
double HSD_comp_gap; // Complementary Gap = sum(x_j * z_j, xu_k * zu_k, tau * kappa)
double HSD_relative_gap; // Relative Gap = |PObj - DObj| / (1 + |DObj|)
// Primal Infeasibility = ||r_p||_2 / (tau + ||x||_2)
// Dual Infeasibility = ||r_d||_2 / (tau + sqrt(sum_LB(z_j^2) + sum_UB(zu_k^2 + yu_k^2) + sum(y_j^2)))
double HSD_primal_infeas, HSD_dual_infeas;
double HSD_infe; // Infe = (tau / kappa) / (tau0 / kappa0)
int HSD_Status; // 0 - OK, 1 - Stalled, 2 - Iteration Limit Exceeded, 3 - Unknown Error

double HSD_D[MAX_COLS], HSD_D_u[MAX_COLS]; // D: X^{-1} * S, the diagonal matrix
double HSD_D_g; // D_g = kappa / tau + D[n_LB + i] * (x[n_LB + i] / tau)^2
double HSD_r_p[MAX_ROWS]; // r_p = b * tau - A * x
double HSD_r_d[MAX_COLS]; // r_d = - c * tau + A^T * y + z; note that r_d is actually -r_d in Section 4!
double HSD_r_g; // r_g = kappa + c^T * x - b^T y
double HSD_mu; // mu = sum(x_j * z_j, xu_k * zu_k, tau * kappa) / (1 + n_LB + n_UB)
double HSD_rfval; // rfval = (1 + n_LB + n_UB) * min(x_j * z_j, xu_k * zu_k, tau * kappa) / sum(x_j * z_j, xu_k * zu_k, tau * kappa)

// D: X^{-1} * Z, the diagonal matrix
// r_p = b * tau - A * x
// r_d = - c * tau + A^T * y + z; note that r_d is actually -r_d in Section 4!
// r_g = kappa + c^T * x - b^T y
// mu = sum(x_j * z_j, xu_k * zu_k, tau * kappa) / (1 + n_LB + n_UB)
// rfval = (1 + n_LB + n_UB) * min(x_j * z_j, xu_k * zu_k, tau * kappa) / sum(x_j * z_j, xu_k * zu_k, tau * kappa)
void HSD_Calc_Newton_Parameters(double* D, double* D_u, double& D_g, double* r_p, double* r_d, double& r_g, double& mu, double& rfval)
{
	double xz_min = HSD_tau * HSD_kappa; // min(x_j * z_j, xu_k * zu_k, tau * kappa)
	double xz_sum = xz_min; // sum(x_j * z_j, xu_k * zu_k, tau * kappa)
	for (int i = 0; i < n_LB; i ++)
	{
		D[i] = HSD_z[i] / HSD_x[i];
		double tmp = HSD_z[i] * HSD_x[i];
		xz_min = min(xz_min, tmp);
		xz_sum += tmp;
		if (D[i] <= 0.0)
		{
			HSD_Status = HSD_STATUS_STALLED;
			printf("HSD_Calc_Newton_Parameters: Negative Iterate.\n");
			return;
		}
	}
	for (int i = 0; i < n_UB; i ++)
	{
		D_u[i] = HSD_zu[i] / HSD_xu[i];
		D[i] += D_u[i];
		double tmp = HSD_zu[i] * HSD_xu[i];
		xz_min = min(xz_min, tmp);
		xz_sum += tmp;
		if (D_u[i] <= 0.0)
		{
			HSD_Status = HSD_STATUS_STALLED;
			printf("HSD_Calc_Newton_Parameters: Negative Iterate.\n");
			return;
		}
	}
	double Tmp = 0;
	for (int i = 0; i < n_FR; i ++)
	{
		D[n_LB + i] = HSD_rho * HSD_xf[i];
		Tmp += HSD_xf[i] * HSD_x[n_LB + i] * HSD_x[n_LB + i];
	}
	D_g = HSD_kappa / HSD_tau + HSD_rho * Tmp / (HSD_tau * HSD_tau); // D_g = kappa / tau + D[n_LB + i] * (x[n_LB + i] / tau)^2

	mu = xz_sum / (1.0 + n_LB + n_UB);
	rfval = (1.0 + n_LB + n_UB) * xz_min / xz_sum;
	
	// r_p = tau * b - A * x
	SetScaledVector(n_Row, HSD_tau, V_RHS, r_p);
	SetATimesVector(false, -1, HSD_x, r_p);

	// r_d = - c * tau + A^T * y + z; note that r_d is actually -r_d in Section 4!
	for (int i = 0; i < n_LB; i ++)
		r_d[i] = HSD_z[i] - V_Cost[i] * HSD_tau;
    for (int i = 0; i < n_UB; i ++)
		r_d[i] += HSD_yu[i];
	for (int i = 0; i < n_FR; i ++)
		r_d[n_LB + i] = -V_Cost[n_LB + i] * HSD_tau;
	for (int i = n_LB + n_FR; i < n_Col; i ++) // Anything Fixed ???
		r_d[i] = 0;
	SetATimesVector(true, 1, HSD_y, r_d);

	// r_g = kappa + c^T * x - b^T * y
	double cTx = DotProduct(n_LB + n_FR, V_Cost, HSD_x);
	double bTy = DotProduct(n_Row, V_RHS, HSD_y) + DotProduct(n_UB, V_UB, HSD_yu);
	r_g = HSD_kappa + cTx - bTy;

	// Other 
	// ||r_p||_2^2 := ||b * tau - A * x||_2^2 + sum((tau * u[k] - xu[k] - x[j])^2)
	double HSD_rp_norm2_sq = DotProduct(n_Row, r_p, r_p);
	for (int i = 0; i < n_UB; i ++)
		HSD_rp_norm2_sq += (HSD_tau * V_UB[i] - HSD_xu[i] - HSD_x[i]) * (HSD_tau * V_UB[i] - HSD_xu[i] - HSD_x[i]);
	// ||r_d||_2^2 := sum((r_d_j)^2) + sum((zu[k] + yu[k])^2)
	double HSD_rd_norm2_sq = DotProduct(n_LB + n_FR, r_d, r_d);
	for (int i = 0; i < n_UB; i ++)
		HSD_rd_norm2_sq += (HSD_zu[i] + HSD_yu[i]) * (HSD_zu[i] + HSD_yu[i]);
	// ||r_g||_2^2
	double HSD_rg_norm2_sq = r_g * r_g;
	// ||r_p||_2^2 + ||r_d||_2^2 + ||r_g||_2^2
	HSD_sum_r_norm = HSD_rp_norm2_sq + HSD_rd_norm2_sq + HSD_rg_norm2_sq;
	// P. Obj = c^T * x / tau
	HSD_primal_obj = V_Cost_Intercept + cTx / HSD_tau;
	// D. Obj = b^T * y / tau
	HSD_dual_obj = V_Cost_Intercept + bTy / HSD_tau;
	// Complementary Gap = sum(x_j * z_j, xu_k * zu_k, tau * kappa)
	HSD_comp_gap = xz_sum;
	// Relative Gap = |PObj - DObj| / (1 + |DObj|)
	HSD_relative_gap = fabs(HSD_primal_obj - HSD_dual_obj) / (1.0 + fabs(HSD_dual_obj));
	// Primal Infeasibility = ||r_p||_2 / (tau + ||x||_2)
	HSD_primal_infeas = sqrt(HSD_rp_norm2_sq) / (HSD_tau + sqrt(DotProduct(n_Col, HSD_x, HSD_x)));
	// Dual Infeasibility = ||r_d||_2 / (tau + sqrt(sum_LB(z_j^2) + sum_UB(zu_k^2 + yu_k^2) + sum(y_j^2)))
	HSD_dual_infeas = sqrt(HSD_rd_norm2_sq) / (HSD_tau + sqrt(DotProduct(n_LB, HSD_z, HSD_z) + DotProduct(n_UB, HSD_zu, HSD_zu) + DotProduct(n_UB, HSD_yu, HSD_yu) + DotProduct(n_Row, HSD_y, HSD_y)));
	// Infe = (tau / kappa) / (tau0 / kappa0)
	HSD_infe = HSD_tau / (HSD_kappa * HSD_inf0);

	/*
	for (int i = 0; i < n_Row; i ++)
		printf("r_p[%d] = %lf\n", i, HSD_r_p[Row_OldToNew[i]]);
	for (int i = 0; i < n_Col; i ++)
		printf("r_d[%d] = %lf\n", i, HSD_r_d[Col_OldToNew[i]]);
	return;
	*/
}

double HSD_SLE_RHS1[MAX_COLS], HSD_SLE_RHS2[MAX_ROWS];
// Call SolveLinearEquation() in Helper.cpp
// Solve [ D  A^T ][ p] = [s], which is equivalent to [ -D  A^T ][p] = [-s]
//       [ A      ][-q]   [t]                         [  A      ][q]   [ t]
void HSD_SolveLinearEquation(double* D, double* D_u, double* s, double* s_u, double* t, double* t_u, double* p, double* p_u, double* q, double* q_u) 
{
	// RHS_j = { s_j + du_k * tu_k - su_k, x_j has upper limit
	//         { s_j                     , otherwise
	for (int i = 0; i < n_LB + n_FR; i ++)
		HSD_SLE_RHS1[i] = s[i];
	for (int i = 0; i < n_UB; i ++)
		HSD_SLE_RHS1[i] += D_u[i] * t_u[i] - s_u[i];
	for (int i = 0; i < n_Row; i ++)
		HSD_SLE_RHS2[i] = t[i];

	/*
	FILE* tmp = fopen("tmp.txt", "w");
	fprintf(tmp, "n = %d;\n", n_LB);
	fprintf(tmp, "m = %d;\n", n_Row);
	fprintf(tmp, "d = zeros(n, 1);\n");
	fprintf(tmp, "rhs1 = zeros(n, 1);\n");
	fprintf(tmp, "rhs2 = zeros(m, 1);\n");
	fprintf(tmp, "A = zeros(m, n);\n");
	for (int i = 0; i < n_LB; i ++)		fprintf(tmp, "d(%d) = %e;\n", i + 1, D[i]);
	for (int i = 0; i < n_LB; i ++)		fprintf(tmp,"rhs1(%d) = %e;\n", i + 1, HSD_SLE_RHS1[i]);
	for (int i = 0; i < n_Row; i ++)	fprintf(tmp,"rhs2(%d) = %e;\n", i + 1, HSD_SLE_RHS2[i]);
	for (int i = 0; i < n_Row; i ++)
		for (int j = V_Matrix_Row_Head[i]; j != -1; j = V_Matrix_Row_Next[j])
			fprintf(tmp,"A(%d, %d) = %e;\n", i + 1, V_Matrix_Col[j] + 1, V_Matrix_Value[j]);
	fprintf(tmp, "rhs = [rhs1; rhs2];\n");
	fprintf(tmp, "lhs = [diag(d) A'; A zeros(m, m)];\n");
	fprintf(tmp, "ans = inv(lhs) * rhs\n");
	fclose(tmp);
	*/
	SolveLinearEquation(D, HSD_SLE_RHS1, HSD_SLE_RHS2, p, q);
	
	/*
	for (int i = 0; i < n_LB + n_FR; i ++)		printf("sol1[%d]\t%e\n", i, p[i]);
	for (int i = 0; i < n_Row; i ++)		printf("sol2[%d]\t%e\n", i, q[i]);
	*/

	for (int i = 0; i < n_Row; i ++)
		q[i] = -q[i];
	// Now [ -D  A^T ][p] = [-s]
	//     [  A      ][q]   [ t]

	// Similarly, solve [ D_u  I ][ p_u] = [s_u]
	//                  [  I     ][-q_u]   [t_u]
	// The solution is: { p_u = t_u
	//                  { q_u = D_u p_u - s_u
	// Some additional term, pu[i] = tu[i] - p[i], qu[i] = pu[i] * Du[i] - su[i]
	for (int i = 0; i < n_UB; i ++)
	{
		p_u[i] = t_u[i] - p[i];
		q_u[i] = D_u[i] * p_u[i] - s_u[i];
	}
}

// (-c; b)^T * (p; q) in Section 5
// Return -sum(c_j p_j) - rho / tau * sum_{free} (xf_k x_j p_j) + b^T q + u^T q_u
double HSD_MixDots(double* p, double* q, double* q_u)
{
	double Ret = - DotProduct(n_LB + n_FR, V_Cost, p);
	double rhoDtau = (HSD_rho / HSD_tau);
	for (int i = 0; i < n_FR; i ++)
		Ret -= rhoDtau * HSD_xf[i] * HSD_x[i + n_LB] * p[i + n_LB];
	Ret += DotProduct(n_Row, V_RHS, q) + DotProduct(n_UB, V_UB, q_u);
	return Ret;
}

// Solve [ -D  A^T ][-p] = [c]
//       [  A      ][-q]   [b]
// DOT = -(-c; b)^T * (p; q)
double HSD_LEBOT_RHS1[MAX_COLS], HSD_LEBOT_RHS1_u[MAX_COLS], HSD_LEBOT_RHS2[MAX_ROWS], HSD_LEBOT_RHS2_u[MAX_COLS];
void HSD_LinearEquation_Bottom(double* D, double* D_u, double& DOT, double* p, double* p_u, double* q, double* q_u)
{
	// RHS1_j = { c_j - rho / tau * xf_k * x_j,  x_j is free
	//          { c_j,                           otherwise
	for (int i = 0; i < n_Col; i ++)
		HSD_LEBOT_RHS1[i] = V_Cost[i];
	double rhoDtau = HSD_rho / HSD_tau;
	for (int i = 0; i < n_FR; i ++)
		HSD_LEBOT_RHS1[n_LB + i] -= rhoDtau * HSD_xf[i] * HSD_x[n_LB + i];
	// RHS2 = -b
	for (int i = 0; i < n_Row; i ++)
		HSD_LEBOT_RHS2[i] = -V_RHS[i];
	// The argumented variable, RHS1_u = 0, RHS2_u = -u
	for (int i = 0; i < n_UB; i ++)
	{
		HSD_LEBOT_RHS1_u[i] = 0;
		HSD_LEBOT_RHS2_u[i] = -V_UB[i];
	}

	// Solve [ -D  A^T ][-p] = [c]
	//       [  A      ][-q]   [b]
	HSD_SolveLinearEquation(D, D_u, HSD_LEBOT_RHS1, HSD_LEBOT_RHS1_u, HSD_LEBOT_RHS2, HSD_LEBOT_RHS2_u, p, p_u, q, q_u);
	/*
  int k;
  for (k = 0; k < n_Col; k ++)
	  printf("p[%d] = %lf\n", k, p[Col_OldToNew[k]]);
  for (k = 0; k < n_Row; k ++)
	  printf("q[%d] = %lf\n", k, q[k]);
  for (k = 0; k < n_UB; k ++)
	  printf("pu[%d] = %lf\n", k, p_u[k]);
  for (k = 0; k < n_UB; k ++)
	  printf("qu[%d] = %lf\n", k, q_u[k]);  
	  */
	// DOT = -(-c; b)^T * (p; q) 
	DOT = HSD_MixDots(p, q, q_u);
}

// (Affine Scaling) Newton Direction. 
// Solve the following linear system. 
// [   A     -b                 ] [  dx  ]   [hat_r_p ]   [        eta * r_p      ]
// [         -c   A^T   I       ] [ dtau ]   [hat_r_d ]   [        eta * r_d      ]
// [ -c^T         b^T        -1 ] [  dy  ] = [hat_r_g ] = [        eta * r_g      ]
// [   Z                X       ] [  dz  ]   [hat_r_xz]   [    -Xz + gamma mu e   ] <- may be changed
// [       kappa            tau ] [dkappa]   [hat_r_tk]   [ -tau kappa + gamma mu ] <- may be changed 
// eta: Input Parameter
// D: -X^{-1} Z, D_u: -X_u^{-1} Z_u, D_g: kappa / tau
// [r_p, r_d, r_g]: Input Parameter, note that r_d is actually -r_d in Section 4!
// DOT: -(-c; b)^T * (p; q) in Section 5
// p, p_u, q, q_u: Solved K * [-p; -q] = [c; b] in HSD_LinearEquation_Bottom
// hat_r_xz, hat_r_xuzu, hat_r_tk: Input Parameter
// [dx, dxu, dtau, dy, dyu, dz, dzu, dkappa]: Output Parameter
// r_p_norm, r_d_norm, r_g_norm: Updated ||r_p||_2^2, ||r_d||_2^2 and ||r_g||_2^2
double HSD_Dir_RHS1[MAX_COLS], HSD_Dir_RHS1_u[MAX_COLS], HSD_Dir_RHS2[MAX_ROWS], HSD_Dir_RHS2_u[MAX_COLS];
double HSD_Dir_p[MAX_COLS], HSD_Dir_p_u[MAX_COLS], HSD_Dir_q[MAX_ROWS], HSD_Dir_q_u[MAX_COLS];
double HSD_Dir_u[MAX_COLS], HSD_Dir_u_u[MAX_COLS], HSD_Dir_v[MAX_ROWS], HSD_Dir_v_u[MAX_COLS];
int HSD_Search_Direction(double eta, double* D, double* D_u, double D_g, double* r_p, double* r_d, double r_g, double DOT, 
						 double* p, double* p_u, double* q, double* q_u, double* hat_r_xz, double* hat_r_xuzu, double hat_r_tk, 
						 double* dx, double* dxu, double& dtau, double* dy, double* dyu, double* dz, double* dzu, double& dkappa, 
						 double& r_p_norm, double& r_d_norm, double& r_g_norm)
{
#ifdef DEBUG_TRACK
printf("In HSD_Search_Direction\n");
#endif
	// Use the relation that
	//     dz = X^{-1} (hat_r_xz - S dx)
	//     dkappa = (hat_r_tk - kappa dtau)
	// The reduced system
	//     [ -X^{-1}Z   A^T        -c      ] [ dx ]   [ hat_r_d - X^{-1} hat_r_xz ]
	//     [     A                 -b      ] [ dy ] = [ hat_r_p                   ]
	//     [   -c^T     b^T  (kappa / tau) ] [dtau]   [ hat_r_g + hat_r_tk / tau  ]
	// Remark that here [dz, dzu, dkappa] represents [hat_r_xz, hat_r_xuzu, hat_r_tk]

	// Solve [ -X^{-1}Z  A^T ][u] = [ hat_r_d - X^{-1} hat_r_xz ]
	//       [    A          ][v]   [ hat_r_p                   ]

	// First, Set RHS: RHS1 = eta * (-r_d);
	//                 RHS2 = eta * r_p;
	//                 RHS1_u = eta * (zu[j] + yu[j]);
	//                 RHS2_u = eta * (tau * u[j] - x[j] - xu[j]);
	SetScaledVector(n_Col, eta, r_d, HSD_Dir_RHS1);
	SetScaledVector(n_Row, eta, r_p, HSD_Dir_RHS2);
	for (int i = 0; i < n_UB; i ++)
	{
		HSD_Dir_RHS1_u[i] = eta * (HSD_zu[i] + HSD_yu[i]);
		HSD_Dir_RHS2_u[i] = eta * (HSD_tau * V_UB[i] - HSD_x[i] - HSD_xu[i]);
	}

	// To overcome numerical difficulty, set tilde_hat_r_xz[i] = { hat_r_xz[i] / x[i], if x[i] >= z[i]
	//                                                           { hat_r_xz[i] / z[i], if x[i] <  z[i]
	for (int i = 0; i < n_LB; i ++)
		hat_r_xz[i] /= max(HSD_x[i], HSD_z[i]);
	// Similarly tilde_hat_r_xuzu
	for (int i = 0; i < n_UB; i ++)
		hat_r_xuzu[i] /= max(HSD_xu[i], HSD_zu[i]);
	// Now hat_r_xz := tilde_hat_r_xz, hat_r_xuzu := tilde_hat_r_xuzu

	// RHS1 should be -(hat_r_d - X^{-1} hat_r_xz)
	// RHS2 should be hat_r_p
	
	// After a variable substitution: tilde_u[i] = { u[i],                      if x[i] >= z[i]
	//                                             { u[i] - hat_r_xz[i] / z[i], if x[i] <  z[i]
	// RHS1: If x[i] >= z[i], -(z[i] / x[i]) tilde_u[i] + Col(A, i)^T * v = hat_r_d[i] - hat_r_xz[i] / x[i]
	//       If x[i] <  z[i], -(z[i] / x[i]) tilde_u[i] + Col(A, i)^T * v = hat_r_d[i]
	// So RHS1[i] should be { -hat_r_d[i] + hat_r_xz / x[i], if x[i] >= z[i]
	//                      { -hat_r_d[i],                   if x[i] <  z[i]
	for (int i = 0; i < n_LB; i ++)
		if (HSD_x[i] >= HSD_z[i])
			HSD_Dir_RHS1[i] += hat_r_xz[i];
	for (int i = 0; i < n_UB; i ++)
		if (HSD_xu[i] >= HSD_zu[i])
			HSD_Dir_RHS1_u[i] += hat_r_xuzu[i];
	// RHS2: For Row i, Row(A, i) * tilde_u = hat_r_p[i] - sum(j: x[j] < z[j]) (A[i, j] * (hat_r_xz[i] / z[i]))
	// RHS2[i] should be hat_r_p[i] - sum(j: x[j] < z[j]) (A[i, j] * (hat_r_xz[i] / z[i]))
	for (int i = 0; i < n_LB; i ++)
		if (HSD_x[i] < HSD_z[i])
		{
			for (int j = V_Matrix_Col_Head[i]; j != -1; j = V_Matrix_Col_Next[j])
				HSD_Dir_RHS2[V_Matrix_Row[j]] -= V_Matrix_Value[j] * hat_r_xz[i];
		}
	for (int i = 0; i < n_UB; i ++)
	{
		if (HSD_x[i] < HSD_z[i])
			HSD_Dir_RHS2_u[i] -= hat_r_xz[i];
		if (HSD_xu[i] < HSD_zu[i])
			HSD_Dir_RHS2_u[i] -= hat_r_xuzu[i];
	}

	// Solve [ -D  A^T ][u] = [ hat_r_d - X^(-1) hat_r_xz ]
	//       [  A      ][v]   [ hat_r_p                   ]
	// which is equivalent to
	// [ -D  A^T ][tilde_u] = [-RHS_1]
	// [  A      ][   v   ]   [ RHS_2]
	
	// Calculating dtau = (hat_r_g + hat_r_tk / tau - (-c; b)^T (u; v)) / (kappa / tau + (-c; b)^T (p, q))
	// dtau = UP / DOWN
	double UP = eta * r_g + hat_r_tk / HSD_tau; // Currently partial UP
	// Adjusting tilde_u to u
	for (int i = 0; i < n_LB; i ++)
		if (HSD_x[i] < HSD_z[i])
			UP += V_Cost[i] * hat_r_xz[i];
	double DOWN = D_g - DOT; // Remark DOT = -(-c; b)^T (p; q)
	HSD_SolveLinearEquation(D, D_u, HSD_Dir_RHS1, HSD_Dir_RHS1_u, HSD_Dir_RHS2, HSD_Dir_RHS2_u, HSD_Dir_u, HSD_Dir_u_u, HSD_Dir_v, HSD_Dir_v_u);
	
	if (DOWN <= 1e-13)
		dtau = 0;
	else
		dtau = (UP - HSD_MixDots(HSD_Dir_u, HSD_Dir_v, HSD_Dir_v_u)) / max(D_g, DOWN);
	
	// Recall that (dx; dy) = (u; v) + (p; q) dtau

	// dx[i] = { (tilde_u[i] + p[i] * dtau),                      x[i] >= z[i]
	//         { (tilde_u[i] + p[i] * dtau) + hat_r_xz[i] / z[i], x[i] <  z[i]
	// dz = X^{-1} (hat_r_xz - Z * dx)
	// dz[i] = { hat_r_xz[i] / x[i] - (z[i] / x[i]) * dx[i],                      x[i] >= z[i]
	//         {                    - (z[i] / x[i]) * (tilde_u[i] + p[i] * dtau), x[i] <  z[i]
	// Note that [p, p_u] is actually [-p, -p_u] in Section 5!
	for (int i = 0; i < n_LB; i ++)
	{
		double Temp = HSD_Dir_u[i] - p[i] * dtau;
		if (HSD_x[i] >= HSD_z[i])
		{
			dx[i] = Temp;
			dz[i] = hat_r_xz[i] - HSD_z[i] * Temp / HSD_x[i];
		}
		else
		{
			dx[i] = hat_r_xz[i] + Temp;
			dz[i] = -HSD_z[i] * Temp / HSD_x[i];
		}
	}
	for (int i = 0; i < n_FR; i ++)
	{
		dx[i + n_LB] = HSD_Dir_u[i + n_LB] - p[i + n_LB] * dtau;
		dz[i + n_LB] = hat_r_xz[i + n_LB];
	}
	// Similarly, dxu and dzu
	for (int i = 0; i < n_UB; i ++)
	{
		double Temp = HSD_Dir_u_u[i] - p_u[i] * dtau;
		if (HSD_xu[i] >= HSD_zu[i])
		{
			dxu[i] = Temp;
			dzu[i] = hat_r_xuzu[i] - HSD_zu[i] * Temp / HSD_xu[i];
		}
		else
		{
			dxu[i] = hat_r_xuzu[i] + Temp;
			dzu[i] = -HSD_zu[i] * Temp / HSD_xu[i];
		}
	}
	// dy = v + q * dtau
	// Note that [q, q_u] is actually [-q, -q_u] in Section 5!
	for (int i = 0; i < n_Row; i ++)
		dy[i] = HSD_Dir_v[i] - q[i] * dtau;
	for (int i = 0; i < n_UB; i ++)
		dyu[i] = HSD_Dir_v_u[i] - q_u[i] * dtau;

	// dkappa = (hat_r_tk - kappa * dtau) / tau
	dkappa = (hat_r_tk - HSD_kappa * dtau) / HSD_tau;

	// Calc (eta + delta) [r_d, r_p, r_g]

	// RHS1 = eta * r_d - c * tau + dz + A^T dy;
	//      = eta * (- c * tau + A^T * y + z) - c * tau + dz + A^T dy;
	SetScaledVector(n_Col, eta, r_d, HSD_Dir_RHS1);
	for (int i = 0; i < n_LB; i ++)
		HSD_Dir_RHS1[i] += -V_Cost[i] * dtau + dz[i];
	double rhoDtau = HSD_rho / HSD_tau;
	for (int i = 0; i < n_FR; i ++)
		HSD_Dir_RHS1[i + n_LB] += -(V_Cost[i] - rhoDtau * HSD_xf[i] * HSD_x[i + n_LB]) * dtau - HSD_rho * HSD_xf[i] * dx[i + n_LB];
	SetATimesVector(true, 1, dy, HSD_Dir_RHS1);

	// RHS2 = eta * r_p + b * dtau - A * dx
	//      = eta * (b * tau - A * x) + b * dtau - A * dx
	for (int i = 0; i < n_Row; i ++)
		HSD_Dir_RHS2[i] = eta * HSD_r_p[i] + V_RHS[i] * dtau;
	SetATimesVector(false, -1, dx, HSD_Dir_RHS2);

	// RHS1_u = eta * (zu + yu) + dyu + dzu;
	// RHS2_u = eta * (tau * u[j] - x[j] - xu[j]) + u[j] * dtau - dxu[j] - dx[j];
	for (int i = 0; i < n_UB; i ++)
	{
		HSD_Dir_RHS1[i] += dyu[i];
		HSD_Dir_RHS1_u[i] = eta * (HSD_zu[i] + HSD_yu[i]) + dyu[i] + dzu[i];
		HSD_Dir_RHS2_u[i] = eta * (HSD_tau * V_UB[i] - HSD_x[i] - HSD_xu[i]) + V_UB[i] * dtau - dxu[i] - dx[i];
	}
	
	// Updated ||r_p||_2^2, ||r_d||_2^2 and ||r_g||_2^2
	r_p_norm = DotProduct(n_Row, HSD_Dir_RHS2, HSD_Dir_RHS2) + DotProduct(n_UB, HSD_Dir_RHS2_u, HSD_Dir_RHS2_u);
	r_d_norm = DotProduct(n_Col, HSD_Dir_RHS1, HSD_Dir_RHS1) + DotProduct(n_UB, HSD_Dir_RHS1_u, HSD_Dir_RHS1_u);
	// r_g = kappa + c^T * x - b^T * y
	double Tmp = 0;
	for (int i = 0; i < n_FR; i ++)
		Tmp += HSD_xf[i] * HSD_x[n_LB + i] * HSD_x[n_LB + i];
	double md = HSD_MixDots(dx, dy, dyu);
	r_g_norm = eta * r_g + dkappa - HSD_MixDots(dx, dy, dyu) - HSD_rho * Tmp * dtau / (HSD_tau * HSD_tau);
	r_g_norm = r_g_norm * r_g_norm;
	
#ifdef DEBUG_TRACK
printf("Out HSD_Search_Direction\n");
#endif
	// Some successful criteria ?????
	if (r_p_norm + r_d_norm <= 1e-10 || r_p_norm + r_d_norm + r_g_norm <= HSD_sum_r_norm * 0.1)
		return 1;
	return 0;
	
	//return 1;
}

// Stepsize in Section 4.1
// Partial stepsize: argmax_{0 <= alpha <= 1} {(x, xu) + alpha (dx, dxu) >= 0}
// Note that free variables are not considered!
double HSD_Stepsize(double* x, double* xu, double* dx, double* dxu)
{
	double Ret = 1.0;
	for (int i = 0; i < n_LB; i ++)
		if (dx[i] < -Variable_Tolerance && x[i] + Ret * dx[i] < 0.0)
			Ret = -x[i] / dx[i];
	for (int i = 0; i < n_UB; i ++)
		if (dxu[i] < -Variable_Tolerance && xu[i] + Ret * dxu[i] < 0.0)
			Ret = -xu[i] / dxu[i];
	return Ret;
}

// Calc Balanceness; alpha's are stepsizes
int HSD_CheckStepsize(double r_g, double eta, double* dx, double* dxu, double dtau, double* dz, double* dzu, double dkappa, 
					   double alpha_p, double alpha_tau, double alpha_d, double alpha_kappa, double beta)
{
	// min((x_j + alpha_p * dx_j) * (z_j + alpha_d * dz_j), (xu_k + alpha_p * dxu_k) * (zu_k + alpha_d * dzu_k), (tau + alpha_tau * dtau) * (kappa + alpha_kappa * dkappa))
	double xz_min = (HSD_tau + alpha_tau * dtau) * (HSD_kappa + alpha_kappa * dkappa);
	// sum((x_j + alpha_p * dx_j) * (z_j + alpha_d * dz_j), (xu_k + alpha_p * dxu_k) * (zu_k + alpha_d * dzu_k), (tau + alpha_tau * dtau) * (kappa + alpha_kappa * dkappa))
	double xz_sum = xz_min;
	for (int i = 0; i < n_LB; i ++)
	{
		double tmp = (HSD_x[i] + alpha_p * dx[i]) * (HSD_z[i] + alpha_d * dz[i]);
		xz_min = min(xz_min, tmp);
		xz_sum += tmp;
	}
	for (int i = 0; i < n_UB; i ++)
	{
		double tmp = (HSD_xu[i] + alpha_p * dxu[i]) * (HSD_zu[i] + alpha_d * dzu[i]);
		xz_min = min(xz_min, tmp);
		xz_sum += tmp;
	}
	// mu = sum(x_j * z_j, xu_k * zu_k, tau * kappa) / (1 + n_LB + n_UB)
	double mu = xz_sum / (1.0 + n_LB + n_UB);
	// LHS = (1 + n_LB + n_UB) * min(x_j * z_j, xu_k * zu_k, tau * kappa) / sum(x_j * z_j, xu_k * zu_k, tau * kappa)
	return (xz_min / mu >= beta);
}

double HSD_hat_r_xz[MAX_COLS], HSD_hat_r_xuzu[MAX_COLS], HSD_hat_r_tk;
double HSD_dx[MAX_COLS], HSD_dxu[MAX_COLS], HSD_dy[MAX_ROWS], HSD_dyu[MAX_COLS], HSD_dz[MAX_COLS], HSD_dzu[MAX_COLS];
double HSD_dtau, HSD_dkappa;
double HSD_dx_b[MAX_COLS], HSD_dxu_b[MAX_COLS], HSD_dy_b[MAX_ROWS], HSD_dyu_b[MAX_ROWS], HSD_dz_b[MAX_COLS], HSD_dzu_b[MAX_COLS];

void HSD_GetStepsize(double& p_size, double& d_size)
{
	p_size = HSD_Stepsize(HSD_x, HSD_xu, HSD_dx, HSD_dxu);
	d_size = HSD_Stepsize(HSD_z, HSD_zu, HSD_dz, HSD_dzu);
	if (HSD_dtau < -Variable_Tolerance)
	{
		p_size = min(p_size, -HSD_tau / HSD_dtau);
		d_size = min(d_size, -HSD_tau / HSD_dtau);
	}
	if (HSD_dkappa < -Variable_Tolerance)
	{
		p_size = min(p_size, -HSD_kappa / HSD_dkappa);
		d_size = min(d_size, -HSD_kappa / HSD_dkappa);
	}
}

void HSD_UpdateStep(double eta, double p_step, double d_step)
{
	double tau_step, kappa_step;
	int PGD = (p_step * HSD_dtau) > (d_step * HSD_dtau);
	if (PGD) // TODO: FOR WHAT REASON?
	{
		tau_step = d_step;
		kappa_step = p_step;
	}
	else
	{
		tau_step = p_step;
		kappa_step = d_step;
	}

	// Self-adjust step length
	double rstep = STEPSIZE_GAMMA;
	double threshold = max(min(sqrt(STEPSIZE_BETA), HSD_rfval * 0.01), STEPSIZE_BETA);
	while (! HSD_CheckStepsize(HSD_r_g, eta, HSD_dx, HSD_dxu, HSD_dtau, HSD_dz, HSD_dzu, HSD_dkappa, 
		rstep * p_step, rstep * tau_step, rstep * d_step, rstep * kappa_step, threshold))
		rstep *= 0.95;
	p_step *= rstep;
	tau_step *= rstep;
	d_step *= rstep;
	kappa_step *= rstep;

	double tau_p = HSD_tau + p_step * HSD_dtau;
	double tau_d = HSD_tau + d_step * HSD_dtau;
	double p_factor, d_factor;
	if (PGD)
	{
		p_factor = tau_d / tau_p;
		d_factor = 1;
	}
	else
	{
		p_factor = 1;
		d_factor = tau_p / tau_d;
	}
	for (int i = 0; i < n_LB; i ++)
	{
		HSD_x[i] = (HSD_x[i] + p_step * HSD_dx[i]) * p_factor;
		HSD_z[i] = (HSD_z[i] + d_step * HSD_dz[i]) * d_factor;
	}
	for (int i = 0; i < n_Row; i ++)
		HSD_y[i] = (HSD_y[i] + d_step * HSD_dy[i]) * d_factor;
	for (int i = 0; i < n_UB; i ++)
	{
		HSD_xu[i] = (HSD_xu[i] + p_step * HSD_dxu[i]) * p_factor;
		HSD_yu[i] = (HSD_yu[i] + d_step * HSD_dyu[i]) * d_factor;
		HSD_zu[i] = (HSD_zu[i] + d_step * HSD_dzu[i]) * d_factor;
	}
	for (int i = 0; i < n_FR; i ++)
		HSD_x[i + n_LB] = (HSD_x[i + n_LB] + p_step * HSD_dx[i + n_LB]) * p_factor;
	
	HSD_tau += tau_step * HSD_dtau;
	HSD_kappa = (HSD_kappa + kappa_step * HSD_dkappa) * p_factor * d_factor;
}

void HSD_GetInitPoint()
{
	// Section 4.4: InitPoint (x, tau, y, z, kappa) = (e, 1, 0, e, 1)
	HSD_Status = HSD_STATUS_OK;
	for (int i = 0; i < n_LB; i ++)
		HSD_x[i] = HSD_z[i] = 1.0;
	for (int i = 0; i < n_UB; i ++)
	{
		HSD_xu[i] = HSD_zu[i] = 1.0;
		HSD_yu[i] = 0.0;
	}
	for (int i = 0; i < n_Row; i ++)
		HSD_y[i] = 0.0;
	HSD_tau = HSD_kappa = 1.0;
	
	// For free variable, (x, z) = (0, 0), xf_i = 1 / (1 + 10 |x_i|)
	for (int i = 0; i < n_FR; i ++)
	{
		HSD_x[n_LB + i] = HSD_z[n_LB + i] = 0.0;
		HSD_xf[i] = 1.0 / (1.0 + 10.0 * fabs(HSD_x[n_LB + i]));
	}

	double n_All = 1.0 + n_LB + n_UB;
	double mu = (DotProduct(n_LB, HSD_x, HSD_z) + DotProduct(n_UB, HSD_xu, HSD_zu) + HSD_tau * HSD_kappa) / n_All;
	HSD_inf0 = HSD_tau / HSD_kappa;

	HSD_rho = 1.0;
	HSD_Calc_Newton_Parameters(HSD_D, HSD_D_u, HSD_D_g, HSD_r_p, HSD_r_d, HSD_r_g, HSD_mu, HSD_rfval);
	if (HSD_Status != HSD_STATUS_OK)
		return;
	// Remember to factorize ADA^T!
	RenewLinearEquation(HSD_D);
	
	double DOT;
	HSD_LinearEquation_Bottom(HSD_D, HSD_D_u, DOT, HSD_Dir_p, HSD_Dir_p_u, HSD_Dir_q, HSD_Dir_q_u);
	double r_p_norm, r_d_norm, r_g_norm;
	
	// Calculate Direction and Stepsize
	double phi = 1.0; // 1.0: Newton Step is correct
	double eta = 1.0;

	// Section 4.4, gamma = mu = 1
	// hat_r_xz = -X * s + mu * e
	// hat_r_tk = -tau * kappa + mu
	for (int i = 0; i < n_LB; i ++)
		HSD_hat_r_xz[i] = -HSD_x[i] * HSD_z[i] + mu;
	for (int i = 0; i < n_FR; i ++)
		HSD_hat_r_xz[i + n_LB] = 0;
	for (int i = 0; i < n_UB; i ++)
		HSD_hat_r_xuzu[i] = -HSD_xu[i] * HSD_zu[i] + mu;
	HSD_hat_r_tk = -HSD_tau * HSD_kappa + mu;

	int Newton_Success = HSD_Search_Direction(eta, HSD_D, HSD_D_u, HSD_D_g, HSD_r_p, HSD_r_d, HSD_r_g, DOT, 
		HSD_Dir_p, HSD_Dir_p_u, HSD_Dir_q, HSD_Dir_q_u, HSD_hat_r_xz, HSD_hat_r_xuzu, HSD_hat_r_tk, 
		HSD_dx, HSD_dxu, HSD_dtau, HSD_dy, HSD_dyu, HSD_dz, HSD_dzu, HSD_dkappa, 
		r_p_norm, r_d_norm, r_g_norm);

	double p_step = 0; // argmax_{0 <= p_step <= 1} {(x, xu, tau) + p_step (dx, dxu, dtau) >= 0}
	double d_step = 0; // argmax_{0 <= d_step <= 1} {(z, zu, kappa) + d_step (dz, dzu, dkappa) >= 0}
	if (Newton_Success)
		HSD_GetStepsize(p_step, d_step);
	else
	{
		printf("HSD_GetInitPoint: Warning: Inaccurate Newton Direction.\n");
		phi = 0.0;
	}

	// Section 4.4, an improvement of gamma = 10 and mu = 1
	// min(p_size, d_size): maximal possible step for this direction
	double gamma = 10.0;
	double alpha_p = max(0.0, 1.0 - min(p_step, d_step)) * gamma * mu; // (1 - alpha) * gamma * mu in (8.23)
	// Below is not consistent to Sect 4.4 but the idea is similar
	// dmu: ignore higher-order information
	double dmu = max((DotProduct(n_LB, HSD_dx, HSD_dz) + DotProduct(n_UB, HSD_dxu, HSD_dzu) + HSD_dtau * HSD_dkappa) / n_All, 0.0);
	for (int i = 0; i < n_LB; i ++)
		HSD_hat_r_xz[i] = -HSD_x[i] * HSD_z[i] + alpha_p - phi * (HSD_dx[i] * HSD_dz[i] - dmu);
	for (int i = 0; i < n_FR; i ++)
		HSD_hat_r_xz[i + n_LB] = 0;
	for (int i = 0; i < n_UB; i ++)
		HSD_hat_r_xuzu[i] = -HSD_xu[i] * HSD_zu[i] + alpha_p - phi * (HSD_dxu[i] * HSD_dzu[i] - dmu);
	HSD_hat_r_tk = -HSD_tau * HSD_kappa + alpha_p - phi * (HSD_dtau * HSD_dkappa - dmu);

	double HSD_dtau_b, HSD_dkappa_b;
	double r_p_norm_b, r_d_norm_b, r_g_norm_b;
	int IMP_Success = HSD_Search_Direction(eta, HSD_D, HSD_D_u, HSD_D_g, HSD_r_p, HSD_r_d, HSD_r_g, DOT, 
		HSD_Dir_p, HSD_Dir_p_u, HSD_Dir_q, HSD_Dir_q_u, HSD_hat_r_xz, HSD_hat_r_xuzu, HSD_hat_r_tk, 
		HSD_dx_b, HSD_dxu_b, HSD_dtau_b, HSD_dy_b, HSD_dyu_b, HSD_dz_b, HSD_dzu_b, HSD_dkappa_b, 
		r_p_norm_b, r_d_norm_b, r_g_norm_b);

	if (! IMP_Success)
	{
		printf("HSD_GetInitPoint: Warning: Inaccurate Improved Direction.\n");
		if (! Newton_Success) // Fail
		{
			printf("HSD_GetInitPoint: Warning: Failed to generate Newton-like starting point.\n");
			return;
		}
	}
	else
	{
		// Use improved direction and dispose the original Newton direction
		memcpy(HSD_dx, HSD_dx_b, sizeof(double) * n_Col);
		memcpy(HSD_dxu, HSD_dxu_b, sizeof(double) * n_UB);
		memcpy(HSD_dy, HSD_dy_b, sizeof(double) * n_Row);
		memcpy(HSD_dyu, HSD_dyu_b, sizeof(double) * n_UB);
		memcpy(HSD_dz, HSD_dz_b, sizeof(double) * n_Col);
		memcpy(HSD_dzu, HSD_dzu_b, sizeof(double) * n_UB);
		HSD_dtau = HSD_dtau_b;
		HSD_dkappa = HSD_dkappa_b;
		r_p_norm = r_p_norm_b;
		r_d_norm = r_d_norm_b;
		r_g_norm = r_g_norm_b;
		HSD_GetStepsize(p_step, d_step);
	}

	if (p_step < Variable_Tolerance && d_step < Variable_Tolerance)
	{
		printf("HSD_GetInitPoint: Warning: Newton-like starting point with 0 stepsize, stop here.\n");
		return;
	}
	HSD_UpdateStep(eta, p_step, d_step);
}

void HSD_Print()
{
	FILE* Out = fopen("HSD_LOG.txt", "w");
	for (int i = 0; i < n_LB + n_FR; i ++)
		fprintf(Out, "HSD_x[%d] = %lf\n", i, HSD_x[i]);
	for (int i = 0; i < n_LB; i ++)
		fprintf(Out, "HSD_z[%d] = %lf\n", i, HSD_z[i]);
	for (int i = 0; i < n_Row; i ++)
		fprintf(Out, "HSD_y[%d] = %lf\n", i, HSD_y[i]);
	for (int i = 0; i < n_UB; i ++)
		fprintf(Out, "HSD_xu[%d] = %lf\n", i, HSD_xu[i]);
	for (int i = 0; i < n_UB; i ++)
		fprintf(Out, "HSD_yu[%d] = %lf\n", i, HSD_yu[i]);
	for (int i = 0; i < n_UB; i ++)
		fprintf(Out, "HSD_zu[%d] = %lf\n", i, HSD_zu[i]);
	fclose(Out);
}

void HSD_Solving()
{
	HSD_GetInitPoint();
	double mu0 = (DotProduct(n_LB, HSD_x, HSD_z) + DotProduct(n_UB, HSD_xu, HSD_zu) + HSD_tau * HSD_kappa) / (1.0 + n_LB + n_UB);
	HSD_mu = mu0;
	
	printf("+------+------------------+------------------+----------+----------+----------+\n");
	printf("| Iter |      P. Obj.     |      D. Obj.     |  Rel.Gap | P.Infeas | D.Infeas |\n");
	printf("+------+------------------+------------------+----------+----------+----------+\n");
	int Iter;
	for (Iter = 0; Iter < Max_Iterations; Iter ++)
	{
#ifdef DEBUG_TRACK
printf("Into Iteration: Iter = %d\n", Iter);
#endif
		HSD_rho = max(HSD_mu / mu0, 1e-12);
		// For free variable, xf_i = 1 / (1 + 10 |x_i|)
		for (int i = 0; i < n_FR; i ++)
			HSD_xf[i] = 1.0 / (1.0 + 10.0 * fabs(HSD_x[n_LB + i]));
		
		HSD_Calc_Newton_Parameters(HSD_D, HSD_D_u, HSD_D_g, HSD_r_p, HSD_r_d, HSD_r_g, HSD_mu, HSD_rfval);
		printf("| %4d | %16.7le | %16.7le | %8.1le | %8.1le | %8.1le |\n", Iter, HSD_primal_obj, HSD_dual_obj, HSD_relative_gap, HSD_primal_infeas, HSD_dual_infeas);
		printf("+------+------------------+------------------+----------+----------+----------+\n");

		// Stopping Criteria
		if (HSD_Status != HSD_STATUS_OK || // Exception
			(HSD_mu <= Mu_Tolerance * mu0 && HSD_infe < Infeasibility_Tolerance) ||
			(HSD_relative_gap <= Gap_Tolerance && HSD_primal_infeas <= Primal_Infeasibility_Tolerance && HSD_dual_infeas <= Dual_Infeasibility_Tolerance))
			break;

		// Factorize ADA^T
		RenewLinearEquation(HSD_D);
		double DOT;
		HSD_LinearEquation_Bottom(HSD_D, HSD_D_u, DOT, HSD_Dir_p, HSD_Dir_p_u, HSD_Dir_q, HSD_Dir_q_u);
		double r_p_norm, r_d_norm, r_g_norm;
		
		// Calculate Direction and Stepsize
		double phi = 1.0; // 1.0: Newton Step is correct

		// Section 4.1, gamma = 0, eta = 1, Pure Newton (Affine Scaling) Direction
		double eta = 1.0;
		// hat_r_xz = -X * s
		// hat_r_tk = -tau * kappa
		for (int i = 0; i < n_LB; i ++)
			HSD_hat_r_xz[i] = -HSD_x[i] * HSD_z[i];
		for (int i = 0; i < n_FR; i ++)
			HSD_hat_r_xz[i + n_LB] = 0;
		for (int i = 0; i < n_UB; i ++)
			HSD_hat_r_xuzu[i] = -HSD_xu[i] * HSD_zu[i];
		HSD_hat_r_tk = -HSD_tau * HSD_kappa;

		int Newton_Success = HSD_Search_Direction(eta, HSD_D, HSD_D_u, HSD_D_g, HSD_r_p, HSD_r_d, HSD_r_g, DOT, 
			HSD_Dir_p, HSD_Dir_p_u, HSD_Dir_q, HSD_Dir_q_u, HSD_hat_r_xz, HSD_hat_r_xuzu, HSD_hat_r_tk, 
			HSD_dx, HSD_dxu, HSD_dtau, HSD_dy, HSD_dyu, HSD_dz, HSD_dzu, HSD_dkappa, 
			r_p_norm, r_d_norm, r_g_norm);

		double p_step = 0.0; // argmax_{0 <= p_step <= 1} {(x, xu, tau) + p_step (dx, dxu, dtau) >= 0}
		double d_step = 0.0; // argmax_{0 <= d_step <= 1} {(z, zu, kappa) + d_step (dz, dzu, dkappa) >= 0}
		double gamma;
		if (Newton_Success)
		{
			HSD_GetStepsize(p_step, d_step);
//printf("%e %e\n", p_step, d_step);
			double new_mu = (HSD_tau + p_step * HSD_dtau) * (HSD_kappa + d_step * HSD_dkappa);
			for (int i = 0; i < n_LB; i ++)
				new_mu += (HSD_x[i] + p_step * HSD_dx[i]) * (HSD_z[i] + d_step * HSD_dz[i]);
			for (int i = 0; i < n_UB; i ++)
				new_mu += (HSD_xu[i] + p_step * HSD_dxu[i]) * (HSD_zu[i] + d_step * HSD_dzu[i]);
			new_mu /= 1.0 + n_LB + n_UB;

			// Not consistent with Section 4.1
			gamma = new_mu / HSD_mu;
			gamma = min(gamma * gamma * gamma, 0.2 * min(gamma, 1.0));
			if (HSD_mu < 1e-3)
				gamma = max(gamma, HSD_mu);
			gamma = max(1e-6, gamma);
		}
		else
		{
			printf("HSD_Solving: Inaccurate Newton Direction, Stalled.\n");
			HSD_Status = HSD_STATUS_STALLED;
			break;
		}

		// The final direction
		eta = 1.0 - gamma;
		double alpha_p = gamma * HSD_mu;
		double dmu = max((DotProduct(n_LB, HSD_dx, HSD_dz) + DotProduct(n_UB, HSD_dxu, HSD_dzu) + HSD_dtau * HSD_dkappa) / (1.0 + n_LB + n_UB), 0.0);
		for (int i = 0; i < n_LB; i ++)
			HSD_hat_r_xz[i] = -HSD_x[i] * HSD_z[i] + alpha_p - phi * (HSD_dx[i] * HSD_dz[i] - dmu);
		for (int i = 0; i < n_FR; i ++)
			HSD_hat_r_xz[i + n_LB] = 0;
		for (int i = 0; i < n_UB; i ++)
			HSD_hat_r_xuzu[i] = -HSD_xu[i] * HSD_zu[i] + alpha_p - phi * (HSD_dxu[i] * HSD_dzu[i] - dmu);
		HSD_hat_r_tk = -HSD_tau * HSD_kappa + alpha_p - phi * (HSD_dtau * HSD_dkappa - dmu);

		double HSD_dtau_b, HSD_dkappa_b;
		double r_p_norm_b, r_d_norm_b, r_g_norm_b;
		int IMP_Success = HSD_Search_Direction(eta, HSD_D, HSD_D_u, HSD_D_g, HSD_r_p, HSD_r_d, HSD_r_g, DOT, 
			HSD_Dir_p, HSD_Dir_p_u, HSD_Dir_q, HSD_Dir_q_u, HSD_hat_r_xz, HSD_hat_r_xuzu, HSD_hat_r_tk, 
			HSD_dx_b, HSD_dxu_b, HSD_dtau_b, HSD_dy_b, HSD_dyu_b, HSD_dz_b, HSD_dzu_b, HSD_dkappa_b, 
			r_p_norm_b, r_d_norm_b, r_g_norm_b);
		if (IMP_Success)
		{
			// Use improved direction and dispose the original Newton direction
			memcpy(HSD_dx, HSD_dx_b, sizeof(double) * n_Col);
			memcpy(HSD_dxu, HSD_dxu_b, sizeof(double) * n_UB);
			memcpy(HSD_dy, HSD_dy_b, sizeof(double) * n_Row);
			memcpy(HSD_dyu, HSD_dyu_b, sizeof(double) * n_UB);
			memcpy(HSD_dz, HSD_dz_b, sizeof(double) * n_Col);
			memcpy(HSD_dzu, HSD_dzu_b, sizeof(double) * n_UB);
			HSD_dtau = HSD_dtau_b;
			HSD_dkappa = HSD_dkappa_b;
			r_p_norm = r_p_norm_b;
			r_d_norm = r_d_norm_b;
			r_g_norm = r_g_norm_b;
			HSD_GetStepsize(p_step, d_step);
		}
//printf("%e %e\n", p_step, d_step);
		if (p_step < Variable_Tolerance && d_step < Variable_Tolerance)
		{
			printf("HSD_Solving: Warning: Newton-like starting point with 0 stepsize, stop here.\n");
			HSD_Status = HSD_STATUS_UNKNOWN_ERROR;
			break;
		}
		HSD_UpdateStep(eta, p_step, d_step);
	}
	
	if (Iter >= Max_Iterations)
		HSD_Status = HSD_STATUS_ITER_LIMIT;
}

int HSD_Main()
{
	if (n_Row == 0 || n_Col == 0)
	{
		printf("No HSD Optimization is required.\n");
		return 0;
	}

	LinearEquation_Construct();
	
	HSD_Solving();

	LinearEquation_Destruct();
	return 0;
}
