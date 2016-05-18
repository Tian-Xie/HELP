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
int HSD_status; // 0 - OK, 1 - Stalled, 2 - Iteration Limit Exceeded, 3 - Unknown Error

double HSD_D[MAX_COLS], HSD_D_u[MAX_COLS]; // D: X^{-1} * S, the diagonal matrix
double HSD_D_g; // D_g = kappa / tau + D[n_LB + i] * (x[n_LB + i] / tau)^2
double HSD_r_p[MAX_ROWS]; // r_p = b * tau - A * x
double HSD_r_d[MAX_COLS]; // r_d = A^T * y + z, add some other terms
double HSD_r_g; // r_g = kappa + c^T * x - b^T y
double HSD_mu; // mu = sum(x_j * z_j, xu_k * zu_k, tau * kappa) / (1 + n_LB + n_UB)
double HSD_rfval; // rfval = (1 + n_LB + n_UB) * min(x_j * z_j, xu_k * zu_k, tau * kappa) / sum(x_j * z_j, xu_k * zu_k, tau * kappa)

inline double DotProduct(int n, double* a, double* b)
{
	double Ret = 0;
	for (int i = 0; i < n; i ++)
		Ret += a[i] * b[i];
	return Ret;
}

void HSD_Calc_Newton_RHS(double* D, double* D_u, double& D_g, double* r_p, double* r_d, double& r_g, double& mu, double& rfval)
{
	// D: X^{-1} * S, the diagonal matrix
	// r_p = b * tau - A * x
	// r_d = A^T * y + z, add some other terms
	// r_g = kappa + c^T * x - b^T y
	// mu = sum(x_j * z_j, xu_k * zu_k, tau * kappa) / (1 + n_LB + n_UB)
	// rfval = (1 + n_LB + n_UB) * min(x_j * z_j, xu_k * zu_k, tau * kappa) / sum(x_j * z_j, xu_k * zu_k, tau * kappa)

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
			HSD_status = HSD_STATUS_STALLED;
			printf("HSD_Calc_Newton_RHS: Negative Iterate.\n");
			return;
		}
	}
	for (int i = 0; i < n_UB; i ++)
	{
		D_u[i] = HSD_zu[i] / HSD_xu[i];
		double tmp = HSD_zu[i] * HSD_xu[i];
		xz_min = min(xz_min, tmp);
		xz_sum += tmp;
		if (D_u[i] <= 0.0)
		{
			HSD_status = HSD_STATUS_STALLED;
			printf("HSD_Calc_Newton_RHS: Negative Iterate.\n");
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

	// r_d = A^T * y + z, add some other terms
	for (int i = 0; i < n_LB; i ++)
		r_d[i] = HSD_z[i] - V_Cost[i] * HSD_tau;
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
}

void HSD_GetInitPoint()
{
	// Sect 4.4: InitPoint (x, tau, y, z, kappa) = (e, 1, 0, e, 1)
	HSD_status = HSD_STATUS_OK;
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

	double d_tau = 0.0;
	double n_All = 1.0 + n_LB + n_UB;
	double mu = (DotProduct(n_LB, HSD_x, HSD_z) + DotProduct(n_UB, HSD_xu, HSD_zu) + HSD_tau * HSD_kappa) / n_All;
	double inf0 = HSD_tau / HSD_kappa;

	HSD_rho = 1.0;
	HSD_Calc_Newton_RHS(HSD_D, HSD_D_u, HSD_D_g, HSD_r_p, HSD_r_d, HSD_r_g, HSD_mu, HSD_rfval);
	if (HSD_status != HSD_STATUS_OK)
		return;

}

int HSD_Main()
{
	CHOLMOD_Construct();
	
	HSD_GetInitPoint();

	CHOLMOD_Destruct();
    return 0;
}
