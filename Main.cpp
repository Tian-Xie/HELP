/*****************************************************************************
*  LPSolver, An experimental implementation of Homogeneous and Self-Dual     *
*  Algorithm for Linear Programming.                                         *
*  Author: Tian Xie (Research Center for Management Science and Information  *
*          Analytics, Shanghai University of Finance and Economics)          *
*  Credits: Fundamental implementation idea originated from COPL_LP.         *
*           (Xiong Zhang and Yinyu Ye)                                       *
*           See http://web.stanford.edu/~yyye/Col.html .                     *
******************************************************************************/

#include "LP.h"

int main(int argc, char* argv[])
{
	CheckError(Prog_Init(), "Initialization Failed!");
	
	// Filename
	if (argc == 1)
	{
		// Input Filename
		strcpy(Filename, "sparse5000");
		//strcpy(Filename, "afiro.mps");
		//strcpy(Filename, "example.mps");
	}
	else
		strcpy(Filename, argv[1]);
	
	double Tm;
	
	Tm = GetTime();
	CheckError(MPS_ReadFile(), "MPS ReadFile Failed!");
	printf("MPS_ReadFile: %.2lf s\n", GetTime() - Tm);
	
	Tm = GetTime();
	CheckError(CRUSH_Main(), "Crushing Failed!");
	printf("Crushing: %.2lf s\n", GetTime() - Tm);

	Tm = GetTime();
	CheckError(Presolve_Main(), "Presolving Failed!");
	printf("Presolving: %.2lf s\n", GetTime() - Tm);

	Tm = GetTime();
	CheckError(HSD_Main(), "Homogeneous and Self-Dual Numerical Solving Failed!");
	printf("Homogeneous and Self-Dual Numerical Solving: %.2lf s\n", GetTime() - Tm);

#ifdef _MSC_VER
	system("pause");
#endif
	return 0;
}
