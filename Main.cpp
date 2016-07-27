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

int main(int argc, char* argv[])
{
	CheckError(Prog_Init(), "Initialization Failed!");
	
	// Filename
	if (argc == 1)
	{
		// Input Filename
		//strcpy(Filename, "LINDEP.mps");
		//strcpy(Filename, "QAP15.SIF");
		//strcpy(Filename, "BRANDY.SIF");
		//strcpy(Filename, "sparse2000");
		strcpy(Filename, "afiro.mps");
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
	if (LP_Status != LP_STATUS_OK)
	{
		Report_Solution();
		return 0;
	}
	//MPS_PrintMatrix(); // Debug only

	Tm = GetTime();
	CheckError(Presolve_Main(), "Presolving Failed!");
	printf("Presolving: %.2lf s\n", GetTime() - Tm);
	if (LP_Status != LP_STATUS_OK)
	{
		Report_Solution();
		return 0;
	}

	Tm = GetTime();
	CheckError(HSD_Main(), "Homogeneous and Self-Dual Numerical Solving Failed!");
	printf("Homogeneous and Self-Dual Numerical Solving: %.2lf s\n", GetTime() - Tm);
	if (LP_Status != LP_STATUS_OK)
	{
		Report_Solution();
		return 0;
	}
	Report_Solution();

#ifdef _MSC_VER
	system("pause");
#endif
	return 0;
}
