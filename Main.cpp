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
#include "PCG.h"

#include "MKL_Util.h"
void MKLTest()
{
	// [1 0 2; 0 3 0; 4 0 5] * [3 5 7]' = [17 15 47]'
	char MKL_matdescra[6] = {0};
	MKL_matdescra[0] = 'G';
	MKL_matdescra[3] = 'C';
	int n_Col = 3;
	int n_Row = 3;
	double csrValA[] = {1, 2, 3, 4, 5};
	int csrColIndA[] = {0, 2, 1, 0, 2};
	int csrRowPtrA_BEG[] = {0, 2, 3};
	int csrRowPtrA_END[] = {2, 3, 5};
	double x[] = {3, 5, 7};
	double y[3];
	mkl_dcsrmv(&CHAR_N, &n_Col, &n_Row, &DOUBLE_ONE, MKL_matdescra, csrValA, csrColIndA, csrRowPtrA_BEG, csrRowPtrA_END, x, &DOUBLE_ZERO, y);
	printf("%lf\n%lf\n%lf\n", y[0], y[1], y[2]);
}

int main(int argc, char* argv[])
{
	//MKLTest();
	//UnitTest();
	//return 0;

	CheckError(Prog_Init(), "Initialization Failed!");
	
	char Buf[100];
	scanf("%s", Buf);

	// Filename
	if (argc == 1)
	{
		// Input Filename
		//strcpy(Filename, "LINDEP.mps");
		//strcpy(Filename, "LINDEPS.mps");
		//strcpy(Filename, "QAP15.SIF");
		//strcpy(Filename, "BRANDY.SIF");
		//strcpy(Filename, "sparse2000");
		//strcpy(Filename, "D:\\NETLIB\\afiro.mps");
		//strcpy(Filename, "D:\\NETLIB\\25FV47.mps");
		//strcpy(Filename, "example3.mps");
		//strcpy(Filename, "example.mps");
		//strcpy(Filename, "D:\\NETLIB\\MODSZK1.MPS");
		//strcpy(Filename, "D:\\NETLIB\\BOEING2.MPS");
		sprintf(Filename, "D:\\NETLIB\\%s.MPS", Buf);
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

	MPS_PrintMatrix(); // Debug only

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
