#include "LP.h"

int main(int argc, char* argv[])
{
	CheckError(Prog_Init(), "Initialization Failed!");
	
	// Filename
	if (argc == 1)
	{
		// Input Filename
		strcpy(Filename, "example.mps");
	}
	else
		strcpy(Filename, argv[1]);
	
	clock_t Tm;
	
	Tm = GetTime();
	CheckError(MPS_ReadFile(), "MPS ReadFile Failed!");
	printf("MPS_ReadFile: %d ms\n", GetTime() - Tm);
	
	Tm = GetTime();
	CheckError(CRUSH_Main(), "Crushing Failed!");
	printf("Crushing: %d ms\n", GetTime() - Tm);

	Tm = GetTime();
	CheckError(Presolve_Main(), "Presolving Failed!");
	printf("Presolving: %d ms\n", GetTime() - Tm);

	Tm = GetTime();
	CheckError(HSD_Main(), "Homogeneous and Self-Dual Numerical Solving Failed!");
	printf("Homogeneous and Self-Dual Numerical Solving: %d ms\n", GetTime() - Tm);

	system("pause");
	return 0;
}