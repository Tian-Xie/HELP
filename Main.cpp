#include "LP.h"

int main(int argc, char* argv[])
{
	CheckError(Prog_Init(), "Initialization Failed!");
	
	// Filename
	if (argc == 1)
	{
		// Input Filename
		strcpy(Filename, "ran-p");
	}
	else
		strcpy(Filename, argv[1]);
	
	CheckError(MPS_ReadFile(), "MPS ReadFile Failed!");
	
	
	
	return 0;
}