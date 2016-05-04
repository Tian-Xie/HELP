#include "LP.h"

void CheckError(int ExitID, char* ErrMsg)
{
	if (ExitID)
	{
		printf("ERROR: Msg = %s\n", ErrMsg);
		exit(ExitID);
	}
}

clock_t GetTime()
{
	return clock();
}