#include "LP.h"

int Temp[MAX_COLS];

void FinalizeModel()
{
	// Counting n_LB, n_UB, n_FR
	n_LB = n_UB = n_FR = 0;
	for (int i = 0; i < n_Col; i ++)
	{
		if (V_LB[i] == 0.0) // With LB
		{
			n_LB ++;
			if (V_UB[i] <= MaxFinite) // With LB and UB
				n_UB ++;
		}
		else
			n_FR ++;
	}
	// Rearranged as: (a) x[0 ~ (n_UB - 1)]: With LB and UB; 
	//                (b) x[n_UB ~ (n_LB - 1)]: With LB only; 
	//                (c) x[n_LB ~ (n_LB + n_FR - 1)]: Free.
	int Ptr1 = 0, Ptr2 = n_UB, Ptr3 = n_LB;
	for (int i = 0; i < n_Col; i ++)
		if (V_LB[i] == 0.0 && V_UB[i] <= MaxFinite) // With LB and UB
			TransOrder[Ptr1 ++] = i;
		else if (V_LB[i] == 0.0) // With LB only
			TransOrder[Ptr2 ++] = i;
		else
			TransOrder[Ptr3 ++] = i;
	// To call original index, just use x[RecoverOrder[i]] instead of x[i]. 
	for (int i = 0; i < n_Col; i ++)
		RecoverOrder[TransOrder[i]] = i;

	// Perform Reordering
	memset(Temp, 0, sizeof(int) * n_Col);
	for (int i = 0; i < n_Col; i ++)
	{
		int j = i;
		int Next = TransOrder[j];
		if (j == Next || Temp[i])
			continue;
		Temp[j] = 1;
		while (true)
		{
			swap(V_Cost[j], V_Cost[Next]);
			swap(V_LB[j], V_LB[Next]);
			swap(V_UB[j], V_UB[Next]);
			swap(V_Matrix_Head[j], V_Matrix_Head[Next]);
			swap(V_Crushing_Times[j], V_Crushing_Times[Next]);
			swap(V_Crushing_Add[j], V_Crushing_Add[Next]);
			Temp[Next] = 1;
			j = Next;
			Next = TransOrder[j];
			if (Temp[Next])
				break;
		}
	}
}

int Presolve_Main()
{
	FinalizeModel();
	return 0;
}