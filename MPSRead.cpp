#include "LP.h"

const int CARD_POS[] = {4, 14, 24, 39, 49};
FILE* Fin;
const int MAX_LINE_LENGTH = 256;
char Buf[MAX_LINE_LENGTH];
unsigned long long Obj_Row;
map <unsigned long long, int> Map_Row, Map_Col;

void MPS_CopyName(char* Dest, char* Src)
{
	int Lenm1 = min(strlen(S), 8) - 1;
	while (Lenm1 >= 0 && (S[Lenm1] == ' ' || S[Lenm1] == '\r' || S[Lenm1] == '\n'))
		Lenm1 --;
	Dest[Lenm1 + 1] = 0;
	for (; Lenm1 >= 0; Lenm1 --)
		Dest[Lenm1] = Src[Lenm1];
}

unsigned long long MPS_HashCode(char* S)
{
	int Lenm1 = min(strlen(S), 8) - 1;
	unsigned long long Ret = 0;
	while (Lenm1 >= 0 && (S[Lenm1] == ' ' || S[Lenm1] == '\r' || S[Lenm1] == '\n'))
		Lenm1 --;
	for (; Lenm1 >= 0; Lenm1 --)
		Ret = (Ret << 8) | ((unsigned char) S[Lenm1]);
	return Ret;
}

void MPS_ReadLine()
{
	while (1)
	{
		if (feof(FIN))
			CheckError(1, "MPS_ReadLine: EOF before ENDATA!");
		if (fgets(Buf, MAX_LINE_LENGTH, Fin) == NULL)
			CheckError(1, "MPS_ReadLine: fgets ERROR!");
		if (Buf[0] != '*' && Buf[0] != '\r' && Buf[0] != '\n')
			break;
	}
}

void MPS_NAME()
{
	if (strncmp(Buf, "NAME", 4))
		CheckError(1, "MPS_ReadFile: Expected NAME!");
	MPS_CopyName(Problem_Name, Buf + 14);
	printf("    NAME SECTION\n");
	printf("        NAME = %s\n", Problem_Name);
	MPS_ReadLine(); // For MPS_ROWS
}

void MPS_ROWS()
{
	if (strncmp(Buf, "ROWS", 4))
		CheckError(1, "MPS_ReadFile: Expected ROWS!");
	printf("    ROWS SECTION\n");
	while (1)
	{
		MPS_ReadLine();
		if (Buf[0] != ' ') // End of ROWS
			break;
		char Type = toupper(Buf + 1);
		if (Type != 'N' && Type != 'E' && Type != 'G' && Type != 'L')
		{
			printf("        Warning: Unknown Type: \"%c\", Ignored.\n", Type);
			continue;
		}
		unsigned long long Hash = MPS_HashCode(Buf + CARD_POS[0]);
		if (Map_Row.find(Hash) != Map_Row.end())
		{
			printf("        Warning: Duplicate Row Name, HashCode = %llu, Ignored.\n", Hash);
			continue;
		}
		if (Type == 'N')
		{
			if (Obj_Row)
			{
				printf("        Warning: Duplicate Objective Row, Ignored.\n");
				continue;
			}
			Obj_Row = Hash; // Store Objective Row's Hash Value
			Map_Row[Hash] = -1; // -1: Objective
			continue;
		}
		Map[Hash] = n_Row;
		Row_Type.push_back(Type);
		V_RHS.push_back(0.0);
		n_Row ++;
	}
	printf("        %d Rows Read\n", n_Row);
	if (! Obj_Row)
		CheckError(1, "MPS_ReadFile: No Objective Row!");
}

void MPS_COLUMNS()
{
	if (strncmp(Buf, "COLUMNS", 7))
		CheckError(1, "MPS_ReadFile: Expected COLUMNS!");
	printf("    COLUMNS SECTION\n");
	while (1)
	{
		MPS_ReadLine();
		if (Buf[0] != ' ') // End of COLUMNS
			break;
		unsigned long long ColHash = MPS_HashCode(Buf + CARD_POS[0]);
		int ColID;
		if (Map_Col.find(ColHash) == Map_Col.end())
		{
			Map_Col[ColHash] = n_Col;
			ColID = n_Col;
			V_Cost.push_back(0.0);
			V_Matrix.push_back(map <int, double> ());
			n_Col ++;
		}
		else
			ColID = Map_Col[ColHash];

		int CurCardLen = strlen(Buf);
		for (int ip = 1; ip <= 3 && CARD_POS[ip] < CurCardLen; ip += 2)
		{
			unsigned long long RowHash = MPS_HashCode(Buf + CARD_POS[ip]);
			int RowID;
			if (Map_Row.find(RowHash) == Map_Row.end())
			{
				printf("        Warning: Row Name Not Found, HashCode = %llu, Ignored.\n", RowHash);
				continue;
			}
			RowID = Map_Row[RowHash];
			double RCValue = atof(Buf + CARD_POS[ip + 1]);
			if (RowID == -1) // Objective, Minimization
			{
				Cost[ColID] = RCValue;
				continue;
			}
			if (fabs(RCValue) < Input_Tolerance)
			{
				printf("        Warning: Zero Element Found, %.10lf, Ignored.\n", RCValue);
				continue;
			}
			if (V_Matrix[ColID].find(RowID) != V_Matrix[ColID].end())
			{
				printf("        Warning: Duplicate Element, RowHashCode = %llu, ColHashCode = %llu, Ignored.\n", RowHash, ColHash);
				continue;
			}
			V_Matrix[ColID][RowID] = RCValue;
			n_Element ++;
		}
	}
	printf("        %d Columns and %d Elements Read\n", n_Col, n_Element);
}

void MPS_RHS()
{
	if (strncmp(Buf, "RHS", 3))
	{
		printf("NULL RHS Section!\n");
		return;
	}
	printf("    RHS SECTION\n");
	unsigned long long Enabled_RHS_Hash = MPS_HashCode(Enabled_RHS);
	int n_RHS = 0;
	while (1)
	{
		MPS_ReadLine();
		if (Buf[0] != ' ') // End of RHS
			break;
		if (Enabled_RHS_Hash != 0 && Enabled_RHS_Hash != MPS_HashCode(Buf + CARD_POS[0]))
			continue;
		
		int CurCardLen = strlen(Buf);
		for (int ip = 1; ip <= 3 && CARD_POS[ip] < CurCardLen; ip += 2)
		{
			unsigned long long RowHash = MPS_HashCode(Buf + CARD_POS[ip]);
			int RowID;
			if (Map_Row.find(RowHash) == Map_Row.end())
			{
				printf("        Warning: Row Name Not Found, HashCode = %llu, Ignored.\n", RowHash);
				continue;
			}
			RowID = Map_Row[RowHash];
			double RCValue = atof(Buf + CARD_POS[ip + 1]);
			if (RowID == -1) // Objective, Minimization
			{
				printf("        Warning: Row Name Cannot Be Cost Row, Ignored.\n");
				continue;
			}
			if (fabs(RCValue) < Input_Tolerance)
			{
				printf("        Warning: Zero Element Found, %.10lf, Ignored.\n", RCValue);
				continue;
			}
			V_RHS[RowID] = RCValue;
			n_RHS ++;
		}
	}
	printf("        %d RHS Read\n", n_RHS);
}

void MPS_RANGES()
{
	V_RHS_r = V_RHS;
	if (strncmp(Buf, "RANGES", 6))
	{
		printf("NULL RANGES Section!\n");
		return;
	}
	printf("    RANGES SECTION\n");
	unsigned long long Enabled_RANGES_Hash = MPS_HashCode(Enabled_RANGES);
	int n_RANGES = 0;
	while (1)
	{
		MPS_ReadLine();
		if (Buf[0] != ' ') // End of RANGES
			break;
		if (Enabled_RANGES_Hash != 0 && Enabled_RANGES_Hash != MPS_HashCode(Buf + CARD_POS[0]))
			continue;
		
		int CurCardLen = strlen(Buf);
		for (int ip = 1; ip <= 3 && CARD_POS[ip] < CurCardLen; ip += 2)
		{
			unsigned long long RowHash = MPS_HashCode(Buf + CARD_POS[ip]);
			int RowID;
			if (Map_Row.find(RowHash) == Map_Row.end())
			{
				printf("        Warning: Row Name Not Found, HashCode = %llu, Ignored.\n", RowHash);
				continue;
			}
			RowID = Map_Row[RowHash];
			double Value = atof(Buf + CARD_POS[ip + 1]);
			if (RowID == -1) // Objective, Minimization
			{
				printf("        Warning: Row Name Cannot Be Cost Row, Ignored.\n");
				continue;
			}
			if (fabs(Value) < Input_Tolerance)
			{
				printf("        Warning: Zero Element Found, %.10lf, Ignored.\n", RCValue);
				continue;
			}
			// TODO, Follow lpguide
			if (Row_Type[RowID] == 'E')
			{
				if (Value < 0)
					V_RHS_r[RowID] += Value;
				else
					V_RHS[RowID] += Value;
			}
			else if (Row_Type[RowID] == 'L')
				V_RHS_r[RowID] -= fabs(Value);
			else if (Row_Type[RowID] == 'G')
				V_RHS[RowID] += fabs(Value);
			else
			{
				printf("        Warning: Duplicate RANGES, HashCode = %llu, Ignored.\n", RowHash);
				continue;
			}
			Row_Type[RowID] = 'R'; // Ranged
		}
	}
	printf("        %d RANGES Read\n", n_RANGES);
	// TODO
}

void MPS_BOUNDS()
{
	if (strncmp(Buf, "BOUNDS", 6))
	{
		printf("NULL BOUNDS Section!\n");
		return;
	}
	V_LB.resize(n_Col);
	V_UB.resize(n_Col);
	for (int i = 0; i < n_Col; i ++)
	{
		V_LB[i] = Var_Lower_Bound;
		V_UB[i] = Var_Upper_Bound;
	}
	
	// TODO
}

void MPS_ENDATA()
{
	if (strncmp(Buf, "ENDATA", 6))
		CheckError(1, "MPS_ReadFile: Expected ENDATA!");
}

int MPS_ReadFile()
{
	n_Row = 0;
	n_Col = 0;
	Map_Row.clear();
	Map_Col.clear();
	Obj_Row = 0;
	n_Element = 0;
	Row_Type.clear();
	V_Cost.clear();
	V_Matrix.clear();
	V_RHS.clear();

	printf("ReadMPSFile BEGIN\n");
	Fin = fopen(Filename, "r");
	MPS_ReadLine();
	
	MPS_NAME();
	MPS_ROWS();
	MPS_COLUMNS();
	MPS_RHS();
	MPS_RANGES();
	MPS_BOUNDS();
	MPS_ENDATA();
	
	fclose(FIN);
	printf("ReadMPSFile END\n");
	return 0;
}