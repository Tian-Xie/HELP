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
#include "Hash.h"

const int CARD_POS[] = {4, 14, 24, 39, 49};
FILE* Fin;
const int MAX_LINE_LENGTH = 256;
char Buf[MAX_LINE_LENGTH];
unsigned long long Obj_Row;
THashTable Hash_Row, Hash_Col;

void MPS_CopyName(char* Dest, char* Src)
{
	int Lenm1 = min((int) strlen(Src), 8) - 1;
	while (Lenm1 >= 0 && (Src[Lenm1] == ' ' || Src[Lenm1] == '\r' || Src[Lenm1] == '\n'))
		Lenm1 --;
	Dest[Lenm1 + 1] = 0;
	for (; Lenm1 >= 0; Lenm1 --)
		Dest[Lenm1] = Src[Lenm1];
}

unsigned long long MPS_HashCode(char* S)
{
	int Lenm1 = min((int) strlen(S), 8) - 1;
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
		if (feof(Fin))
			CheckError(1, "MPS_ReadLine: EOF before ENDATA!");
		if (fgets(Buf, MAX_LINE_LENGTH, Fin) == NULL)
			CheckError(1, "MPS_ReadLine: fgets ERROR!");
		if (Buf[0] != '*' && Buf[0] != '\r' && Buf[0] != '\n')
			break;
	}
	for (int i = 0; Buf[i]; i ++)
		if (Buf[i] == '\r' || Buf[i] == '\n')
		{
			Buf[i] = 0;
			break;
		}
}

void MPS_NAME()
{
	if (strncmp(Buf, "NAME", 4))
		CheckError(1, "MPS_NAME: Expected NAME!");
	MPS_CopyName(Problem_Name, Buf + 14);
	printf("    NAME SECTION\n");
	printf("        NAME = \"%s\"\n", Problem_Name);
	MPS_ReadLine(); // For MPS_ROWS
}

void MPS_ROWS()
{
	if (strncmp(Buf, "ROWS", 4))
		CheckError(1, "MPS_ROWS: Expected ROWS!");
	printf("    ROWS SECTION\n");
	while (1)
	{
		MPS_ReadLine();
		if (Buf[0] != ' ') // End of ROWS
			break;
		char Type = toupper(Buf[1]);
		if (Type != 'N' && Type != 'E' && Type != 'G' && Type != 'L')
		{
			printf("        Warning: Unknown Type: \"%c\", Ignored.\n", Type);
			continue;
		}
		unsigned long long Hash = MPS_HashCode(Buf + CARD_POS[0]);
		if (Type == 'N')
		{
			if (Obj_Row)
			{
				printf("        Warning: Duplicate Objective Row, Ignored.\n");
				continue;
			}
			Obj_Row = Hash; // Store Objective Row's Hash Value
			Hash_Row.Insert(Hash, HASH_OBJECTIVE); // Objective
			continue;
		}
		if (Hash_Row.Insert(Hash, n_Row) == HASH_INSERT_FOUND)
		{
			printf("        Warning: Duplicate Row Name, HashCode = %llu, Ignored.\n", Hash);
			continue;
		}
		Row_Type[n_Row] = Type;
		V_RHS[n_Row] = 0.0;
		n_Row ++;
	}
	printf("        %d Rows Read (Objective Row Excluded)\n", n_Row);
	if (! Obj_Row)
		CheckError(1, "MPS_ROWS: No Objective Row!");
}

void MPS_COLUMNS()
{
	if (strncmp(Buf, "COLUMNS", 7))
		CheckError(1, "MPS_COLUMNS: Expected COLUMNS!");
	printf("    COLUMNS SECTION\n");
	while (1)
	{
		MPS_ReadLine();
		if (Buf[0] != ' ') // End of COLUMNS
			break;
		unsigned long long ColHash = MPS_HashCode(Buf + CARD_POS[0]);
		int ColID = Hash_Col.Find(ColHash);
		if (ColID == HASH_NOT_FOUND)
		{
			Hash_Col.Insert(ColHash, n_Col);
			ColID = n_Col;
			V_Cost[n_Col] = 0.0;
			V_Matrix_Col_Head[n_Col] = -1;
			n_Col ++;
		}

		int CurCardLen = strlen(Buf);
		for (int ip = 1; ip <= 3 && CARD_POS[ip] < CurCardLen; ip += 2)
		{
			unsigned long long RowHash = MPS_HashCode(Buf + CARD_POS[ip]);
			int RowID = Hash_Row.Find(RowHash);
			if (RowID == HASH_NOT_FOUND)
			{
				printf("        Warning: Row Name Not Found, HashCode = %llu, Ignored.\n", RowHash);
				continue;
			}
			double RCValue = atof(Buf + CARD_POS[ip + 1]);
			if (RowID == HASH_OBJECTIVE) // Objective, Minimization
			{
				V_Cost[ColID] = RCValue;
				continue;
			}
			if (fabs(RCValue) < Input_Tolerance)
			{
				printf("        Warning: Zero Element Found, %.10lf, Ignored.\n", RCValue);
				continue;
			}
			// Cannot check duplicate element
			V_Matrix_Row[n_Element] = RowID;
			V_Matrix_Value[n_Element] = RCValue;
			V_Matrix_Col_Next[n_Element] = V_Matrix_Col_Head[ColID];
			V_Matrix_Col_Head[ColID] = n_Element;
			n_Element ++;
		}
	}
	printf("        %d Columns and %d Elements Read (Objective Row Excluded)\n", n_Col, n_Element);
}

void MPS_RHS()
{
	if (strncmp(Buf, "RHS", 3))
	{
		printf("    NULL RHS Section!\n");
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
			int RowID = Hash_Row.Find(RowHash);
			if (RowID == HASH_NOT_FOUND)
			{
				printf("        Warning: Row Name Not Found, HashCode = %llu, Ignored.\n", RowHash);
				continue;
			}
			double RCValue = atof(Buf + CARD_POS[ip + 1]);
			if (RowID == HASH_OBJECTIVE) // Objective, Minimization
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
	for (int i = 0; i < n_Row; i ++)
		V_RHS_r[i] = V_RHS[i];
	if (strncmp(Buf, "RANGES", 6))
	{
		printf("    NULL RANGES Section!\n");
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
			int RowID = Hash_Row.Find(RowHash);
			if (RowID == HASH_NOT_FOUND)
			{
				printf("        Warning: Row Name Not Found, HashCode = %llu, Ignored.\n", RowHash);
				continue;
			}
			double Value = atof(Buf + CARD_POS[ip + 1]);
			if (RowID == HASH_OBJECTIVE) // Objective, Minimization
			{
				printf("        Warning: Row Name Cannot Be Cost Row, Ignored.\n");
				continue;
			}
			if (fabs(Value) < Input_Tolerance)
			{
				printf("        Warning: Zero Element Found, %.10lf, Ignored.\n", Value);
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
}

void MPS_BOUNDS()
{
	for (int i = 0; i < n_Col; i ++)
	{
		V_LB[i] = Var_Lower_Bound;
		V_UB[i] = Var_Upper_Bound;
	}
	
	if (strncmp(Buf, "BOUNDS", 6))
	{
		printf("    NULL BOUNDS Section!\n");
		return;
	}
	printf("    BOUNDS SECTION\n");
	unsigned long long Enabled_BOUNDS_Hash = MPS_HashCode(Enabled_BOUNDS);
	int n_BOUNDS = 0;
	while (1)
	{
		MPS_ReadLine();
		if (Buf[0] != ' ') // End of BOUNDS
			break;
		if (Enabled_BOUNDS_Hash != 0 && Enabled_BOUNDS_Hash != MPS_HashCode(Buf + CARD_POS[0]))
			continue;
		
		unsigned long long ColHash = MPS_HashCode(Buf + CARD_POS[1]);
		int ColID = Hash_Col.Find(ColHash);
		if (ColID == HASH_NOT_FOUND)
		{
			printf("        Warning: Column Name Not Found, HashCode = %llu, Ignored.\n", ColHash);
			continue;
		}
		double Value = 0.0;
		if (strlen(Buf) > 24)
			Value = atof(Buf + CARD_POS[2]);

		if (Buf[1] == 'L' && Buf[2] == 'O' && Buf[3] == ' ')
			V_LB[ColID] = Value;
		else if (Buf[1] == 'U' && Buf[2] == 'P' && Buf[3] == ' ')
			V_UB[ColID] = Value;
		else if (Buf[1] == 'F' && Buf[2] == 'X' && Buf[3] == ' ')
			V_LB[ColID] = V_UB[ColID] = Value;
		else if (Buf[1] == 'F' && Buf[2] == 'R' && Buf[3] == ' ')
		{
			V_LB[ColID] = -MaxPositive;
			V_UB[ColID] = MaxPositive;
		}
		else if (Buf[1] == 'M' && Buf[2] == 'I' && Buf[3] == ' ')
		{
			V_LB[ColID] = -MaxPositive;
			V_UB[ColID] = Value;
		}
		else if (Buf[1] == 'P' && Buf[2] == 'L' && Buf[3] == ' ')
		{
			V_LB[ColID] = Value;
			V_UB[ColID] = MaxPositive;
		}
		n_BOUNDS ++;
	}
	printf("        %d BOUNDS Read\n", n_BOUNDS);
}

void MPS_ENDATA()
{
	if (strncmp(Buf, "ENDATA", 6))
		CheckError(1, "MPS_ENDATA: Expected ENDATA!");
}

int MPS_ReadFile()
{
	n_Row = 0;
	n_Col = 0;
	Hash_Row.Init(MAX_ROWS, 999983);
	Hash_Col.Init(MAX_COLS, 999983);
	Obj_Row = 0;
	n_Element = 0;
	V_Cost_Intercept = 0;
	
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
	
	fclose(Fin);
	Hash_Row.Release();
	Hash_Col.Release();
	printf("ReadMPSFile END\n");
	return 0;
}
