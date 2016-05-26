/*****************************************************************************
*  LPSolver, An experimental implementation of Homogeneous and Self-Dual     *
*  Algorithm for Linear Programming.                                         *
*  Author: Tian Xie (Research Center for Management Science and Information  *
*          Analytics, Shanghai University of Finance and Economics)          *
*  Credits: (1) Fundamental implementation idea originated from COPL_LP.     *
*               (Xiong Zhang and Yinyu Ye)                                   *
*               See http://web.stanford.edu/~yyye/Col.html .                 *
*           (2) Sparse Cholesky Decomposition is supported by CHOLMOD.       *
*               (Timothy A. Davis)                                           *
******************************************************************************/

#include "LP.h"
#include "Hash.h"

THashTable::THashTable()
{
	this -> SIZE = 0;
	this -> HashMod = 999983;
	this -> Key = 0;
	this -> Row = 0;
	this -> Next = 0;
	this -> Head = 0;
}

int THashTable::Init(int _LENGTH, int _HashMod)
{
	Release();
	this -> Key = new unsigned long long[_LENGTH];
	this -> Row = new long[_LENGTH];
	this -> Next = new long[_LENGTH];
	this -> Head = new long[_HashMod];
	this -> HashMod = _HashMod;
	this -> SIZE = 0;
	memset(this -> Head, -1, sizeof(long) * _HashMod);
	return 0;
}

int THashTable::Release()
{
	if (this -> Key)
		delete this -> Key;
	this -> Key = 0;
	if (this -> Row)
		delete this -> Row;
	this -> Row = 0;
	if (this -> Next)
		delete this -> Next;
	this -> Next = 0;
	if (this -> Head)
		delete this -> Head;
	this -> Head = 0;
	this -> SIZE = 0;
	return 0;
}

long THashTable::Find(unsigned long long x)// -2 = Not Found, -1 = Objective
{
	for (int i = this -> Head[x % this -> HashMod]; i != -1; i = this -> Next[i])
		if (this -> Key[i] == x)
			return this -> Row[i];
	return HASH_NOT_FOUND;
}

int THashTable::Insert(unsigned long long x, long _Row)
{
	long HashKey = x % this -> HashMod;
	for (int i = this -> Head[HashKey]; i != -1; i = this -> Next[i])
		if (this -> Key[i] == x)
			return HASH_INSERT_FOUND; // Found
	this -> Key[this -> SIZE] = x;
	this -> Row[this -> SIZE] = _Row;
	this -> Next[this -> SIZE] = this -> Head[HashKey];
	this -> Head[HashKey] = this -> SIZE;
	this -> SIZE ++;
	return 0;
}
