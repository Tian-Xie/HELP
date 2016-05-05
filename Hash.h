#ifndef _HASH_H

#define _HASH_H

#include "LP.h"

const long HASH_OBJECTIVE = -1;
const long HASH_NOT_FOUND = -2;
const int HASH_INSERT_FOUND = 1;

class THashTable
{
private:
	int SIZE;
	int HashMod;
	unsigned long long* Key;
	long* Row;
	long* Next;
	long* Head;
public:
	THashTable();
	int Init(int _LENGTH, int _HashMod);
	int Release();
	long Find(unsigned long long x); // HASH_NOT_FOUND = Not Found, HASH_OBJECTIVE = Objective
	int Insert(unsigned long long x, long _Row);
};

#endif _HASH_H