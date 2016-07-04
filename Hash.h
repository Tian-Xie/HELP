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

#ifndef _HASH_H

#define _HASH_H

#pragma once

#include "LP.h"

const int HASH_OBJECTIVE = -1;
const int HASH_NOT_FOUND = -2;
const int HASH_INSERT_FOUND = 1;

class THashTable
{
private:
	int SIZE;
	int HashMod;
	int* Row;
	int* Next;
	int* Head;
	unsigned long long* Key;
public:
	THashTable();
	int Init(int _LENGTH, int _HashMod);
	int Release();
	int Find(unsigned long long x); // HASH_NOT_FOUND = Not Found, HASH_OBJECTIVE = Objective
	int Insert(unsigned long long x, int _Row);
};

#endif
