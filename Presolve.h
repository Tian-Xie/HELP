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
 
#ifndef _PRESOLVE_H

#define _PRESOLVE_H

extern int Presolve_Modified;
extern double Row_1Norm[MAX_ROWS], Col_1Norm[MAX_COLS];
extern int Row_Disable[MAX_ROWS], Col_Disable[MAX_COLS];
extern int Row_Element_Count[MAX_ROWS], Col_Element_Count[MAX_COLS];
extern int Presolve_Linked_List_Head, Presolve_Linked_List_Tail, Presolve_Linked_List_Next[MAX_ROWS + MAX_COLS];

#endif
