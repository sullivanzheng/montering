#ifndef CHAINGBVAR_H
#define CHAINGBVAR_H	

#include "chainconst.h"
#include <vector>

extern int maxnum, totsegnum;
extern int crank_min_length,crank_max_length;
extern int VEcutoff;
extern double g_mainchain,g_localchain;//(k=150)//g=2.402948861(k=10),k=1,
extern double bpperunit;
extern double maxRotAng;//crankshaft
extern double P_SMALLROTATION;//Probability for small angle crankshaft rotation
						//supposed to create local fluctuation.
extern double P_REPT; //Probability to make a reptation move.
extern int protect_list[maxa];
extern int reptation_maxlen;
extern int reptation_minlen;
extern int rept_move_range;
#endif /* CHAINGBVARH */