#ifndef CHAINGBVAR_H
#define CHAINGBVAR_H	

#include "chainconst.h"
#include <vector>

extern long maxnum, totsegnum;
extern long crank_min_length,crank_max_length;
extern long VEcutoff;
extern double g;//(k=150)//g=2.402948861(k=10),k=1,
extern double bpperseg;
extern double maxRotAng;//crankshaft
extern double P_SMALLROTATION;//Probability for small angle crankshaft rotation
						//supposed to create local fluctuation.
extern double P_REPT; //Probability to make a reptation move.
extern long protect_list[maxa];
extern long reptation_maxlen;
extern long reptation_minlen;
extern long rept_move_range;
#endif /* CHAINGBVARH */