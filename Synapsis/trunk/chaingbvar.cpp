#include "chaingbvar.h"
int maxnum, totsegnum;
int crank_min_length, crank_max_length;
int VEcutoff;
double g;//(k=150)//g=2.402948861(k=10),k=1,
double bpperseg;
double maxRotAng;//crankshaft
double P_SMALLROTATION;//Probability for small angle crankshaft rotation
double P_REPT;
int protect_list[maxa];
int reptation_maxlen;
int reptation_minlen;
int rept_move_range;


