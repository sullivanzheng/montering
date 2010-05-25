#include "chaingbvar.h"
long maxnum, totsegnum;
long crank_min_length, crank_max_length;
long VEcutoff;
double g;//(k=150)//g=2.402948861(k=10),k=1,
double bpperseg;
double maxRotAng;//crankshaft
double P_SMALLROTATION;//Probability for small angle crankshaft rotation
double P_REPT;
long protect_list[maxa];
long reptation_maxlen;
long reptation_minlen;
long rept_move_range;


