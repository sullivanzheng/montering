#include "chaingbvar.h"
int maxnum, totsegnum;
int crank_min_length, crank_max_length;
int VEcutoff;
double g;//(k=150)//g=2.402948861(k=10),k=1,
double bpperseg;
double maxRotAng;//crankshaft
double P_SMALLROTATION;//Probability for small angle crankshaft rotation
std::vector<int> protect_list;
