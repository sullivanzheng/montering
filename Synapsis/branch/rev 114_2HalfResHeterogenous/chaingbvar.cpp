#include "chaingbvar.h"
long maxnum, totsegnum;
long crank_min_length, crank_max_length;
long VEcutoff;
double bpperunit;
double maxRotAng;//crankshaft
double P_SMALLROTATION;//Probability for small angle crankshaft rotation
double P_REPT;
long protect_list[maxa];
double specify_rigidity_list[maxa];
long reptation_maxlen;
long reptation_minlen;
long rept_move_range;
long VolEx_cutoff_rigidbody;


