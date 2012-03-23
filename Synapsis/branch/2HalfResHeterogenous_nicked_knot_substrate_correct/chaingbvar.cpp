#include "chaingbvar.h"
long maxnum, totsegnum;
long crank_min_length, crank_max_length;
long VEcutoff;
double bpperunit;
double maxRotAng;//crankshaft
double P_SMALLROTATION;//Probability for small angle crankshaft rotation
double P_REPT;
double P_REPT_SIMP;
double P_TREADMILL;
long protect_list[maxa];
double specify_rigidity_list[maxa];
long reptation_maxlen;
long reptation_minlen;
long rept_move_range;
double rept_min_seglength;
long VolEx_cutoff_rigidbody;
double VolEx_R_rigid;


long RBAUS_COLLECT_ENABLED;
long RBAUS_LOAD_LAST;
BiasingPotential U(-10.0,30.0,100,100000);

double initial_guess_siteII_umbrella_energy;
double initial_guess_siteI_umbrella_energy;