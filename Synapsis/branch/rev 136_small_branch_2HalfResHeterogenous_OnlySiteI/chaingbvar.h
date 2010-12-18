#ifndef CHAINGBVAR_H
#define CHAINGBVAR_H	

#include "chainconst.h"
#include "bias.h"
#include <vector>

extern long maxnum, totsegnum;
extern long crank_min_length,crank_max_length;
extern long VEcutoff;
extern double bpperunit;
extern double maxRotAng;//crankshaft
extern double P_SMALLROTATION;//Probability for small angle crankshaft rotation
						//supposed to create local fluctuation.
extern double P_REPT; //Probability to make a reptation move.
extern double P_REPT_SIMP;
extern double P_TREADMILL;//Probability to make a treadmill move.
extern long protect_list[maxa];
extern double specify_rigidity_list[maxa];
extern long reptation_maxlen;
extern long reptation_minlen;
extern long rept_move_range;
extern double rept_min_seglength;
extern long VolEx_cutoff_rigidbody;
extern BiasingPotential U;

extern long RBAUS_COLLECT_ENABLED;
extern long RBAUS_LOAD_LAST;

#endif /* CHAINGBVARH */