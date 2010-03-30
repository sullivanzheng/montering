#ifndef CHAINGBVAR_H
#define CHAINGBVAR_H	

#include "chainconst.h"

extern double P_SMALLROTATION;//Probability for small angle crankshaft rotation
							//supposed to create local fluctuation.
extern double theta_k,h;//properties of energy curve.
extern double DELTA_TW_K;//unwinding angle of a single disruption.
extern double g,//(k=150)//g=2.402948861(k=10),k=1,
	maxRotAng,//crankshaft
    const_rotpercentage;//halfChain
extern  int maxnum;//SHOULD BE CHANGED ACCORDING TO THE LENGTH OF THE CHAIN.
extern  int totsegnum;
extern double bpperseg;
extern  int crank_max_length;//Crank length [5,20)

#endif /* CHAINGBVARH */