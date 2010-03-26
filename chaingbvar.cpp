#include "chaingbvar.h"

double P_SMALLROTATION=0.80;//Probability for small angle crankshaft rotation
							//supposed to create local fluctuation.
double DELTA_TW_K=0.0;//unwinding angle of a single disruption.
double g=0,//(k=150)//g=2.402948861(k=10),k=1,//k segments per Kuhn length.150bp for B-DNA.
	maxRotAng=40.0/180.0*PI,
    const_rotpercentage=0.1;
int maxnum=0;//SHOULD BE CHANGED ACCORDING TO THE LENGTH OF THE CHAIN.
int totsegnum=0;
double bpperseg=0;
int crank_max_length=0;//Crank length [5,20)