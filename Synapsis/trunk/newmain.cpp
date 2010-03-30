#include <iostream>
#include "MCbox.h"

int maindo(){
	crank_max_length=0;
    g=7.054815521; //100 seg, 28.571429	seg/Kuhn*4kuhn. 200seg.10.5bp/seg
	MCbox_circular sim("200seg7-10,103-106",200,24);
	sim.performMetropolisCircularCrankOnly(1e7);	
	return 0;
}

int main(){
	maindo();
	return 0;
}