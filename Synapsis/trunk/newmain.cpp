#include <iostream>
#include "MCbox.h"

int maindo(){
	crank_max_length=0;
    g=6.161278446; //100 seg, 25seg/Kuhn*4kuhn.
	MCbox_circular sim(100,24);
	sim.performMetropolisCircularCrankOnly(10000000);	
	return 0;
}

int main(){
	maindo();
	return 0;
}