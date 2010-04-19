#include <iostream>
#include "MCbox.h"

int maindo(){
	MCbox_circular sim("_config.txt");
	sim.performMetropolisCircularCrankOnly(1e8);	
	return 0;
}

int main(){
	maindo();
	return 0;
}