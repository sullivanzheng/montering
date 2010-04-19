#include <iostream>
#include "MCbox.h"

int maindo(){
	MCbox_circular sim("_config.cfg");
	sim.performMetropolisCircularCrankRept(1e8);	
	return 0;
}

int main(){
	maindo();
	return 0;
}