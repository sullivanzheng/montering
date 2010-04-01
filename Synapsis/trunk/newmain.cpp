#include <iostream>
#include "MCbox.h"

int maindo(){
	MCbox_circular sim("config.txt",24);
	sim.performMetropolisCircularCrankOnly(1e7);	
	return 0;
}

int main(){
	maindo();
	return 0;
}