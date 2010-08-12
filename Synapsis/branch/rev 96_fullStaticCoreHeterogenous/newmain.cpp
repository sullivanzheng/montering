#include <iostream>
#include "MCbox.h"

long maindo(){
	MCbox_circular sim("_config.cfg");
	sim.performMetropolisCircularCrankRept(1e9);
	return 0;
}

long main(){
	maindo();
	return 0;
}
