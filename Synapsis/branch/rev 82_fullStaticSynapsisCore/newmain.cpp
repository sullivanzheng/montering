#include <iostream>
#include "MCbox.h"

long maindo(){
	MCbox_circular sim("_config.cfg");
	sim.performMetropolisCircularCrankRept(1e8);
	return 0;
}

long main(){
	//maindo();
	using namespace std;
	for (long i=0;i<100;i++){
		cout<<drand(1.0)<<" ";
	}
}
