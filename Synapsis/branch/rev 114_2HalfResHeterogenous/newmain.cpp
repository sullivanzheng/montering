#include <iostream>
#include "MCbox.h"
#include <cmath>


long maindo(){
	MCbox_circular sim("_config.cfg");
	sim.performMetropolisCircularCrankRept(1e9);
	return 0;
}

inline double tempE(double x){
	if (x<-5 || x>5)
		return 3*(fabs(x)-5)*(fabs(x)-5)+10;
	return 15*sin(x);
}
int maindo1(){
	MTRand53 mt(1234);
	BiasingPotential Ut(-5,5,100);
	unsigned long n=0;
	double x=0;
	double E=tempE(x)+ Ut.getBiasingE(x);
	for (n=0;n<1000000;n++){
		double dx=(mt()-.5)*2.*0.1;
		double newx=x+dx;
		double newE=tempE(newx)+ Ut.getBiasingE(newx);
		if (newE<E){
			x=newx;E=newE;
		}
		else{
			double temp=mt();
			if (temp<exp(-newE+E)){
				x=newx;
				E=newE;
			}
		}
		Ut.collect(x);
		if (n % 100000==0)
			Ut.pickle("testAF.txt");
	}
	return 0;
}

int main(){
	maindo();
	return 0;
}


