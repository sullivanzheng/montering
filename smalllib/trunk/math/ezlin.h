#ifndef EZLIN_H
#define EZLIN_H
#include <iostream>
#include <cmath>
#include <vector>
//facility functions for matrix and vector manipulation
/*
inline double dot_product(double a[3],double b[3]);
inline double modu(double a,double b,double c);
inline long mat33addOW(double M1[3][3],double M2[3][3]);
inline long mat33mulvec3(double M[3][3],double v[3],double v_o[3]);
inline long mat33mulscalOW(double M[][3],double a);
inline long mat33disp(double M[3][3]);*/

inline double modu(double a,double b,double c){
	return sqrt(a*a+b*b+c*c);
}

inline double modu2(double a,double b,double c){
	return a*a+b*b+c*c;
}

inline double modu2V(double a[3]){
	return a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
}

inline double dot_product(double a[3],double b[3]){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline double betaArray12(double a[3],double b[3]){
	return acos(dot_product(a,b)/modu(a[0],a[1],a[2])/modu(b[0],b[1],b[2]));
}

inline double betaVec1Vec2(std::vector<double> a, std::vector<double> b){
	double _a[3],_b[3];
    for (long i=0;i<3;i++) {_a[i]=a[i];_b[i]=b[i];}
	return betaArray12(_a,_b);
}

inline long mat33addOW(double M1[3][3],double M2[3][3]){
	//add M2 to M1, OW stands for overwrite.
    for(long j=0;j<3;j++) M1[0][j]+=M2[0][j];
    for(long j=0;j<3;j++) M1[1][j]+=M2[1][j];
    for(long j=0;j<3;j++) M1[2][j]+=M2[2][j];
    return 0;
}
inline long mat33mulvec3(double M[3][3],double v[3],double v_o[3]) {
	v_o[0]=dot_product(M[0],v);
	v_o[1]=dot_product(M[1],v);
	v_o[2]=dot_product(M[2],v);
	return 0;
}
inline long mat33mulscalOW(double M[3][3],double a){
	for (long j=0;j<3;j++) M[0][j]*=a;
    for (long j=0;j<3;j++) M[1][j]*=a;
    for (long j=0;j<3;j++) M[2][j]*=a;
return 0;
}

inline long mat33disp(double M[3][3]){
	for (long i=0;i<3;i++){
		for (long j=0;j<3;j++){
			std::cout <<i+1<<','<<j+1<<':'<<M[i][j]<<'\t';
		}
		std::cout <<std::endl;
	}
	return 0;
}

inline long Xprod(double a[3], double b[3], double o[3]){
	o[0]= a[1]*b[2]-b[1]*a[2];
	o[1]=-a[0]*b[2]+b[0]*a[2];
	o[2]= a[0]*b[1]+b[0]*a[1];
	return 0;
}

inline long Xprod_exp(double a0, double a1, double a2,
						 double b0, double b1, double b2,
						 double o[3]){
	o[0]= a1*b2-b1*a2;
	o[1]=-a0*b2+b0*a2;
	o[2]= a0*b1+b0*a1;
	return 0;
}

#endif /* EZLIN_H */