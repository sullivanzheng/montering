#ifndef EZLIN_H
#define EZLIN_H
#include <iostream>
#include <cmath>
#include <vector>
//facility functions for matrix and vector manipulation
/*
inline double dot_product(double a[3],double b[3]);
inline double modu(double a,double b,double c);
inline int mat33addOW(double M1[3][3],double M2[3][3]);
inline int mat33mulvec3(double M[3][3],double v[3],double v_o[3]);
inline int mat33mulscalOW(double M[][3],double a);
inline int mat33disp(double M[3][3]);*/

inline double modu(double a,double b,double c){
	return sqrt(a*a+b*b+c*c);
}

inline double modu2(double a,double b,double c){
	return a*a+b*b+c*c;
}

inline double dot_product(double a[3],double b[3]){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline double betaArray12(double a[3],double b[3]){
	return acos(dot_product(a,b)/modu(a[0],a[1],a[2])/modu(b[0],b[1],b[2]));
}

inline double betaVec1Vec2(std::vector<double> a, std::vector<double> b){
	double _a[3],_b[3];
    for (int i=0;i<3;i++) {_a[i]=a[i];_b[i]=b[i];}
	return betaArray12(_a,_b);
}

inline int mat33addOW(double M1[3][3],double M2[3][3]){
	//add M2 to M1, OW stands for overwrite.
    for(int j=0;j<3;j++) M1[0][j]+=M2[0][j];
    for(int j=0;j<3;j++) M1[1][j]+=M2[1][j];
    for(int j=0;j<3;j++) M1[2][j]+=M2[2][j];
    return 0;
}
inline int mat33mulvec3(double M[3][3],double v[3],double v_o[3]) {
	v_o[0]=dot_product(M[0],v);
	v_o[1]=dot_product(M[1],v);
	v_o[2]=dot_product(M[2],v);
	return 0;
}
inline int mat33mulscalOW(double M[3][3],double a){
	for (int j=0;j<3;j++) M[0][j]*=a;
    for (int j=0;j<3;j++) M[1][j]*=a;
    for (int j=0;j<3;j++) M[2][j]*=a;
return 0;
}

inline int mat33disp(double M[3][3]){
	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			std::cout <<i+1<<','<<j+1<<':'<<M[i][j]<<'\t';
		}
		std::cout <<std::endl;
	}
	return 0;
}

#endif /* EZLIN_H */