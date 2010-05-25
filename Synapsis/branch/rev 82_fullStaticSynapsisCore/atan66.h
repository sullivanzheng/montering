#ifndef atan_66_H
#define atan_66_H
//from: http://www.ganssle.com/approx/approx.pdf

//
// This is the main arctangent approximation "driver"
// It reduces the input argument's range to [0, pi/12],
// and then calls the approximator.
//
//


inline double atan_66s(double x)
{
const double c1=1.6867629106;
const double c2=0.4378497304;
const double c3=1.6867633134;
double x2; // The input argument squared
x2=x * x;
return (x*(c1 + x2*c2)/(c3 + x2));
}

double atan66(double x){
double y; // return from atan__s function
static bool const FALSE=false,TRUE=true;
int complement= FALSE; // true if arg was >1
int region= FALSE; // true depending on region arg is in
int sign= FALSE; // true if arg was < 0
static double const 
		halfpi=1.5707963267948966,
		sixthpi=0.52359877559829882,
		tansixthpi=0.57735026918962573,
		tantwelfthpi=0.2679491924311227;
if (x <0 ){
x=-x;
sign=TRUE; // arctan(-x)=-arctan(x)
}
if (x > 1.0){
x=1.0/x; // keep arg between 0 and 1
complement=TRUE;
}
if (x > tantwelfthpi){
x = (x-tansixthpi)/(1+tansixthpi*x); // reduce arg to under tan(pi/12)
region=TRUE;
}
y=atan_66s(x); // run the approximation
if (region) y+=sixthpi; // correct for region we're in
if (complement)y=halfpi-y; // correct for 1/x if we did that
if (sign)y=-y; // correct for negative arg
return (y);
}


// *********************************************************
// ***
// *** Routines to compute arctangent to 6.6 digits
// *** of accuracy.
// ***
// *********************************************************
//
// atan_66s computes atan(x)
//
// Accurate to about 6.6 decimal digits over the range [0, pi/12].
//
// Algorithm:
// atan(x)= x(c1 + c2*x**2)/(c3 + x**2)
//
#endif