#ifndef CHAINCONST_H
#define CHAINCONST_H
const double PI=3.1415926535;
const long maxa=910;
const long maxcross=500;

//DNA torsional constant.
//C=2.4~3e-19 erg.cm=2.4~3e-28 J.m. Change this unit to kT.basepairlength
const double C_t = 2.4e-28 / (1.3806503e-23 * 300) / 3.4e-10;
#endif /* CHAINCONST_H */
