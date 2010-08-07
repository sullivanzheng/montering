#include "chain.h"
#include <cmath>

long CircularChain::getBranchNumber(){
	
	//Max subcript in Fortran convention.
	long jr1=maxnum+1;

	//dLk=enkuhn*282./10.5*sigma 
	//enkuhn is number of kuhn length in the circle.
	double sigma=dLk/(bpperseg*totsegnum/10.5);

	//ek=jr1/enkuhn  which is num of seg per kuhn in my program.
	double ek=282.686/bpperseg;

	//The constant. If curve local writhe is larger than wrc,
	//this part of curve will be considered as a branch.
	static const double wrc=-1.15; //TODO this number will change if large circle is used?
	double s100=100*sigma;

	//Length of the local curve examined.
    static const long loop=long((8.*s100+90)*ek/10.0);

	//Once a peak is detected, the program will jump forward a distance
	//of loopd to get away from this peak, which eliminate multiple counting
	//of the same peak.
	static const long loopd=long(0.667*loop)+1;
	
	double wrt;

	long nl=0;
    long lf=0;
	long il=1;
	long ik;
L871:
	ik=il+loop;
	wrt=this->_wrfun(il-1,ik-1);//il-1 because il is fortran convention here.
	if (wrt < wrc) {
          nl=nl+1;
          if (il<loopd) lf=loopd-il;
          il=il+loopd;
	}
	else{
          il=il+1;
	}
    if (il <= jr1-lf) goto L871;
    if (nl < 2) nl=2;
	return nl;
}

long CircularChain::scanBranch(char* filename){
	using namespace std;
	
	//Max subcript in Fortran convention.
	long jr1=maxnum+1;

	//dLk=enkuhn*282./10.5*sigma 
	//enkuhn is number of kuhn length in the circle.
	double sigma=dLk/(bpperseg*totsegnum/10.5);

	//ek=jr1/enkuhn  which is num of seg per kuhn in my program.
	double ek=282.686/bpperseg;

	//The constant. If curve local writhe is larger than wrc,
	//this part of curve will be considered as a branch.
	static const double wrc=-1.15; //TODO this number will change if large circle is used?
	double s100=100*sigma;

	//Length of the local curve examined.
    static const long loop=long((8.*s100+90)*ek/10.0);

	//Once a peak is detected, the program will jump forward a distance
	//of loopd to get away from this peak, which eliminate multiple counting
	//of the same peak.
	static const long loopd=long(0.667*loop)+1;

	ofstream fp(filename);
	
	double wrt;

	long m, n;
	for (m=0;m<=maxnum;m++){
		n=m+loop;
		wrt=this->_wrfun(m,n);
		fp<<m<<" "<<n<<" "<<wrt<<endl;
	}
	fp.close();
	return 0;
}
double CircularChain::_wrfun(long m, long n)
{	
//Calculate writhing number from m to n.

	/* is,it,jr1 input parameter */
	/* System generated locals */
    long i__1, i__2;
    double ret_val;

    /* Local variables */
    static long i, k;
    static double x[maxa], y[maxa], z[maxa];
    static double f1, f2, f3, f4, ab, ai, ak;
    static long ip, kp;
    static double dx[maxa], dy[maxa], dz[maxa];
    static long is, it;
    static double rx, ry, rz;
    static long jr1;
    static double zn1, zn2, zn3, zn4, abc, cal, abq, sql, www, bllb, blli, bllk;


/*     common /x/ x(jrm),y(jrm),z(jrm) */
/*     common /d/ dx(jrm),dy(jrm),dz(jrm) */
/*     writhing number calculation for the piece [is,it] of the chain. */
/*     cartesian coordinates of chain segments are x(i),y(i),z(i). */
/*     dx, dy, and dz are arrays of bond vectors, ie dx(i) = x(i+1)-x(i). */

	//Importing x and dx
	for (long tempi=0;tempi<=maxnum;tempi++){
		dx[tempi]=C[tempi].dx; dy[tempi]=C[tempi].dy; dz[tempi]=C[tempi].dz;
		x[tempi]=C[tempi].x; y[tempi]=C[tempi].y; z[tempi]=C[tempi].z;
	}
	x[maxnum+1]=dx[maxnum]+x[maxnum];
	y[maxnum+1]=dy[maxnum]+y[maxnum];
	z[maxnum+1]=dz[maxnum]+z[maxnum];
	is=m+1; it=n+1; jr1=maxnum+1;


    www = 0.f;
    i__1 = it - 2;
    for (ip = is; ip <= i__1; ++ip) {
	i = ip;
	if (i > jr1) {
	    i -= jr1;
	}
	i__2 = it;
	for (kp = ip + 2; kp <= i__2; ++kp) {
	    k = kp;
	    if (k > jr1) {
		k -= jr1;
	    }
	    if (k - i > jr1 - 2) {
		goto L32;
	    }
/* L100: */
	    bllb = (x[i - 1] - x[k - 1]) * (dy[k - 1] * dz[i - 1] - dy[
		    i - 1] * dz[k - 1]) + (y[i - 1] - y[k - 1]) * (dz[k - 
		    1] * dx[i - 1] - dz[i - 1] * dx[k - 1]) + (z[i - 
		    1] - z[k - 1]) * (dx[k - 1] * dy[i - 1] - dx[i - 1] 
		    * dy[k - 1]);
	    cal = dx[k - 1] * dx[i - 1] + dy[k - 1] * dy[i - 1] + dz[k - 
		    1] * dz[i - 1];
	    sql = 1.f - cal * cal;
	    if (sql <= 1e-4f) {
		goto L32;
	    }
	    rx = x[k] - x[i];
	    ry = y[k] - y[i];
	    rz = z[k] - z[i];
	    zn1 = sqrt(rx * rx + ry * ry + rz * rz);
	    rx = x[k] - x[i - 1];
	    ry = y[k] - y[i - 1];
	    rz = z[k] - z[i - 1];
	    zn2 = sqrt(rx * rx + ry * ry + rz * rz);
	    rx = x[k - 1] - x[i];
	    ry = y[k - 1] - y[i];
	    rz = z[k - 1] - z[i];
	    zn3 = sqrt(rx * rx + ry * ry + rz * rz);
	    rx = x[k - 1] - x[i - 1];
	    ry = y[k - 1] - y[i - 1];
	    rz = z[k - 1] - z[i - 1];
	    zn4 = sqrt(rx * rx + ry * ry + rz * rz);
	    bllk = (x[i - 1] - x[k - 1]) * dx[k - 1] + (y[i - 1] - y[k - 
		    1]) * dy[k - 1] + (z[i - 1] - z[k - 1]) * dz[k - 1];
	    blli = (x[i - 1] - x[k - 1]) * dx[i - 1] + (y[i - 1] - y[k 
		    - 1]) * dy[i - 1] + (z[i - 1] - z[k - 1]) * dz[
		    i - 1];
	    ak = (bllk - blli * cal) / sql;
	    ai = (bllk * cal - blli) / sql;
	    ab = bllb / sql;
	    abq = ab * ab;
	    abc = abq * cal;
	    f1 = -atan(((1.f - ak) * (1.f - ai) + abc) / (ab * zn1));
	    f2 = -atan(((1.f - ak) * (-ai) + abc) / (ab * zn2));
	    f3 = -atan((-ak * (1.f - ai) + abc) / (ab * zn3));
	    f4 = -atan((-ak * (-ai) + abc) / (ab * zn4));
	    www = www + f1 - f2 - f3 + f4;
L32:
	    ;
	}
/* L31: */
    }
    ret_val = www / PI / 2;
    return ret_val;
} 

