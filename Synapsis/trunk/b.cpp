#include "chain.h"
#include "f2c.h"

double CircularChain::_bwr(int m,int n){
    if ((m<0||n<0)||(n>maxnum || m>maxnum))
	{
		std::cout<<"[in _bwr(int m,int n)]Illegal values of m and n ("<<m<<','<<n<<')';
		exit(EXIT_FAILURE);
	}
	if (m>=n)
	{
		std::cout<<"[in _bwr(int m,int n)]Illegal values of m and n, m should be smaller than n ("<<m<<','<<n<<')';
		exit(EXIT_FAILURE);
	}

	/* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, m__;
    static double w[maxa];
    static double a1, a2;
    static int m1;
    static double t1, t2;
    static int ii, ik, in;
    static int mm;
    static double dx[maxa], dy[maxa], dz[maxa];
    static double si;
    static double sm, rp;
    static int jr1;
    static double sm1;
    static double pi12;
    static double rxy;
    static int ibad[maxa];
    static double cchi;
    static double beta, cphi[maxa];
    static double schi;
    static double seta, sphi[maxa], rxyp, rxyp2, perpx, perpy, perpz;
    static double perpzr;



/* Function parameter */
/* The returning parameter */
/*     *  eps=1.d-9, eps1=1.d-11, one=1.d0-1.d-12 these are actually constants. */
/*      common /x/ x(maxa),y(maxa),z(maxa) */
/*      common /d/ dx(maxa),dy(maxa),dz(maxa) */
/*      common /w/ w(maxa) */
/* convert constants to variables and change them to const in C++. */
    const static double pi2 = PI * 2, pil2 = PI / 2;
	const static double eps = 1e-9, eps1 = 1e-11,  one = .99999999999900002;
/* ------------------------------------------------ */

/* import C[i].dX to dX*///--------------------------------
	for (int tempi=0;tempi<=maxnum;tempi++){
		dx[tempi]=C[tempi].dx; dy[tempi]=C[tempi].dy; dz[tempi]=C[tempi].dz;		
	}
	for (int tempi=0;tempi<=maxnum;tempi++){
		w[tempi]=0;//TODO w or as we can see from the program, it is not necessary.
	}
	jr1=maxnum+1;
	in=m+1;ik=n+1;
//----------------------------------------------------------------

    beta = 0.f;
    i__1 = ik + 1;
    for (ii = in - 1; ii <= i__1; ++ii) {
	i__ = ii;
	if (i__ == 0) {
	    i__ = jr1;
	}
	if (i__ > jr1) {
	    i__ -= jr1;
	}
	rxy = sqrt(dx[i__ - 1] * dx[i__ - 1] + dy[i__ - 1] * dy[i__ - 1]);
	if (rxy >= eps) {
	    cphi[i__ - 1] = dx[i__ - 1] / rxy;
	    sphi[i__ - 1] = dy[i__ - 1] / rxy;
	    ibad[i__ - 1] = 0;
	} else {
	    cphi[i__ - 1] = 1.f;
	    sphi[i__ - 1] = 0.f;
	    ibad[i__ - 1] = 1;
	}
/* L11: */
    }
/* $DIR NO_RECURRENCE */
    i__1 = ik;
    for (mm = in - 1; mm <= i__1; ++mm) {
	m__ = mm;
	m1 = m__ + 1;
	if (m__ == 0) {
	    m__ = jr1;
	}
	if (m__ > jr1) {
	    m__ -= jr1;
	}
	if (m1 > jr1) {
	    m1 -= jr1;
	}
	if (ibad[m__ - 1] == 1) {
	    perpx = -dz[m__ - 1] * dy[m1 - 1];
	    perpy = dz[m__ - 1] * dx[m1 - 1] - eps * dz[m1 - 1];
	    perpz = eps * dy[m1 - 1];
	} else if (ibad[m1 - 1] == 1) {
	    perpx = dy[m__ - 1] * dz[m1 - 1];
	    perpy = dz[m__ - 1] * eps - dx[m__ - 1] * dz[m1 - 1];
	    perpz = -dy[m__ - 1] * eps;
	} else {
	    perpx = dy[m__ - 1] * dz[m1 - 1] - dz[m__ - 1] * dy[m1 - 1];
	    perpy = dz[m__ - 1] * dx[m1 - 1] - dx[m__ - 1] * dz[m1 - 1];
	    perpz = dx[m__ - 1] * dy[m1 - 1] - dy[m__ - 1] * dx[m1 - 1];
	}
	if (abs(perpz) < eps1) {
	    w[m1 - 1] = 0.f;
	    goto L1;
	}
	perpzr = perpz;
	si = perpzr>=0?1.0:-1.0;
	rxyp2 = perpx * perpx + perpy * perpy;
	rxyp = sqrt(rxyp2);
	if (rxyp < eps1) {
	    w[m1 - 1] = 0.f;
	    goto L1;
	}
	rp = sqrt(rxyp2 + perpz * perpz);
	seta = rxyp / rp;
	cchi = perpx / rxyp * si;
	schi = perpy / rxyp * si;
	sm1 = sphi[m1 - 1] * cchi - cphi[m1 - 1] * schi;
	sm = sphi[m__ - 1] * cchi - cphi[m__ - 1] * schi;
	a1 = seta * sm1;
	if (a1 > one) {
	    t1 = pi12;
	} else if (a1 < -one) {
	    t1 = -pi12;
	} else {
	    t1 = asin(a1);
	}
	a2 = seta * sm;
	if (a2 > one) {
	    t2 = pi12;
	} else if (a2 < -one) {
	    t2 = -pi12;
	} else {
	    t2 = asin(a2);
	}
	w[m1 - 1] = t1 - t2;
	beta += w[m1 - 1];
L1:
	;
    }
    beta /= pi2;
    return beta;
} /* bwr_ */