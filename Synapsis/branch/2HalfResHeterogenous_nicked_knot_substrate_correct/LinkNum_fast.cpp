#include <iostream>
#include <fstream> //For debugging
#include <cmath>
#include "chain.h"
using namespace std;

long CircularChain::productLk_fast(long vertM, long vertN)
//Fast because it will discard overcrossing judgement 
//Since we have no idea of the length of cross over segments, this subroutine should be used with greatest caution.
//In test this subroutine is only marginally faster (<5%) than unaccelerated one.
{	
	if ((vertM<0||vertN<0)||(vertN>maxnum || vertM>maxnum))
	{
		cout<<"[CircularChain::productLk]"
			"Illegal values of m and n ("<<vertM<<','<<vertN<<')';
		exit(EXIT_FAILURE);
	}
	if (vertM>vertN){
		long temp;
		temp=vertM;vertM=vertN;vertN=temp;
	}

	long ret_val, i__1, i__2, i__3;
    double d__1,d__2;
	static const long jrm=910,icrm=910;
	static const double eps=1.0e-5;

    /* Local variables */
    static double d;
    static long k, m, n;
    static double r;
    static double x[jrm], y[jrm], z[jrm];
    static double h1[jrm], h2[jrm], h3[jrm];
    static long j1,n1, j2, n2, i1;
    static double r1, s1, s2, r2;
    static long m1;
    static double da[MAXMatrixDet*MAXMatrixDet];
    static double h11, h21, h31;
    static long kb, l11, l21, n11, mj;
    static double cx[icrm];
    static long ir, mp, ks, ix[icrm], ns, nv;
    static double rx, xr;
    static long ic1[icrm], ic2[icrm], jr3, jr4, jr5;
    static double rl1, rl2;
    static long nv1;
    static double rz1, rz2, dkn;
    static long ipv;

/*	L1 and L2, are the numbers of segments in the first and second contours. 
The coordinates of the first contour are written in in the first (L1+1) elements 
of arrays x, y, z; the elements from L1+2 till L1+L2+2 contain the second contour. Note that
	x((L1+1)=x(1)
	y((L1+1)=y(1)
	z((L1+1)=z(1)
	x((L1+L2+2)=x(L1+2)
	y((L1+L2+2)=y(L1+2)
	z((L1+L2+2)=z(L1+2) */

	long l1=vertM+(maxnum-vertN+1);//Number of seg in the 1st circle
	long l2=vertN-vertM;           //Number of seg in the 2nd circle

//Assume strand exchange happened at the vertices m and n (m<n)
//Therefore the circle will be broken to two circles and they
//are arranged into X in the following fashion:
//C[0~m-1,n~maxnum] maped to X[0~l1-1] and X[l1] to C[0]
//C[m~n-1] maped to X[l1+1~l1+l2] and X[l1+l2+1] to C[m]

	long i,j;
	for (i=0,j=0;i<=vertM-1;i++,j++){
		x[j]=C[i].x;y[j]=C[i].y;z[j]=C[i].z;
	}

	for (i=vertN,j=vertM;i<=maxnum;i++,j++){
		x[j]=C[i].x;y[j]=C[i].y;z[j]=C[i].z;
	}
	
	x[l1]=C[0].x;y[l1]=C[0].y;z[l1]=C[0].z;

	for (i=vertM,j=l1+1;i<=vertN-1;i++,j++){
		x[j]=C[i].x;y[j]=C[i].y;z[j]=C[i].z;
	}
	x[l1+l2+1]=C[vertM].x;
	y[l1+l2+1]=C[vertM].y;
	z[l1+l2+1]=C[vertM].z;

/*//-=--------Debug
	ofstream fp("product.txt");
	for (i=0;i<=l1+l2+1;i++){
		fp<<x[i]<<" "<<y[i]<<" "<<z[i]<<endl;
	}
	fp.close();
//--------------*/


/*      COMMON /x/ X(jrm),Y(jrm),Z(jrm) */
/*      COMMON /F/ DA(200,200) */
    l11 = l1 + 1;
    l21 = l2 + 1;
    jr3 = l1 + l2 - 1;
    jr4 = l1 + l2 + 1;
    jr5 = l1 + l2 + 3;
    ipv = 0;
L297:
    n = 1;
    m = 1;
    i__1 = jr4;
    for (n1 = 1; n1 <= i__1; ++n1) {
	if (n1 == l11) {
	    goto L108;
	}
	s1 = x[n1] - x[n1 - 1];
	s2 = y[n1] - y[n1 - 1];
	h1[n1 - 1] = -s2;
	h2[n1 - 1] = s1;
	h3[n1 - 1] = -x[n1 - 1] * s2 + y[n1 - 1] * s1;
L108:
	;
    }
    i__1 = jr3;
    for (n1 = 1; n1 <= i__1; ++n1) {
	if (n1 == l11) {
	    goto L110;
	}
	j1 = n1 + 2;
	j2 = jr4;
	if (n1 == l1 + 2) {
	    j2 = jr4 - 1;
	}
	h11 = h1[n1 - 1];
	h21 = h2[n1 - 1];
	h31 = h3[n1 - 1];
	i__2 = j2;
	for (n2 = j1; n2 <= i__2; ++n2) {
	    if (n2 == l11) {
		goto L111;
	    }
	    if (n1 != 1) {
		goto L112;
	    }
	    if (n2 == l1) {
		goto L111;
	    }
	    if ((d__1 = x[n1 - 1] - x[n2 - 1], fabs(d__1)) > this->max_seglength * 2 || 
			(d__2 = y[n1 - 1] - y[n2 - 1], fabs(d__2)) > this->max_seglength * 2) {
		goto L111;
	    }
L112:
	    d = h11 * h2[n2 - 1] - h21 * h1[n2 - 1];
	    if (fabs(d) < eps) {
		goto L111;
	    }
	    rx = (h31 * h2[n2 - 1] - h21 * h3[n2 - 1]) / d;
	    if ((x[n1 - 1] - rx) * (x[n1] - rx) >= 0.f) {
		goto L111;
	    }
	    if ((x[n2 - 1] - rx) * (x[n2] - rx) >= 0.f) {
		goto L111;
	    }
	    rz1 = (rx - x[n1 - 1]) / h21 * (z[n1] - z[n1 - 1]) + z[n1 - 
		    1];
	    rz2 = (rx - x[n2 - 1]) / h2[n2 - 1] * (z[n2] - z[n2 - 1]) + 
		    z[n2 - 1];
	    cx[n - 1] = rx;
	    if (rz1 - rz2 >= 0.f) {
		goto L401;
	    }
	    ic1[n - 1] = n1;
	    ic2[n - 1] = n2;
	    if (n1 <= l1) {
		++m;
	    }
	    goto L402;
L401:
	    ic1[n - 1] = n2;
	    ic2[n - 1] = n1;
	    if (n2 <= l1) {
		++m;
	    }
L402:
	    ++n;
	    if (n <= icrm) {
		goto L111;
	    }
	    if (ipv < 2) {
		goto L291;
	    }
	    ret_val = 10000;
	    return ret_val;
L111:
	    ;
	}
L110:
	;
    }
    --m;
    --n;
    ns = n - 1;
    i__1 = ns;
    for (n1 = 1; n1 <= i__1; ++n1) {
	i1 = ic1[n1 - 1];
	nv1 = n1;
	n11 = n1 + 1;
	i__2 = n;
	for (n2 = n11; n2 <= i__2; ++n2) {
	    if ((i__3 = ic1[n2 - 1] - i1) < 0) {
		goto L405;
	    } else if (i__3 == 0) {
		goto L406;
	    } else {
		goto L404;
	    }
L405:
	    i1 = ic1[n2 - 1];
	    nv1 = n2;
	    goto L404;
L406:
	    rl1 = (d__1 = x[i1 - 1] - cx[nv1 - 1], abs(d__1));
	    rl2 = (d__1 = x[i1 - 1] - cx[n2 - 1], abs(d__1));
	    if (rl1 < rl2) {
		goto L404;
	    }
	    nv1 = n2;
L404:
	    ;
	}
	if (n1 >= nv1) {
	    goto L403;
	}
	ir = ic1[n1 - 1];
	ic1[n1 - 1] = ic1[nv1 - 1];
	ic1[nv1 - 1] = ir;
	ir = ic2[n1 - 1];
	ic2[n1 - 1] = ic2[nv1 - 1];
	ic2[nv1 - 1] = ir;
	xr = cx[n1 - 1];
	cx[n1 - 1] = cx[nv1 - 1];
	cx[nv1 - 1] = xr;
L403:
	;
    }
    i__1 = n;
    for (n1 = 1; n1 <= i__1; ++n1) {
	nv = ic2[n1 - 1];
	n2 = 0;
L72:
	++n2;
	if (nv < l1 + 1 && n2 >= m + 1) {
	    goto L711;
	}
	if (n2 >= n + 1) {
	    goto L712;
	}
/* L70: */
	if ((i__2 = ic1[n2 - 1] - nv) < 0) {
	    goto L72;
	} else if (i__2 == 0) {
	    goto L73;
	} else {
	    goto L71;
	}
L73:
	r1 = (d__1 = x[nv - 1] - cx[n1 - 1], abs(d__1));
/* L74: */
	r2 = (d__1 = x[nv - 1] - cx[n2 - 1], abs(d__1));
	if (r1 < r2) {
	    goto L71;
	}
	goto L72;
L711:
	ix[n1 - 1] = 1;
	goto L19;
L712:
	ix[n1 - 1] = m + 1;
	goto L19;
L71:
	ix[n1 - 1] = n2;
L19:
	;
    }
    if (m <= 0) {
	goto L500;
    }
    if (n - m <= 0) {
	goto L500;
    }
L230:
    mj = 0;
    i = 0;
L231:
    ++i;
L251:
    if (i == m) {
	goto L240;
    }
    if (i == n) {
	goto L250;
    }
/* L254: */
    if (ix[i - 1] != i && ix[i - 1] != i + 1) {
	goto L231;
    }
    ++mj;
    --n;
    if (i < m) {
	--m;
    }
    if (m <= 0) {
	goto L500;
    }
    if (n - m <= 0) {
	goto L500;
    }
    i__1 = n;
    for (k = i; k <= i__1; ++k) {
	ix[k - 1] = ix[k];
/* L233: */
    }
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] > i) {
	    --ix[k - 1];
	}
/* L234: */
    }
    goto L251;
L240:
    if (ix[m - 1] != m && ix[m - 1] != 1) {
	goto L231;
    }
    ++mj;
    --n;
    if (m <= 1) {
	goto L500;
    }
    i__1 = n;
    for (k = i; k <= i__1; ++k) {
	ix[k - 1] = ix[k];
/* L252: */
    }
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] == m) {
	    ix[k - 1] = 1;
	}
/* L242: */
    }
    --m;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] > i) {
	    --ix[k - 1];
	}
/* L253: */
    }
    goto L251;
L250:
    if (ix[n - 1] != n && ix[n - 1] != m + 1) {
	goto L232;
    }
    ++mj;
    if (n - m <= 1) {
	goto L500;
    }
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] >= n) {
	    ix[k - 1] = m + 1;
	}
/* L255: */
    }
    --n;
L232:
    if (mj > 0) {
	goto L230;
    }
    if (m <= 1 || n - m <= 1) {
	goto L505;
    }
    i = 0;
L235:
    ++i;
L265:
    if (i == m) {
	++i;
    }
    if (i >= n - 1) {
	goto L266;
    }
    if (ix[i - 1] != ix[i]) {
	goto L235;
    }
    i1 = i + 1;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] == i1) {
	    goto L235;
	}
/* L237: */
    }
    mp = m - 1;
    if (i > m) {
	goto L271;
    }
    if (m <= 2) {
	goto L500;
    }
    m += -2;
    goto L372;
L271:
    if (n - m <= 2) {
	goto L500;
    }
L372:
    n += -2;
    ++mj;
    i__1 = n;
    for (k = i; k <= i__1; ++k) {
	ix[k - 1] = ix[k + 1];
/* L238: */
    }
    if (i != mp) {
	goto L339;
    }
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] == i) {
	    ix[k - 1] = 1;
	}
/* L340: */
    }
L339:
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] > i) {
	    ix[k - 1] += -2;
	}
/* L239: */
    }
    if (m <= 1 || n - m <= 1) {
	goto L230;
    }
    goto L265;
L266:
    if (ix[n - 2] != ix[n - 1]) {
	goto L236;
    }
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] == n) {
	    goto L236;
	}
/* L272: */
    }
    if (n - m <= 2) {
	goto L500;
    }
    ++mj;
    n += -2;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] > n) {
	    ix[k - 1] = m + 1;
	}
/* L273: */
    }
    if (m <= 1 || n - m <= 1) {
	goto L230;
    }
L236:
    if (mj > 0) {
	goto L230;
    }
    i__1 = m;
    for (ks = 1; ks <= i__1; ++ks) {
	if (ix[ks - 1] > m) {
	    goto L490;
	}
/* L503: */
    }
L500:
    dkn = 0.f;
    goto L550;
L505:
    dkn = 1.f;
    goto L550;
L490:
    m1 = m + 1;
    i__1 = n;
    for (ks = m1; ks <= i__1; ++ks) {
	if (ix[ks - 1] < m1) {
	    goto L504;
	}
/* L491: */
    }
    goto L500;
L504:
    --n;
    if (n <= MAXMatrixDet) {
	goto L124;
    }
    if (ipv < 2) {
	goto L291;
    }
    ret_val = 100000;
    return ret_val;
L124:
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = n;
	for (i = 1; i <= i__2; ++i) {
	    da[k + i * MAXMatrixDet - (MAXMatrixDet+1)] = 0.f;
/* L261: */
	}
    }
    i__2 = n;
    for (k = 1; k <= i__2; ++k) {
	i = ix[k - 1];
	da[k + k * MAXMatrixDet - (MAXMatrixDet+1)] = 1.f;
	da[k + i * MAXMatrixDet - (MAXMatrixDet+1)] = -2.f;
	if (k != m) {
	    goto L98;
	}
	da[k - 1] = 1.f;
	goto L262;
L98:
	da[k + (k + 1) * MAXMatrixDet - (MAXMatrixDet+1)] = 1.f;
L262:
	;
    }
    dkn = _det(n, da) / 2.f;
L550:
    ret_val = long(fabs(dkn) + 0.1);
    return ret_val;
L291:
    ++ipv;
    if (ipv >= 2) {
	goto L295;
    }
    i__2 = jr4;
    for (kb = 1; kb <= i__2; ++kb) {
	r = y[kb - 1];
	y[kb - 1] = z[kb - 1];
/* L296: */
	z[kb - 1] = r;
    }
    goto L297;
L295:
    i__2 = jr4;
    for (kb = 1; kb <= i__2; ++kb) {
	r = x[kb - 1];
	x[kb - 1] = z[kb - 1];
/* L298: */
	z[kb - 1] = r;
    }
    goto L297;
}

