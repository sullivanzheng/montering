#include "chain.h"
#include "f2c.h"

//Calculate integer part of writhe and AlexPoly(-1)
int CircularChain::_kndwr_topl_update(double &topl,int &ierr){
	const int isi=801;
	const int range1=4001;

    /* System generated locals */
    int i__1, i__2, i__3;
    double r__1;
    double d__1, d__2;
    /* Local variables */
    static double c__;
    static double d__;
    static int i__, j, k, l;
    static double x[maxa], y[maxa], z[maxa];
    static int i1, n1, n2, m4, n4;
	
    static double r1, r2, da[isi*isi];
    static int n11, n21, n41, jj;
    static double dr;
    static int mj;
    static double cx[range1];
    static double dx[maxa], dy[maxa], dz[maxa];
    static int ir, ks, ix[range1], js, nv;
    static double rx, xr;
    static int ic1[range1], ic2[range1], jr1, jr2, jr3, jr4;
    static double rl1, rl2;
    static int nv1;
    static double px1, py1, px2, py2;
    static double rz1, rz2;
    static double  px11, px21, drx;
    static int jwr;
    static double zet;
    static double pdx1, pdy1, pdx2, pdy2;
    static int jmin;
    static double pdx12;
    static int kmax, iint, iiwr;

    const double deps = 1e-10;
    const double eps = 1e-4f;

/* import C[i].dX to dX*///--------------------------------
	for (int tempi=0;tempi<=maxnum;tempi++){
		x[tempi]=C[tempi].x; y[tempi]=C[tempi].y; z[tempi]=C[tempi].z;		
		dx[tempi]=C[tempi].dx; dy[tempi]=C[tempi].dy; dz[tempi]=C[tempi].dz;		
	}
	jr1=maxnum+1;
	x[maxnum+1]=C[0].x;y[maxnum+1]=C[0].y;z[maxnum+1]=C[0].z;
	ierr=0;
//_________________________________________________________

/* inventorization of intersections */
/* n4 - number of intersections */
/* cx(i) - x-value of i-th intersection */
/* ic1(i) - number of undergoing segment for i-th intersection */
/* ic2(i) - number of overgoing segment for i-th intersection */
    jr2 = jr1 - 2;
    n4 = 1;
    i__1 = jr2;
    for (n1 = 1; n1 <= i__1; ++n1) {
	n11 = n1 + 1;
	jr3 = jr1;
	if (n1 < 2) {
	    --jr3;
	}
	jr4 = n1 + 2;
	pdx1 = dx[n1 - 1];
	pdy1 = dy[n1 - 1];
	px1 = x[n1 - 1];
	py1 = y[n1 - 1];
	px11 = x[n11 - 1];
	i__2 = jr3;
	for (n2 = jr4; n2 <= i__2; ++n2) {
	    px2 = x[n2 - 1];
	    py2 = y[n2 - 1];
	    if ((d__1 = px1 - px2, abs(d__1)) > 2.01) continue;
		else if ((d__2 = py1 - py2, abs(d__2)) > 2.01) continue;
	    pdx2 = dx[n2 - 1];
	    pdy2 = dy[n2 - 1];
	    n21 = n2 + 1;
	    px21 = x[n21 - 1];
	    d__ = pdx2 * pdy1 - pdx1 * pdy2;
	    if (abs(d__) < deps) {
		goto L111;
	    }
	    pdx12 = pdx1 * pdx2;
	    drx = (py2 * pdx12 - px2 * pdx1 * pdy2 - py1 * pdx12 + px1 * pdx2 
		    * pdy1) / d__;
	    if (px1 <= px11) {
		if (drx < px1 || drx >= px11) {
		    goto L111;
		}
	    } else {
		if (drx <= px11 || drx > px1) {
		    goto L111;
		}
	    }
	    if (px2 <= px21) {
		if (drx < px2 || drx >= px21) {
		    goto L111;
		}
	    } else {
		if (drx <= px21 || drx > px2) {
		    goto L111;
		}
	    }
	    rx = drx;
	    if ((double) n4 > range1) {
		goto L1002;
	    }
	    rz1 = (rx - x[n1 - 1]) / dx[n1 - 1] * dz[n1 - 1] + z[n1 - 1];
	    rz2 = (rx - x[n2 - 1]) / dx[n2 - 1] * dz[n2 - 1] + z[n2 - 1];
	    cx[n4 - 1] = rx;
	    if (rz1 - rz2 >= 0.) {
		goto L401;
	    }
	    ic1[n4 - 1] = n1;
	    ic2[n4 - 1] = n2;
	    goto L402;
L401:
	    ic1[n4 - 1] = n2;
	    ic2[n4 - 1] = n1;
L402:
	    ++n4;
L111:
	    ;
	}
/* L110: */
    }
    --n4;
/* calculation of the directional writhing number */
    jwr = 0;
    i__1 = n4;
    for (iint = 1; iint <= i__1; ++iint) {
	zet = dx[ic1[iint - 1] - 1] * dy[ic2[iint - 1] - 1] - dx[ic2[iint - 1]
		 - 1] * dy[ic1[iint - 1] - 1];
	if (zet < 0.0) {
	    iiwr = 1;
	} else {
	    iiwr = -1;
	}
	jwr += iiwr;
/* L1301: */
    }
/* renumeration of intersections in the proper order */
    if (n4 <= 2) {
	goto L31;
    }
    n41 = n4 - 1;
    i__1 = n41;
    for (n1 = 1; n1 <= i__1; ++n1) {
	i1 = ic1[n1 - 1];
	nv1 = n1;
	n11 = n1 + 1;
	i__2 = n4;
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
/*   determination of the number ix(i) of overgoing  part of line */
/*   for i-th intersection */
    i__1 = n4;
    for (n1 = 1; n1 <= i__1; ++n1) {
	nv = ic2[n1 - 1];
	n2 = 0;
L72:
	++n2;
	if (n2 >= n4 + 1) {
	    goto L771;
	}
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
L771:
	ix[n1 - 1] = 1;
	goto L19;
L71:
	ix[n1 - 1] = n2;
L19:
	;
    }
/* attempts for some detanglement */
L230:
    mj = 0;
    i__ = 0;
L231:
    ++i__;
    if (i__ >= n4) {
	goto L240;
    }
    if (ix[i__ - 1] != i__ && ix[i__ - 1] != i__ + 1) {
	goto L231;
    }
    ++mj;
    --n4;
    if (n4 <= 2) {
	goto L31;
    }
    i__1 = n4;
    for (k = i__; k <= i__1; ++k) {
	ix[k - 1] = ix[k];
/* L233: */
    }
    i__1 = n4;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] > i__) {
	    --ix[k - 1];
	}
/* L234: */
    }
    --i__;
    goto L231;
L240:
    if (ix[n4 - 1] != n4 && ix[n4 - 1] != 1) {
	goto L232;
    }
    i__1 = n4;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] >= n4) {
	    ix[k - 1] = 1;
	}
/* L242: */
    }
    ++mj;
    --n4;
    if (n4 <= 2) {
	goto L31;
    }
L232:
    if (mj > 0) {
	goto L230;
    }
    i__ = 0;
L235:
    ++i__;
    if (i__ >= n4 - 1) {
	goto L244;
    }
    if (ix[i__ - 1] != ix[i__]) {
	goto L235;
    }
    i1 = i__ + 1;
    i__1 = n4;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] == i1) {
	    goto L235;
	}
/* L237: */
    }
    ++mj;
    n4 += -2;
    if (n4 <= 2) {
	goto L31;
    }
    i__1 = n4;
    for (k = i__; k <= i__1; ++k) {
	ix[k - 1] = ix[k + 1];
/* L238: */
    }
    i__1 = n4;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] > i__) {
	    ix[k - 1] += -2;
	}
/* L239: */
    }
    goto L235;
L244:
    if (ix[n4 - 2] != ix[n4 - 1]) {
	goto L236;
    }
    i__1 = n4;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] == n4) {
	    goto L236;
	}
/* L246: */
    }
    ++mj;
    n4 += -2;
    if (n4 <= 2) {
	goto L31;
    }
    i__1 = n4;
    for (k = 1; k <= i__1; ++k) {
	if (ix[k - 1] == n4 + 1) {
	    ix[k - 1] = 1;
	}
/* L247: */
    }
L236:
    if (mj > 0) {
	goto L230;
    }
/* array assignments */
    m4 = n4 - 1;
    if (m4 > isi) {
	goto L1001;
    }
    i__1 = m4;
    for (ks = 1; ks <= i__1; ++ks) {
	i__2 = m4;
	for (js = 1; js <= i__2; ++js) {
	    da[ks + js * isi - (isi+1)] = 0.0;
/* L501: */
	}
	da[ks + ks * isi - (isi+1)] = 1.0;
	da[ks + (ks + 1) * isi - (isi+1)] = 1.0;
	js = ix[ks - 1];
	da[ks + js * isi - (isi+1)] = -2.0;
/* L500: */
    }
/* calculation of determinant */
    c__ = 1.0;
    kmax = m4 - 1;
    i__1 = kmax;
    for (k = 1; k <= i__1; ++k) {
	if ((r__1 = da[k + k * isi - (isi+1)], dabs(r__1)) < eps) {
	    jj = k + 1;
	    c__ = -c__;
L50:
	    if (jj > m4) {
		c__ = 0.0;
		goto L90;
	    }
	    if ((r__1 = da[jj + k * isi - (isi+1)], dabs(r__1)) > eps) {
		i__2 = m4;
		for (l = 1; l <= i__2; ++l) {
		    dr = da[k + l * isi - (isi+1)];
		    da[k + l * isi - (isi+1)] = da[jj + l * isi - (isi+1)];
/* L80: */
		    da[jj + l * isi - (isi+1)] = dr;
		}
	    } else {
		++jj;
		goto L50;
	    }
	}
	jmin = k + 1;
	i__2 = m4;
	for (j = jmin; j <= i__2; ++j) {
	    if ((r__1 = da[k + j * isi - (isi+1)], dabs(r__1)) > eps) {
		dr = da[k + j * isi - (isi+1)] / da[k + k * isi - (isi+1)];
		i__3 = m4;
		for (i__ = k; i__ <= i__3; ++i__) {
		    da[i__ + j * isi - (isi+1)] -= dr * da[i__ + k * isi - (isi+1)];
/* L22: */
		}
	    }
/* L21: */
	}
/* L20: */
    }
    i__1 = m4;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
	c__ *= da[i__ + i__ * isi - (isi+1)];
    }
L90:
    c__ = dabs(c__);
/* last assignments */
    topl = dabs(c__) + .1f;
    goto L6;
L31:
    topl = 1;
L6:
    return jwr;
L1001:
    ierr = 3;
    topl = 2;
    return 0;
L1002:
    ierr = 2;
    topl = 2;
    return 0;
} /* _kndwr_topl_update */


//Calculate integer part of writhe only.
int CircularChain::_kndwr(int &ierr){
	const int isi=801;
	const int range1=4001;

    /* System generated locals */
    int i__1, i__2, i__3;
    double r__1;
    double d__1, d__2;
    /* Local variables */
    static double c__;
    static double d__;
    static int i__, j, k, l;
    static double x[maxa], y[maxa], z[maxa];
    static int i1, n1, n2, m4, n4;
	
    static double r1, r2, da[isi*isi];
    static int n11, n21, n41, jj;
    static double dr;
    static int mj;
    static double cx[range1];
    static double dx[maxa], dy[maxa], dz[maxa];
    static int ir, ks, ix[range1], js, nv;
    static double rx, xr;
    static int ic1[range1], ic2[range1], jr1, jr2, jr3, jr4;
    static double rl1, rl2;
    static int nv1;
    static double px1, py1, px2, py2;
    static double rz1, rz2;
    static double  px11, px21, drx;
    static int jwr;
    static double zet;
    static double pdx1, pdy1, pdx2, pdy2;
    static int jmin;
    static double pdx12;
    static int kmax, iint, iiwr;

    const double deps = 1e-10;
    const double eps = 1e-4f;

	ierr=0;
/* import C[i].dX to dX*///--------------------------------
	for (int tempi=0;tempi<=maxnum;tempi++){
		x[tempi]=C[tempi].x; y[tempi]=C[tempi].y; z[tempi]=C[tempi].z;		
		dx[tempi]=C[tempi].dx; dy[tempi]=C[tempi].dy; dz[tempi]=C[tempi].dz;		
	}
	jr1=maxnum+1;
	x[maxnum+1]=C[0].x;y[maxnum+1]=C[0].y;z[maxnum+1]=C[0].z;

//_________________________________________________________

/* inventorization of intersections */
/* n4 - number of intersections */
/* cx(i) - x-value of i-th intersection */
/* ic1(i) - number of undergoing segment for i-th intersection */
/* ic2(i) - number of overgoing segment for i-th intersection */
    jr2 = jr1 - 2;
    n4 = 1;
    i__1 = jr2;
    for (n1 = 1; n1 <= i__1; ++n1) {
	n11 = n1 + 1;
	jr3 = jr1;
	if (n1 < 2) {
	    --jr3;
	}
	jr4 = n1 + 2;
	pdx1 = dx[n1 - 1];
	pdy1 = dy[n1 - 1];
	px1 = x[n1 - 1];
	py1 = y[n1 - 1];
	px11 = x[n11 - 1];
	i__2 = jr3;
	for (n2 = jr4; n2 <= i__2; ++n2) {
	    px2 = x[n2 - 1];
	    py2 = y[n2 - 1];
		if ((d__1 = px1 - px2, abs(d__1)) > 2.01 || (d__2 = py1 - py2, abs(
			d__2)) > 2.01) goto L111;
	    pdx2 = dx[n2 - 1];
	    pdy2 = dy[n2 - 1];
	    n21 = n2 + 1;
	    px21 = x[n21 - 1];
	    d__ = pdx2 * pdy1 - pdx1 * pdy2;
	    if (abs(d__) < deps) {
		goto L111;
	    }
	    pdx12 = pdx1 * pdx2;
	    drx = (py2 * pdx12 - px2 * pdx1 * pdy2 - py1 * pdx12 + px1 * pdx2 
		    * pdy1) / d__;
	    if (px1 <= px11) {
		if (drx < px1 || drx >= px11) {
		    goto L111;
		}
	    } else {
		if (drx <= px11 || drx > px1) {
		    goto L111;
		}
	    }
	    if (px2 <= px21) {
		if (drx < px2 || drx >= px21) {
		    goto L111;
		}
	    } else {
		if (drx <= px21 || drx > px2) {
		    continue;
		}
	    }
	    rx = drx;
	    if ((double) n4 > range1) {
		goto L1002;
	    }
	    rz1 = (rx - x[n1 - 1]) / dx[n1 - 1] * dz[n1 - 1] + z[n1 - 1];
	    rz2 = (rx - x[n2 - 1]) / dx[n2 - 1] * dz[n2 - 1] + z[n2 - 1];
	    cx[n4 - 1] = rx;
	    if (rz1 - rz2 >= 0.) {
		goto L401;
	    }
	    ic1[n4 - 1] = n1;
	    ic2[n4 - 1] = n2;
	    goto L402;
L401:
	    ic1[n4 - 1] = n2;
	    ic2[n4 - 1] = n1;
L402:
	    ++n4;
L111:
	    ;
	}
/* L110: */
    }
    --n4;
/* calculation of the directional writhing number */
    jwr = 0;
    i__1 = n4;
    for (iint = 1; iint <= i__1; ++iint) {
		zet = dx[ic1[iint - 1] - 1] * dy[ic2[iint - 1] - 1] - dx[ic2[iint - 1]
			 - 1] * dy[ic1[iint - 1] - 1];
		if (zet < 0.0) {
			iiwr = 1;
		} else {
			iiwr = -1;
		}
		jwr += iiwr;
	}
	ierr = 0;
	return jwr;

L1002:
    ierr = 2;
    return -9999;
} /* _kndwr */