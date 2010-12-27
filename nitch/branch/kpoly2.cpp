/* subroutine from Alex V. */
/* calculate I(-1) and I(-2) Alexander Polynomial--put at ial(1) ial(2) */
#include "chain.h"
int CircularChain::kpoly2(long ial[2],long ierr)
{
    /* Initialized data */

    static double const deps = 1.0000e-7;

    /* System generated locals */
    long i__1, i__2, i__3;
    double d__1, d__2;
	const long isi=801;
    /* Local variables */
    static double c__, d__;
    static long i__, j, k, l;
    static double t, x[maxa], y[maxa], z__[maxa];
    static long i1, i2, n1, n2, m4, n4;
    static double r1, r2, da[isi*isi]	/* was [isi][isi] */;
    static long id[isi], n11, n21, n41, jj;
    static double dr;
    static long mj;
    static double cx[isi], dx[maxa], dy[maxa], dz[maxa];
    static long ir, ks, ix[isi], js, nv;
    static double rx, xr;
    static long ic1[isi], ic2[isi], jr1, jr2, jr3, jr4, n4m;
    static double rl1, rl2;
    static long nv1;
    static double px1, py1, px2, py2, rz1, rz2;
    static long ikn;
    static double px11, px21, drx, pdx1, pdy1, pdx2, pdy2;
    static long jmin;
    static double pdx12;
    static long kmax;


	jr1=totsegnum;
	for(int iii=0;iii<=maxnum;iii++){
		x[iii] =C[iii].x;
		y[iii] =C[iii].y;
		z__[iii] =C[iii].z;
		dx[iii]=C[iii].dx;
		dy[iii]=C[iii].dy;
		dz[iii]=C[iii].dz;
	}
	x[maxnum+1]=C[0].x;
	y[maxnum+1]=C[0].y;
	z__[maxnum+1]=C[0].z;
	dx[maxnum+1]=C[0].dx;
	dy[maxnum+1]=C[0].dy;
	dz[maxnum+1]=C[0].dz;
	ierr=0;

    jr2 = jr1 - 2;
/* inventarization of intersections */
/* n4 - number of intersections */
/* cx(i) - x-value of i-th intersection */
/* ic1(i) - number of undergoing segment for i-th intersection */
/* ic2(i) - number of overgoing segment for i-th intersection */
    n4 = 1;
    t = -1.f;
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
		//This waiver is no more valid for heterogenous segmentation.
	    /*if ((d__1 = px2 - px1, fabs(d__1)) > 2.f || (d__2 = py2 - py1, fabs(
		    d__2)) > 2.f) {
		goto L111;
	    }*/
		//This waiver can be restored if the maximum segmentation is known.
	    pdx2 = dx[n2 - 1];
	    pdy2 = dy[n2 - 1];
	    n21 = n2 + 1;
	    px21 = x[n21 - 1];
	    d__ = pdx2 * pdy1 - pdx1 * pdy2;
	    if (fabs(d__) < deps) {
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
	    if (n4 > isi) {
		goto L1002;
	    }
	    rz1 = (rx - x[n1 - 1]) / dx[n1 - 1] * dz[n1 - 1] + z__[n1 - 1];
	    rz2 = (rx - x[n2 - 1]) / dx[n2 - 1] * dz[n2 - 1] + z__[n2 - 1];
	    cx[n4 - 1] = rx;
	    if (rz1 - rz2 >= 0.f) {
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
    n4m = n4;
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
	    rl1 = (d__1 = x[i1 - 1] - cx[nv1 - 1], fabs(d__1));
	    rl2 = (d__1 = x[i1 - 1] - cx[n2 - 1], fabs(d__1));
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
/* Determination of the crossing type */
L801:
    n4 = n4m;
    i__1 = n4;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i1 = ic1[i__ - 1];
	i2 = ic2[i__ - 1];
	if (dx[i1 - 1] * dy[i2 - 1] - dx[i2 - 1] * dy[i1 - 1] > 0.f) {
	    id[i__ - 1] = 1;
	} else {
	    id[i__ - 1] = -1;
	}
/* L59: */
    }
/* determination of the number ix(i) of overgoing  part of line */
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
	r1 = (d__1 = x[nv - 1] - cx[n1 - 1], fabs(d__1));
/* L74: */
	r2 = (d__1 = x[nv - 1] - cx[n2 - 1], fabs(d__1));
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
    if (t < -1.5f) {
	goto L802;
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
L802:
    m4 = n4 - 1;
    i__1 = m4;
    for (ks = 1; ks <= i__1; ++ks) {
	i__2 = m4;
	for (js = 1; js <= i__2; ++js) {
	    da[ks + js * isi - (isi+1)] = 0.;
/* L501: */
	}
	js = ix[ks - 1];
	if (js == ks || js == ks + 1) {
	    da[ks + ks * isi - (isi+1)] = -1.f;
	    da[ks + (ks + 1) * isi - (isi+1)] = 1.f;
	} else {
	    if (id[ks - 1] > 0) {
		da[ks + ks * isi - (isi+1)] = 1.f;
		da[ks + (ks + 1) * isi - (isi+1)] = -t;
	    } else {
		da[ks + ks * isi - (isi+1)] = -t;
		da[ks + (ks + 1) * isi - (isi+1)] = 1.f;
	    }
	    da[ks + js * isi - (isi+1)] = t - 1.f;
	}
/* L500: */
    }
/* calculation of determinant */
    c__ = 1.f;
    kmax = m4 - 1;
    i__1 = kmax;
    for (k = 1; k <= i__1; ++k) {
	if ((d__1 = da[k + k * isi - (isi+1)], fabs(d__1)) < deps) {
	    jj = k + 1;
	    c__ = -c__;
L50:
	    if (jj > m4) {
		c__ = 0.f;
		goto L90;
	    }
	    if ((d__1 = da[jj + k * isi - (isi+1)], fabs(d__1)) > deps) {
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
	    if ((d__1 = da[k + j * isi - (isi+1)], fabs(d__1)) > deps) {
		dr = da[k + j * isi - (isi+1)] / da[k + k * isi - (isi+1)];
		i__3 = m4;
		for (i__ = k; i__ <= i__3; ++i__) {
		    da[i__ + j * isi - (isi+1)] -= dr * da[i__ + k * isi - 
			    (isi+1)];
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
    c__ = fabs(c__);
L790:
    if (c__ > 1e7) {
	c__ /= 2.f;
	goto L790;
    }
    ikn = (long) (c__ + .1f);
    if (t == -2.f) {
L789:
	if (ikn / 2 << 1 < ikn) {
	    goto L788;
	}
	ikn /= 2;
	goto L789;
L788:
	ial[1] = ikn;
	return 0;
    }
    if (ikn == 1) {
	ial[0] = 1;
	ial[1] = 1;
	return 0;
    } else if (ikn == 3) {
	ial[0] = 3;
	ial[1] = 7;
	return 0;
    } else if (ikn == 5 || ikn == 7) {
	ial[0] = ikn;
	t = -2.f;
	goto L801;
    } else if (ikn > 8) {
	ial[0] = ikn;
	return 0;
    }
L31:
    ial[0] = 1;
    ial[1] = 1;
    return 0;
L1002:
    ierr = 2;
    return 0;
} /* kpoly_ */

