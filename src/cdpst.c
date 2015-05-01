/* cdpst.f -- translated by f2c (version 12.02.01).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include <stdlib.h> /* For exit() */
#include <f2c.h>

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b4 = .875;
static doublereal c_b5 = .75;
static doublereal c_b6 = -.15;
static doublereal c_b7 = .78;

/* DECK CDPST */
/* Subroutine */ int cdpst_(real *el, S_fp f, S_fp fa, real *h__, integer *
	impl, S_fp jacobn, integer *matdim, integer *miter, integer *ml, 
	integer *mu, integer *n, integer *nde, integer *nq, complex *save2, 
	real *t, S_fp users, complex *y, complex *yh, complex *ywt, real *
	uround, integer *nfe, integer *nje, complex *a, complex *dfdy, 
	complex *fac, logical *ier, integer *ipvt, complex *save1, integer *
	iswflg, real *bnd, integer *jstate)
{
    /* System generated locals */
    integer a_dim1, a_offset, dfdy_dim1, dfdy_offset, yh_dim1, yh_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    real r__1, r__2, r__3, r__4;
    doublereal d__1;
    complex q__1, q__2;

    /* Local variables */
    static integer i__, j, k, j2;
    static real bl, bp, br;
    static complex dy, yj;
    static integer mw;
    static complex ys;
    static real diff;
    static integer info, imax;
    static real zmin, zmax;
    extern /* Subroutine */ int cgbfa_(complex *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), cgefa_(complex *, 
	    integer *, integer *, integer *, integer *);
    static integer iflag;
    static real scale;
    static complex cfctr;
    extern doublereal scnrm2_(integer *, complex *, integer *);
    static real facmin, factor, dfdymx;

/* ***BEGIN PROLOGUE  CDPST */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine CDPST evaluates the Jacobian matrix of the right */
/*            hand side of the differential equations. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      COMPLEX (SDPST-S, DDPST-D, CDPST-C) */
/* ***AUTHOR  Kahaner, D. K., (NIST) */
/*             National Institute of Standards and Technology */
/*             Gaithersburg, MD  20899 */
/*           Sutherland, C. D., (LANL) */
/*             Mail Stop D466 */
/*             Los Alamos National Laboratory */
/*             Los Alamos, NM  87545 */
/* ***DESCRIPTION */

/*  If MITER is 1, 2, 4, or 5, the matrix */
/*  P = I - L(0)*H*Jacobian is stored in DFDY and subjected to LU */
/*  decomposition, with the results also stored in DFDY. */

/* ***ROUTINES CALLED  CGBFA, CGEFA, SCNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  CDPST */
/* ***FIRST EXECUTABLE STATEMENT  CDPST */
    /* Parameter adjustments */
    el -= 14;
    dfdy_dim1 = *matdim;
    dfdy_offset = 1 + dfdy_dim1;
    dfdy -= dfdy_offset;
    a_dim1 = *matdim;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    yh_dim1 = *n;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --save2;
    --y;
    --ywt;
    --fac;
    --ipvt;
    --save1;

    /* Function Body */
    ++(*nje);
    *ier = FALSE_;
    if (*miter == 1 || *miter == 2) {
	if (*miter == 1) {
	    (*jacobn)(n, t, &y[1], &dfdy[dfdy_offset], matdim, ml, mu);
	    if (*n == 0) {
		*jstate = 8;
		return 0;
	    }
	    if (*iswflg == 3) {
		i__1 = *n * *n;
		*bnd = scnrm2_(&i__1, &dfdy[dfdy_offset], &c__1);
	    }
	    factor = -el[*nq * 13 + 1] * *h__;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* L110: */
		    i__3 = i__ + j * dfdy_dim1;
		    i__4 = i__ + j * dfdy_dim1;
		    q__1.r = factor * dfdy[i__4].r, q__1.i = factor * dfdy[
			    i__4].i;
		    dfdy[i__3].r = q__1.r, dfdy[i__3].i = q__1.i;
		}
	    }
	} else if (*miter == 2) {
	    d__1 = (doublereal) (*uround);
	    br = pow_dd(&d__1, &c_b4);
	    d__1 = (doublereal) (*uround);
	    bl = pow_dd(&d__1, &c_b5);
	    d__1 = (doublereal) (*uround);
	    bp = pow_dd(&d__1, &c_b6);
	    d__1 = (doublereal) (*uround);
	    facmin = pow_dd(&d__1, &c_b7);
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
		if (c_abs(&y[j]) > c_abs(&ywt[j])) {
		    i__4 = j;
		    ys.r = y[i__4].r, ys.i = y[i__4].i;
		} else {
		    i__4 = j;
		    ys.r = ywt[i__4].r, ys.i = ywt[i__4].i;
		}
L120:
		i__4 = j;
		q__1.r = fac[i__4].r * ys.r - fac[i__4].i * ys.i, q__1.i = 
			fac[i__4].r * ys.i + fac[i__4].i * ys.r;
		dy.r = q__1.r, dy.i = q__1.i;
		if (dy.r == 0.f && dy.i == 0.f) {
		    i__4 = j;
		    if (fac[i__4].r < .5f) {
			i__4 = j;
/* Computing MIN */
			i__2 = j;
			r__2 = fac[i__2].r * 100.f;
			r__1 = dmin(r__2,.5f);
			fac[i__4].r = r__1, fac[i__4].i = 0.f;
			goto L120;
		    } else {
			dy.r = ys.r, dy.i = ys.i;
		    }
		}
		i__4 = j;
		q__2.r = y[i__4].r + dy.r, q__2.i = y[i__4].i + dy.i;
		i__2 = j;
		q__1.r = q__2.r - y[i__2].r, q__1.i = q__2.i - y[i__2].i;
		dy.r = q__1.r, dy.i = q__1.i;
		i__4 = j;
		yj.r = y[i__4].r, yj.i = y[i__4].i;
		i__4 = j;
		i__2 = j;
		q__1.r = y[i__2].r + dy.r, q__1.i = y[i__2].i + dy.i;
		y[i__4].r = q__1.r, y[i__4].i = q__1.i;
		(*f)(n, t, &y[1], &save1[1]);
		if (*n == 0) {
		    *jstate = 6;
		    return 0;
		}
		i__4 = j;
		y[i__4].r = yj.r, y[i__4].i = yj.i;
		r__1 = -el[*nq * 13 + 1] * *h__;
		q__2.r = r__1, q__2.i = 0.f;
		c_div(&q__1, &q__2, &dy);
		cfctr.r = q__1.r, cfctr.i = q__1.i;
		i__4 = *n;
		for (i__ = 1; i__ <= i__4; ++i__) {
/* L140: */
		    i__2 = i__ + j * dfdy_dim1;
		    i__1 = i__;
		    i__5 = i__;
		    q__2.r = save1[i__1].r - save2[i__5].r, q__2.i = save1[
			    i__1].i - save2[i__5].i;
		    q__1.r = q__2.r * cfctr.r - q__2.i * cfctr.i, q__1.i = 
			    q__2.r * cfctr.i + q__2.i * cfctr.r;
		    dfdy[i__2].r = q__1.r, dfdy[i__2].i = q__1.i;
		}
/*                                                                 Step 1 */
		q__1.r = save2[1].r - save1[1].r, q__1.i = save2[1].i - save1[
			1].i;
		diff = c_abs(&q__1);
		imax = 1;
		i__2 = *n;
		for (i__ = 2; i__ <= i__2; ++i__) {
		    i__1 = i__;
		    i__5 = i__;
		    q__1.r = save2[i__1].r - save1[i__5].r, q__1.i = save2[
			    i__1].i - save1[i__5].i;
		    if (c_abs(&q__1) > diff) {
			imax = i__;
			i__1 = i__;
			i__5 = i__;
			q__1.r = save2[i__1].r - save1[i__5].r, q__1.i = 
				save2[i__1].i - save1[i__5].i;
			diff = c_abs(&q__1);
		    }
/* L150: */
		}
/*                                                                 Step 2 */
/* Computing MIN */
		r__1 = c_abs(&save2[imax]), r__2 = c_abs(&save1[imax]);
		if (dmin(r__1,r__2) > 0.f) {
/* Computing MAX */
		    r__1 = c_abs(&save2[imax]), r__2 = c_abs(&save1[imax]);
		    scale = dmax(r__1,r__2);
/*                                                                 Step 3 */
		    if (diff > scale * .5f) {
			i__2 = j;
/* Computing MAX */
			i__1 = j;
			r__2 = facmin, r__3 = fac[i__1].r * .5f;
			r__1 = dmax(r__2,r__3);
			fac[i__2].r = r__1, fac[i__2].i = 0.f;
		    } else if (br * scale <= diff && diff <= bl * scale) {
			i__2 = j;
/* Computing MIN */
			i__1 = j;
			r__2 = fac[i__1].r * 2.f;
			r__1 = dmin(r__2,.5f);
			fac[i__2].r = r__1, fac[i__2].i = 0.f;
/*                                                                 Step 4 */
		    } else if (diff < br * scale) {
			i__2 = j;
/* Computing MIN */
			i__1 = j;
			r__2 = bp * fac[i__1].r;
			r__1 = dmin(r__2,.5f);
			fac[i__2].r = r__1, fac[i__2].i = 0.f;
		    }
		}
/* L170: */
	    }
	    if (*iswflg == 3) {
		i__3 = *n * *n;
		*bnd = scnrm2_(&i__3, &dfdy[dfdy_offset], &c__1) / (-el[*nq * 
			13 + 1] * *h__);
	    }
	    *nfe += *n;
	}
	if (*impl == 0) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L190: */
		i__2 = i__ + i__ * dfdy_dim1;
		i__1 = i__ + i__ * dfdy_dim1;
		q__1.r = dfdy[i__1].r + 1.f, q__1.i = dfdy[i__1].i;
		dfdy[i__2].r = q__1.r, dfdy[i__2].i = q__1.i;
	    }
	} else if (*impl == 1) {
	    (*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* L210: */
		    i__3 = i__ + j * dfdy_dim1;
		    i__5 = i__ + j * dfdy_dim1;
		    i__4 = i__ + j * a_dim1;
		    q__1.r = dfdy[i__5].r + a[i__4].r, q__1.i = dfdy[i__5].i 
			    + a[i__4].i;
		    dfdy[i__3].r = q__1.r, dfdy[i__3].i = q__1.i;
		}
	    }
	} else if (*impl == 2) {
	    (*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__3 = *nde;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L230: */
		i__5 = i__ + i__ * dfdy_dim1;
		i__4 = i__ + i__ * dfdy_dim1;
		i__1 = i__ + a_dim1;
		q__1.r = dfdy[i__4].r + a[i__1].r, q__1.i = dfdy[i__4].i + a[
			i__1].i;
		dfdy[i__5].r = q__1.r, dfdy[i__5].i = q__1.i;
	    }
	} else if (*impl == 3) {
	    (*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__5 = *nde;
	    for (j = 1; j <= i__5; ++j) {
		i__4 = *nde;
		for (i__ = 1; i__ <= i__4; ++i__) {
/* L220: */
		    i__1 = i__ + j * dfdy_dim1;
		    i__3 = i__ + j * dfdy_dim1;
		    i__2 = i__ + j * a_dim1;
		    q__1.r = dfdy[i__3].r + a[i__2].r, q__1.i = dfdy[i__3].i 
			    + a[i__2].i;
		    dfdy[i__1].r = q__1.r, dfdy[i__1].i = q__1.i;
		}
	    }
	}
	cgefa_(&dfdy[dfdy_offset], matdim, n, &ipvt[1], &info);
	if (info != 0) {
	    *ier = TRUE_;
	}
    } else if (*miter == 4 || *miter == 5) {
	if (*miter == 4) {
	    (*jacobn)(n, t, &y[1], &dfdy[*ml + 1 + dfdy_dim1], matdim, ml, mu)
		    ;
	    if (*n == 0) {
		*jstate = 8;
		return 0;
	    }
	    factor = -el[*nq * 13 + 1] * *h__;
	    mw = *ml + *mu + 1;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
		i__3 = *ml + 1, i__2 = mw + 1 - j;
/* Computing MIN */
		i__5 = mw + *n - j, i__6 = mw + *ml;
		i__4 = min(i__5,i__6);
		for (i__ = max(i__3,i__2); i__ <= i__4; ++i__) {
/* L260: */
		    i__3 = i__ + j * dfdy_dim1;
		    i__2 = i__ + j * dfdy_dim1;
		    q__1.r = factor * dfdy[i__2].r, q__1.i = factor * dfdy[
			    i__2].i;
		    dfdy[i__3].r = q__1.r, dfdy[i__3].i = q__1.i;
		}
	    }
	} else if (*miter == 5) {
	    d__1 = (doublereal) (*uround);
	    br = pow_dd(&d__1, &c_b4);
	    d__1 = (doublereal) (*uround);
	    bl = pow_dd(&d__1, &c_b5);
	    d__1 = (doublereal) (*uround);
	    bp = pow_dd(&d__1, &c_b6);
	    d__1 = (doublereal) (*uround);
	    facmin = pow_dd(&d__1, &c_b7);
	    mw = *ml + *mu + 1;
	    j2 = min(mw,*n);
	    i__3 = j2;
	    for (j = 1; j <= i__3; ++j) {
		i__2 = *n;
		i__4 = mw;
		for (k = j; i__4 < 0 ? k >= i__2 : k <= i__2; k += i__4) {
		    if (c_abs(&y[k]) > c_abs(&ywt[k])) {
			i__1 = k;
			ys.r = y[i__1].r, ys.i = y[i__1].i;
		    } else {
			i__1 = k;
			ys.r = ywt[i__1].r, ys.i = ywt[i__1].i;
		    }
L280:
		    i__1 = k;
		    q__1.r = fac[i__1].r * ys.r - fac[i__1].i * ys.i, q__1.i =
			     fac[i__1].r * ys.i + fac[i__1].i * ys.r;
		    dy.r = q__1.r, dy.i = q__1.i;
		    if (dy.r == 0.f && dy.i == 0.f) {
			i__1 = k;
			if (fac[i__1].r < .5f) {
			    i__1 = k;
/* Computing MIN */
			    i__5 = k;
			    r__2 = fac[i__5].r * 100.f;
			    r__1 = dmin(r__2,.5f);
			    fac[i__1].r = r__1, fac[i__1].i = 0.f;
			    goto L280;
			} else {
			    dy.r = ys.r, dy.i = ys.i;
			}
		    }
		    i__1 = k;
		    q__2.r = y[i__1].r + dy.r, q__2.i = y[i__1].i + dy.i;
		    i__5 = k;
		    q__1.r = q__2.r - y[i__5].r, q__1.i = q__2.i - y[i__5].i;
		    dy.r = q__1.r, dy.i = q__1.i;
		    i__1 = mw + k * dfdy_dim1;
		    i__5 = k;
		    dfdy[i__1].r = y[i__5].r, dfdy[i__1].i = y[i__5].i;
/* L290: */
		    i__1 = k;
		    i__5 = k;
		    q__1.r = y[i__5].r + dy.r, q__1.i = y[i__5].i + dy.i;
		    y[i__1].r = q__1.r, y[i__1].i = q__1.i;
		}
		(*f)(n, t, &y[1], &save1[1]);
		if (*n == 0) {
		    *jstate = 6;
		    return 0;
		}
		i__1 = *n;
		i__5 = mw;
		for (k = j; i__5 < 0 ? k >= i__1 : k <= i__1; k += i__5) {
		    i__4 = k;
		    i__2 = mw + k * dfdy_dim1;
		    q__1.r = y[i__4].r - dfdy[i__2].r, q__1.i = y[i__4].i - 
			    dfdy[i__2].i;
		    dy.r = q__1.r, dy.i = q__1.i;
		    i__4 = k;
		    i__2 = mw + k * dfdy_dim1;
		    y[i__4].r = dfdy[i__2].r, y[i__4].i = dfdy[i__2].i;
		    r__1 = -el[*nq * 13 + 1] * *h__;
		    q__2.r = r__1, q__2.i = 0.f;
		    c_div(&q__1, &q__2, &dy);
		    cfctr.r = q__1.r, cfctr.i = q__1.i;
/* Computing MAX */
		    i__4 = *ml + 1, i__2 = mw + 1 - k;
/* Computing MIN */
		    i__7 = mw + *n - k, i__8 = mw + *ml;
		    i__6 = min(i__7,i__8);
		    for (i__ = max(i__4,i__2); i__ <= i__6; ++i__) {
/* L300: */
			i__4 = i__ + k * dfdy_dim1;
			i__2 = i__ + k - mw;
			i__7 = i__ + k - mw;
			q__2.r = save1[i__2].r - save2[i__7].r, q__2.i = 
				save1[i__2].i - save2[i__7].i;
			q__1.r = cfctr.r * q__2.r - cfctr.i * q__2.i, q__1.i =
				 cfctr.r * q__2.i + cfctr.i * q__2.r;
			dfdy[i__4].r = q__1.r, dfdy[i__4].i = q__1.i;
		    }
/*                                                                 Step 1 */
/* Computing MAX */
		    i__4 = 1, i__2 = k - *mu;
		    imax = max(i__4,i__2);
		    i__4 = imax;
		    i__2 = imax;
		    q__1.r = save2[i__4].r - save1[i__2].r, q__1.i = save2[
			    i__4].i - save1[i__2].i;
		    diff = c_abs(&q__1);
/* Computing MAX */
		    i__4 = 1, i__2 = k - *mu;
/* Computing MIN */
		    i__6 = k + *ml;
		    i__7 = min(i__6,*n);
		    for (i__ = max(i__4,i__2) + 1; i__ <= i__7; ++i__) {
			i__4 = i__;
			i__2 = i__;
			q__1.r = save2[i__4].r - save1[i__2].r, q__1.i = 
				save2[i__4].i - save1[i__2].i;
			if (c_abs(&q__1) > diff) {
			    imax = i__;
			    i__4 = i__;
			    i__2 = i__;
			    q__1.r = save2[i__4].r - save1[i__2].r, q__1.i = 
				    save2[i__4].i - save1[i__2].i;
			    diff = c_abs(&q__1);
			}
/* L310: */
		    }
/*                                                                 Step 2 */
/* Computing MIN */
		    r__1 = c_abs(&save2[imax]), r__2 = c_abs(&save1[imax]);
		    if (dmin(r__1,r__2) > 0.f) {
/* Computing MAX */
			r__1 = c_abs(&save2[imax]), r__2 = c_abs(&save1[imax])
				;
			scale = dmax(r__1,r__2);
/*                                                                 Step 3 */
			if (diff > scale * .5f) {
			    i__7 = j;
/* Computing MAX */
			    i__4 = j;
			    r__2 = facmin, r__3 = fac[i__4].r * .5f;
			    r__1 = dmax(r__2,r__3);
			    fac[i__7].r = r__1, fac[i__7].i = 0.f;
			} else if (br * scale <= diff && diff <= bl * scale) {
			    i__7 = j;
/* Computing MIN */
			    i__4 = j;
			    r__2 = fac[i__4].r * 2.f;
			    r__1 = dmin(r__2,.5f);
			    fac[i__7].r = r__1, fac[i__7].i = 0.f;
/*                                                                 Step 4 */
			} else if (diff < br * scale) {
			    i__7 = k;
/* Computing MIN */
			    i__4 = k;
			    r__2 = bp * fac[i__4].r;
			    r__1 = dmin(r__2,.5f);
			    fac[i__7].r = r__1, fac[i__7].i = 0.f;
			}
		    }
/* L330: */
		}
/* L340: */
	    }
	    *nfe += j2;
	}
	if (*iswflg == 3) {
	    dfdymx = 0.f;
	    i__3 = *n;
	    for (j = 1; j <= i__3; ++j) {
/* Computing MAX */
		i__5 = *ml + 1, i__1 = mw + 1 - j;
/* Computing MIN */
		i__4 = mw + *n - j, i__2 = mw + *ml;
		i__7 = min(i__4,i__2);
		for (i__ = max(i__5,i__1); i__ <= i__7; ++i__) {
/* Computing MAX */
		    i__5 = i__ + j * dfdy_dim1;
		    r__3 = (r__1 = dfdy[i__5].r, dabs(r__1)), r__4 = (r__2 = 
			    r_imag(&dfdy[i__ + j * dfdy_dim1]), dabs(r__2));
		    zmax = dmax(r__3,r__4);
/* Computing MIN */
		    i__5 = i__ + j * dfdy_dim1;
		    r__3 = (r__1 = dfdy[i__5].r, dabs(r__1)), r__4 = (r__2 = 
			    r_imag(&dfdy[i__ + j * dfdy_dim1]), dabs(r__2));
		    zmin = dmin(r__3,r__4);
		    if (zmax != 0.f) {
/* Computing MAX */
/* Computing 2nd power */
			r__3 = zmin / zmax;
			r__1 = dfdymx, r__2 = zmax * sqrt(r__3 * r__3 + 1.f);
			dfdymx = dmax(r__1,r__2);
		    }
/* L345: */
		}
	    }
	    *bnd = 0.f;
	    if (dfdymx != 0.f) {
		i__7 = *n;
		for (j = 1; j <= i__7; ++j) {
/* Computing MAX */
		    i__3 = *ml + 1, i__5 = mw + 1 - j;
/* Computing MIN */
		    i__4 = mw + *n - j, i__2 = mw + *ml;
		    i__1 = min(i__4,i__2);
		    for (i__ = max(i__3,i__5); i__ <= i__1; ++i__) {
			i__3 = i__ + j * dfdy_dim1;
/* Computing 2nd power */
			r__1 = dfdy[i__3].r / dfdymx;
/* Computing 2nd power */
			r__2 = r_imag(&dfdy[i__ + j * dfdy_dim1]) / dfdymx;
			*bnd = *bnd + r__1 * r__1 + r__2 * r__2;
/* L350: */
		    }
		}
		*bnd = dfdymx * sqrt(*bnd) / (-el[*nq * 13 + 1] * *h__);
	    }
	}
	if (*impl == 0) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* L360: */
		i__7 = mw + j * dfdy_dim1;
		i__3 = mw + j * dfdy_dim1;
		q__1.r = dfdy[i__3].r + 1.f, q__1.i = dfdy[i__3].i;
		dfdy[i__7].r = q__1.r, dfdy[i__7].i = q__1.i;
	    }
	} else if (*impl == 1) {
	    (*fa)(n, t, &y[1], &a[*ml + 1 + a_dim1], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__7 = *n;
	    for (j = 1; j <= i__7; ++j) {
/* Computing MAX */
		i__3 = *ml + 1, i__1 = mw + 1 - j;
/* Computing MIN */
		i__4 = mw + *n - j, i__2 = mw + *ml;
		i__5 = min(i__4,i__2);
		for (i__ = max(i__3,i__1); i__ <= i__5; ++i__) {
/* L380: */
		    i__3 = i__ + j * dfdy_dim1;
		    i__1 = i__ + j * dfdy_dim1;
		    i__4 = i__ + j * a_dim1;
		    q__1.r = dfdy[i__1].r + a[i__4].r, q__1.i = dfdy[i__1].i 
			    + a[i__4].i;
		    dfdy[i__3].r = q__1.r, dfdy[i__3].i = q__1.i;
		}
	    }
	} else if (*impl == 2) {
	    (*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__3 = *nde;
	    for (j = 1; j <= i__3; ++j) {
/* L400: */
		i__1 = mw + j * dfdy_dim1;
		i__4 = mw + j * dfdy_dim1;
		i__5 = j + a_dim1;
		q__1.r = dfdy[i__4].r + a[i__5].r, q__1.i = dfdy[i__4].i + a[
			i__5].i;
		dfdy[i__1].r = q__1.r, dfdy[i__1].i = q__1.i;
	    }
	} else if (*impl == 3) {
	    (*fa)(n, t, &y[1], &a[*ml + 1 + a_dim1], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__1 = *nde;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
		i__4 = *ml + 1, i__5 = mw + 1 - j;
/* Computing MIN */
		i__7 = mw + *nde - j, i__2 = mw + *ml;
		i__3 = min(i__7,i__2);
		for (i__ = max(i__4,i__5); i__ <= i__3; ++i__) {
/* L390: */
		    i__4 = i__ + j * dfdy_dim1;
		    i__5 = i__ + j * dfdy_dim1;
		    i__7 = i__ + j * a_dim1;
		    q__1.r = dfdy[i__5].r + a[i__7].r, q__1.i = dfdy[i__5].i 
			    + a[i__7].i;
		    dfdy[i__4].r = q__1.r, dfdy[i__4].i = q__1.i;
		}
	    }
	}
	cgbfa_(&dfdy[dfdy_offset], matdim, n, ml, mu, &ipvt[1], &info);
	if (info != 0) {
	    *ier = TRUE_;
	}
    } else if (*miter == 3) {
	iflag = 1;
	(*users)(&y[1], &yh[(yh_dim1 << 1) + 1], &ywt[1], &save1[1], &save2[1]
		, t, h__, &el[*nq * 13 + 1], impl, n, nde, &iflag);
	if (iflag == -1) {
	    *ier = TRUE_;
	    return 0;
	}
	if (*n == 0) {
	    *jstate = 10;
	    return 0;
	}
    }
    return 0;
} /* cdpst_ */

