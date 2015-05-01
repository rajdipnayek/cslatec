/* sdpst.f -- translated by f2c (version 12.02.01).
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

/* DECK SDPST */
/* Subroutine */ int sdpst_(real *el, S_fp f, S_fp fa, real *h__, integer *
	impl, S_fp jacobn, integer *matdim, integer *miter, integer *ml, 
	integer *mu, integer *n, integer *nde, integer *nq, real *save2, real 
	*t, S_fp users, real *y, real *yh, real *ywt, real *uround, integer *
	nfe, integer *nje, real *a, real *dfdy, real *fac, logical *ier, 
	integer *ipvt, real *save1, integer *iswflg, real *bnd, integer *
	jstate)
{
    /* System generated locals */
    integer a_dim1, a_offset, dfdy_dim1, dfdy_offset, yh_dim1, yh_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    real r__1, r__2, r__3, r__4;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, j2;
    static real bl, bp, br, dy, yj;
    static integer mw;
    static real ys, diff;
    static integer info, imax;
    extern doublereal snrm2_(integer *, real *, integer *);
    static integer iflag;
    extern /* Subroutine */ int sgbfa_(real *, integer *, integer *, integer *
	    , integer *, integer *, integer *), sgefa_(real *, integer *, 
	    integer *, integer *, integer *);
    static real scale, facmin, factor, dfdymx;

/* ***BEGIN PROLOGUE  SDPST */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine SDPST evaluates the Jacobian matrix of the right */
/*            hand side of the differential equations. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      SINGLE PRECISION (SDPST-S, DDPST-D, CDPST-C) */
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
/* ***ROUTINES CALLED  SGBFA, SGEFA, SNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  SDPST */
/* ***FIRST EXECUTABLE STATEMENT  SDPST */
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
		*bnd = snrm2_(&i__1, &dfdy[dfdy_offset], &c__1);
	    }
	    factor = -el[*nq * 13 + 1] * *h__;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* L110: */
		    dfdy[i__ + j * dfdy_dim1] = factor * dfdy[i__ + j * 
			    dfdy_dim1];
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
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
/* Computing MAX */
		r__3 = (r__1 = ywt[j], dabs(r__1)), r__4 = (r__2 = y[j], dabs(
			r__2));
		ys = dmax(r__3,r__4);
L120:
		dy = fac[j] * ys;
		if (dy == 0.f) {
		    if (fac[j] < .5f) {
/* Computing MIN */
			r__1 = fac[j] * 100.f;
			fac[j] = dmin(r__1,.5f);
			goto L120;
		    } else {
			dy = ys;
		    }
		}
		if (*nq == 1) {
		    dy = r_sign(&dy, &save2[j]);
		} else {
		    dy = r_sign(&dy, &yh[j + yh_dim1 * 3]);
		}
		dy = y[j] + dy - y[j];
		yj = y[j];
		y[j] += dy;
		(*f)(n, t, &y[1], &save1[1]);
		if (*n == 0) {
		    *jstate = 6;
		    return 0;
		}
		y[j] = yj;
		factor = -el[*nq * 13 + 1] * *h__ / dy;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* L140: */
		    dfdy[i__ + j * dfdy_dim1] = (save1[i__] - save2[i__]) * 
			    factor;
		}
/*                                                                 Step 1 */
		diff = (r__1 = save2[1] - save1[1], dabs(r__1));
		imax = 1;
		i__1 = *n;
		for (i__ = 2; i__ <= i__1; ++i__) {
		    if ((r__1 = save2[i__] - save1[i__], dabs(r__1)) > diff) {
			imax = i__;
			diff = (r__1 = save2[i__] - save1[i__], dabs(r__1));
		    }
/* L150: */
		}
/*                                                                 Step 2 */
/* Computing MIN */
		r__3 = (r__1 = save2[imax], dabs(r__1)), r__4 = (r__2 = save1[
			imax], dabs(r__2));
		if (dmin(r__3,r__4) > 0.f) {
/* Computing MAX */
		    r__3 = (r__1 = save2[imax], dabs(r__1)), r__4 = (r__2 = 
			    save1[imax], dabs(r__2));
		    scale = dmax(r__3,r__4);
/*                                                                 Step 3 */
		    if (diff > scale * .5f) {
/* Computing MAX */
			r__1 = facmin, r__2 = fac[j] * .5f;
			fac[j] = dmax(r__1,r__2);
		    } else if (br * scale <= diff && diff <= bl * scale) {
/* Computing MIN */
			r__1 = fac[j] * 2.f;
			fac[j] = dmin(r__1,.5f);
/*                                                                 Step 4 */
		    } else if (diff < br * scale) {
/* Computing MIN */
			r__1 = bp * fac[j];
			fac[j] = dmin(r__1,.5f);
		    }
		}
/* L170: */
	    }
	    if (*iswflg == 3) {
		i__2 = *n * *n;
		*bnd = snrm2_(&i__2, &dfdy[dfdy_offset], &c__1) / (-el[*nq * 
			13 + 1] * *h__);
	    }
	    *nfe += *n;
	}
	if (*impl == 0) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L190: */
		dfdy[i__ + i__ * dfdy_dim1] += 1.f;
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
		    dfdy[i__ + j * dfdy_dim1] += a[i__ + j * a_dim1];
		}
	    }
	} else if (*impl == 2) {
	    (*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__1 = *nde;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L230: */
		dfdy[i__ + i__ * dfdy_dim1] += a[i__ + a_dim1];
	    }
	} else if (*impl == 3) {
	    (*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__1 = *nde;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *nde;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* L220: */
		    dfdy[i__ + j * dfdy_dim1] += a[i__ + j * a_dim1];
		}
	    }
	}
	sgefa_(&dfdy[dfdy_offset], matdim, n, &ipvt[1], &info);
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
	    i__2 = *n;
	    for (j = 1; j <= i__2; ++j) {
/* Computing MAX */
		i__1 = *ml + 1, i__3 = mw + 1 - j;
/* Computing MIN */
		i__5 = mw + *n - j, i__6 = mw + *ml;
		i__4 = min(i__5,i__6);
		for (i__ = max(i__1,i__3); i__ <= i__4; ++i__) {
/* L260: */
		    dfdy[i__ + j * dfdy_dim1] = factor * dfdy[i__ + j * 
			    dfdy_dim1];
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
	    i__4 = j2;
	    for (j = 1; j <= i__4; ++j) {
		i__2 = *n;
		i__1 = mw;
		for (k = j; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
/* Computing MAX */
		    r__3 = (r__1 = ywt[k], dabs(r__1)), r__4 = (r__2 = y[k], 
			    dabs(r__2));
		    ys = dmax(r__3,r__4);
L280:
		    dy = fac[k] * ys;
		    if (dy == 0.f) {
			if (fac[k] < .5f) {
/* Computing MIN */
			    r__1 = fac[k] * 100.f;
			    fac[k] = dmin(r__1,.5f);
			    goto L280;
			} else {
			    dy = ys;
			}
		    }
		    if (*nq == 1) {
			dy = r_sign(&dy, &save2[k]);
		    } else {
			dy = r_sign(&dy, &yh[k + yh_dim1 * 3]);
		    }
		    dy = y[k] + dy - y[k];
		    dfdy[mw + k * dfdy_dim1] = y[k];
/* L290: */
		    y[k] += dy;
		}
		(*f)(n, t, &y[1], &save1[1]);
		if (*n == 0) {
		    *jstate = 6;
		    return 0;
		}
		i__1 = *n;
		i__2 = mw;
		for (k = j; i__2 < 0 ? k >= i__1 : k <= i__1; k += i__2) {
		    y[k] = dfdy[mw + k * dfdy_dim1];
/* Computing MAX */
		    r__3 = (r__1 = ywt[k], dabs(r__1)), r__4 = (r__2 = y[k], 
			    dabs(r__2));
		    ys = dmax(r__3,r__4);
		    dy = fac[k] * ys;
		    if (dy == 0.f) {
			dy = ys;
		    }
		    if (*nq == 1) {
			dy = r_sign(&dy, &save2[k]);
		    } else {
			dy = r_sign(&dy, &yh[k + yh_dim1 * 3]);
		    }
		    dy = y[k] + dy - y[k];
		    factor = -el[*nq * 13 + 1] * *h__ / dy;
/* Computing MAX */
		    i__3 = *ml + 1, i__5 = mw + 1 - k;
/* Computing MIN */
		    i__7 = mw + *n - k, i__8 = mw + *ml;
		    i__6 = min(i__7,i__8);
		    for (i__ = max(i__3,i__5); i__ <= i__6; ++i__) {
/* L300: */
			dfdy[i__ + k * dfdy_dim1] = factor * (save1[i__ + k - 
				mw] - save2[i__ + k - mw]);
		    }
/*                                                                 Step 1 */
/* Computing MAX */
		    i__6 = 1, i__3 = k - *mu;
		    imax = max(i__6,i__3);
		    diff = (r__1 = save2[imax] - save1[imax], dabs(r__1));
/* Computing MAX */
		    i__6 = 1, i__3 = k - *mu;
/* Computing MIN */
		    i__7 = k + *ml;
		    i__5 = min(i__7,*n);
		    for (i__ = max(i__6,i__3) + 1; i__ <= i__5; ++i__) {
			if ((r__1 = save2[i__] - save1[i__], dabs(r__1)) > 
				diff) {
			    imax = i__;
			    diff = (r__1 = save2[i__] - save1[i__], dabs(r__1)
				    );
			}
/* L310: */
		    }
/*                                                                 Step 2 */
/* Computing MIN */
		    r__3 = (r__1 = save2[imax], dabs(r__1)), r__4 = (r__2 = 
			    save1[imax], dabs(r__2));
		    if (dmin(r__3,r__4) > 0.f) {
/* Computing MAX */
			r__3 = (r__1 = save2[imax], dabs(r__1)), r__4 = (r__2 
				= save1[imax], dabs(r__2));
			scale = dmax(r__3,r__4);
/*                                                                 Step 3 */
			if (diff > scale * .5f) {
/* Computing MAX */
			    r__1 = facmin, r__2 = fac[j] * .5f;
			    fac[j] = dmax(r__1,r__2);
			} else if (br * scale <= diff && diff <= bl * scale) {
/* Computing MIN */
			    r__1 = fac[j] * 2.f;
			    fac[j] = dmin(r__1,.5f);
/*                                                                 Step 4 */
			} else if (diff < br * scale) {
/* Computing MIN */
			    r__1 = bp * fac[k];
			    fac[k] = dmin(r__1,.5f);
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
	    i__4 = *n;
	    for (j = 1; j <= i__4; ++j) {
/* Computing MAX */
		i__2 = *ml + 1, i__1 = mw + 1 - j;
/* Computing MIN */
		i__6 = mw + *n - j, i__3 = mw + *ml;
		i__5 = min(i__6,i__3);
		for (i__ = max(i__2,i__1); i__ <= i__5; ++i__) {
/* L345: */
/* Computing MAX */
		    r__2 = dfdymx, r__3 = (r__1 = dfdy[i__ + j * dfdy_dim1], 
			    dabs(r__1));
		    dfdymx = dmax(r__2,r__3);
		}
	    }
	    *bnd = 0.f;
	    if (dfdymx != 0.f) {
		i__5 = *n;
		for (j = 1; j <= i__5; ++j) {
/* Computing MAX */
		    i__4 = *ml + 1, i__2 = mw + 1 - j;
/* Computing MIN */
		    i__6 = mw + *n - j, i__3 = mw + *ml;
		    i__1 = min(i__6,i__3);
		    for (i__ = max(i__4,i__2); i__ <= i__1; ++i__) {
/* L350: */
/* Computing 2nd power */
			r__1 = dfdy[i__ + j * dfdy_dim1] / dfdymx;
			*bnd += r__1 * r__1;
		    }
		}
		*bnd = dfdymx * sqrt(*bnd) / (-el[*nq * 13 + 1] * *h__);
	    }
	}
	if (*impl == 0) {
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* L360: */
		dfdy[mw + j * dfdy_dim1] += 1.f;
	    }
	} else if (*impl == 1) {
	    (*fa)(n, t, &y[1], &a[*ml + 1 + a_dim1], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
		i__5 = *ml + 1, i__4 = mw + 1 - j;
/* Computing MIN */
		i__6 = mw + *n - j, i__3 = mw + *ml;
		i__2 = min(i__6,i__3);
		for (i__ = max(i__5,i__4); i__ <= i__2; ++i__) {
/* L380: */
		    dfdy[i__ + j * dfdy_dim1] += a[i__ + j * a_dim1];
		}
	    }
	} else if (*impl == 2) {
	    (*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__2 = *nde;
	    for (j = 1; j <= i__2; ++j) {
/* L400: */
		dfdy[mw + j * dfdy_dim1] += a[j + a_dim1];
	    }
	} else if (*impl == 3) {
	    (*fa)(n, t, &y[1], &a[*ml + 1 + a_dim1], matdim, ml, mu, nde);
	    if (*n == 0) {
		*jstate = 9;
		return 0;
	    }
	    i__2 = *nde;
	    for (j = 1; j <= i__2; ++j) {
/* Computing MAX */
		i__1 = *ml + 1, i__5 = mw + 1 - j;
/* Computing MIN */
		i__6 = mw + *nde - j, i__3 = mw + *ml;
		i__4 = min(i__6,i__3);
		for (i__ = max(i__1,i__5); i__ <= i__4; ++i__) {
/* L390: */
		    dfdy[i__ + j * dfdy_dim1] += a[i__ + j * a_dim1];
		}
	    }
	}
	sgbfa_(&dfdy[dfdy_offset], matdim, n, ml, mu, &ipvt[1], &info);
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
} /* sdpst_ */

