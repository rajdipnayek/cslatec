/* cdcor.f -- translated by f2c (version 12.02.01).
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
static integer c__0 = 0;

/* DECK CDCOR */
/* Subroutine */ int cdcor_(complex *dfdy, real *el, S_fp fa, real *h__, 
	integer *ierror, integer *impl, integer *ipvt, integer *matdim, 
	integer *miter, integer *ml, integer *mu, integer *n, integer *nde, 
	integer *nq, real *t, S_fp users, complex *y, complex *yh, complex *
	ywt, logical *evalfa, complex *save1, complex *save2, complex *a, 
	real *d__, integer *jstate)
{
    /* System generated locals */
    integer a_dim1, a_offset, dfdy_dim1, dfdy_offset, yh_dim1, yh_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    real r__1, r__2, r__3;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    static integer i__, j, mw, iflag;
    extern /* Subroutine */ int cgbsl_(complex *, integer *, integer *, 
	    integer *, integer *, integer *, complex *, integer *), cgesl_(
	    complex *, integer *, integer *, integer *, complex *, integer *);
    extern doublereal scnrm2_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CDCOR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine CDCOR computes corrections to the Y array. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      COMPLEX (SDCOR-S, DDCOR-D, CDCOR-C) */
/* ***AUTHOR  Kahaner, D. K., (NIST) */
/*             National Institute of Standards and Technology */
/*             Gaithersburg, MD  20899 */
/*           Sutherland, C. D., (LANL) */
/*             Mail Stop D466 */
/*             Los Alamos National Laboratory */
/*             Los Alamos, NM  87545 */
/* ***DESCRIPTION */

/*  In the case of functional iteration, update Y directly from the */
/*  result of the last call to F. */
/*  In the case of the chord method, compute the corrector error and */
/*  solve the linear system with that as right hand side and DFDY as */
/*  coefficient matrix, using the LU decomposition if MITER is 1, 2, 4, */
/*  or 5. */

/* ***ROUTINES CALLED  CGBSL, CGESL, SCNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  CDCOR */
/* ***FIRST EXECUTABLE STATEMENT  CDCOR */
    /* Parameter adjustments */
    el -= 14;
    --ipvt;
    a_dim1 = *matdim;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    dfdy_dim1 = *matdim;
    dfdy_offset = 1 + dfdy_dim1;
    dfdy -= dfdy_offset;
    yh_dim1 = *n;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --y;
    --ywt;
    --save1;
    --save2;

    /* Function Body */
    if (*miter == 0) {
	if (*ierror == 1 || *ierror == 5) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L100: */
		i__2 = i__;
		i__3 = i__;
		q__4.r = *h__ * save2[i__3].r, q__4.i = *h__ * save2[i__3].i;
		i__4 = i__ + (yh_dim1 << 1);
		q__3.r = q__4.r - yh[i__4].r, q__3.i = q__4.i - yh[i__4].i;
		i__5 = i__;
		q__2.r = q__3.r - save1[i__5].r, q__2.i = q__3.i - save1[i__5]
			.i;
		c_div(&q__1, &q__2, &ywt[i__]);
		save1[i__2].r = q__1.r, save1[i__2].i = q__1.i;
	    }
	} else {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__3 = i__;
		i__4 = i__;
		q__4.r = *h__ * save2[i__4].r, q__4.i = *h__ * save2[i__4].i;
		i__5 = i__ + (yh_dim1 << 1);
		q__3.r = q__4.r - yh[i__5].r, q__3.i = q__4.i - yh[i__5].i;
		i__1 = i__;
		q__2.r = q__3.r - save1[i__1].r, q__2.i = q__3.i - save1[i__1]
			.i;
/* Computing MAX */
		r__2 = c_abs(&y[i__]), r__3 = c_abs(&ywt[i__]);
		r__1 = dmax(r__2,r__3);
		q__1.r = q__2.r / r__1, q__1.i = q__2.i / r__1;
		save1[i__3].r = q__1.r, save1[i__3].i = q__1.i;
/* L102: */
	    }
	}
	*d__ = scnrm2_(n, &save1[1], &c__1) / sqrt((real) (*n));
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L105: */
	    i__3 = i__;
	    i__4 = i__;
	    q__2.r = *h__ * save2[i__4].r, q__2.i = *h__ * save2[i__4].i;
	    i__5 = i__ + (yh_dim1 << 1);
	    q__1.r = q__2.r - yh[i__5].r, q__1.i = q__2.i - yh[i__5].i;
	    save1[i__3].r = q__1.r, save1[i__3].i = q__1.i;
	}
    } else if (*miter == 1 || *miter == 2) {
	if (*impl == 0) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L130: */
		i__4 = i__;
		i__5 = i__;
		q__3.r = *h__ * save2[i__5].r, q__3.i = *h__ * save2[i__5].i;
		i__2 = i__ + (yh_dim1 << 1);
		q__2.r = q__3.r - yh[i__2].r, q__2.i = q__3.i - yh[i__2].i;
		i__1 = i__;
		q__1.r = q__2.r - save1[i__1].r, q__1.i = q__2.i - save1[i__1]
			.i;
		save2[i__4].r = q__1.r, save2[i__4].i = q__1.i;
	    }
	} else if (*impl == 1) {
	    if (*evalfa) {
		(*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
		if (*n == 0) {
		    *jstate = 9;
		    return 0;
		}
	    } else {
		*evalfa = TRUE_;
	    }
	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/* L150: */
		i__5 = i__;
		i__2 = i__;
		q__1.r = *h__ * save2[i__2].r, q__1.i = *h__ * save2[i__2].i;
		save2[i__5].r = q__1.r, save2[i__5].i = q__1.i;
	    }
	    i__5 = *n;
	    for (j = 1; j <= i__5; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* L160: */
		    i__4 = i__;
		    i__1 = i__;
		    i__3 = i__ + j * a_dim1;
		    i__6 = j + (yh_dim1 << 1);
		    i__7 = j;
		    q__3.r = yh[i__6].r + save1[i__7].r, q__3.i = yh[i__6].i 
			    + save1[i__7].i;
		    q__2.r = a[i__3].r * q__3.r - a[i__3].i * q__3.i, q__2.i =
			     a[i__3].r * q__3.i + a[i__3].i * q__3.r;
		    q__1.r = save2[i__1].r - q__2.r, q__1.i = save2[i__1].i - 
			    q__2.i;
		    save2[i__4].r = q__1.r, save2[i__4].i = q__1.i;
		}
	    }
	} else if (*impl == 2) {
	    if (*evalfa) {
		(*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
		if (*n == 0) {
		    *jstate = 9;
		    return 0;
		}
	    } else {
		*evalfa = TRUE_;
	    }
	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/* L180: */
		i__1 = i__;
		i__3 = i__;
		q__2.r = *h__ * save2[i__3].r, q__2.i = *h__ * save2[i__3].i;
		i__6 = i__ + a_dim1;
		i__7 = i__ + (yh_dim1 << 1);
		i__2 = i__;
		q__4.r = yh[i__7].r + save1[i__2].r, q__4.i = yh[i__7].i + 
			save1[i__2].i;
		q__3.r = a[i__6].r * q__4.r - a[i__6].i * q__4.i, q__3.i = a[
			i__6].r * q__4.i + a[i__6].i * q__4.r;
		q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
		save2[i__1].r = q__1.r, save2[i__1].i = q__1.i;
	    }
	} else if (*impl == 3) {
	    if (*evalfa) {
		(*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
		if (*n == 0) {
		    *jstate = 9;
		    return 0;
		}
	    } else {
		*evalfa = TRUE_;
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L140: */
		i__3 = i__;
		i__6 = i__;
		q__1.r = *h__ * save2[i__6].r, q__1.i = *h__ * save2[i__6].i;
		save2[i__3].r = q__1.r, save2[i__3].i = q__1.i;
	    }
	    i__3 = *nde;
	    for (j = 1; j <= i__3; ++j) {
		i__6 = *nde;
		for (i__ = 1; i__ <= i__6; ++i__) {
/* L170: */
		    i__1 = i__;
		    i__7 = i__;
		    i__2 = i__ + j * a_dim1;
		    i__4 = j + (yh_dim1 << 1);
		    i__5 = j;
		    q__3.r = yh[i__4].r + save1[i__5].r, q__3.i = yh[i__4].i 
			    + save1[i__5].i;
		    q__2.r = a[i__2].r * q__3.r - a[i__2].i * q__3.i, q__2.i =
			     a[i__2].r * q__3.i + a[i__2].i * q__3.r;
		    q__1.r = save2[i__7].r - q__2.r, q__1.i = save2[i__7].i - 
			    q__2.i;
		    save2[i__1].r = q__1.r, save2[i__1].i = q__1.i;
		}
	    }
	}
	cgesl_(&dfdy[dfdy_offset], matdim, n, &ipvt[1], &save2[1], &c__0);
	if (*ierror == 1 || *ierror == 5) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__7 = i__;
		i__2 = i__;
		i__4 = i__;
		q__1.r = save1[i__2].r + save2[i__4].r, q__1.i = save1[i__2]
			.i + save2[i__4].i;
		save1[i__7].r = q__1.r, save1[i__7].i = q__1.i;
/* L200: */
		i__7 = i__;
		c_div(&q__1, &save2[i__], &ywt[i__]);
		save2[i__7].r = q__1.r, save2[i__7].i = q__1.i;
	    }
	} else {
	    i__7 = *n;
	    for (i__ = 1; i__ <= i__7; ++i__) {
		i__1 = i__;
		i__2 = i__;
		i__4 = i__;
		q__1.r = save1[i__2].r + save2[i__4].r, q__1.i = save1[i__2]
			.i + save2[i__4].i;
		save1[i__1].r = q__1.r, save1[i__1].i = q__1.i;
/* L205: */
		i__1 = i__;
		i__2 = i__;
/* Computing MAX */
		r__2 = c_abs(&y[i__]), r__3 = c_abs(&ywt[i__]);
		r__1 = dmax(r__2,r__3);
		q__1.r = save2[i__2].r / r__1, q__1.i = save2[i__2].i / r__1;
		save2[i__1].r = q__1.r, save2[i__1].i = q__1.i;
	    }
	}
	*d__ = scnrm2_(n, &save2[1], &c__1) / sqrt((real) (*n));
    } else if (*miter == 4 || *miter == 5) {
	if (*impl == 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L230: */
		i__2 = i__;
		i__7 = i__;
		q__3.r = *h__ * save2[i__7].r, q__3.i = *h__ * save2[i__7].i;
		i__4 = i__ + (yh_dim1 << 1);
		q__2.r = q__3.r - yh[i__4].r, q__2.i = q__3.i - yh[i__4].i;
		i__5 = i__;
		q__1.r = q__2.r - save1[i__5].r, q__1.i = q__2.i - save1[i__5]
			.i;
		save2[i__2].r = q__1.r, save2[i__2].i = q__1.i;
	    }
	} else if (*impl == 1) {
	    if (*evalfa) {
		(*fa)(n, t, &y[1], &a[*ml + 1 + a_dim1], matdim, ml, mu, nde);
		if (*n == 0) {
		    *jstate = 9;
		    return 0;
		}
	    } else {
		*evalfa = TRUE_;
	    }
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L250: */
		i__7 = i__;
		i__4 = i__;
		q__1.r = *h__ * save2[i__4].r, q__1.i = *h__ * save2[i__4].i;
		save2[i__7].r = q__1.r, save2[i__7].i = q__1.i;
	    }
	    mw = *ml + 1 + *mu;
	    i__7 = *n;
	    for (j = 1; j <= i__7; ++j) {
/* Computing MAX */
		i__4 = *ml + 1, i__2 = mw + 1 - j;
/* Computing MIN */
		i__1 = mw + *n - j, i__6 = mw + *ml;
		i__5 = min(i__1,i__6);
		for (i__ = max(i__4,i__2); i__ <= i__5; ++i__) {
		    i__4 = i__ + j - mw;
		    i__2 = i__ + j - mw;
		    i__1 = i__ + j * a_dim1;
		    i__6 = j + (yh_dim1 << 1);
		    i__3 = j;
		    q__3.r = yh[i__6].r + save1[i__3].r, q__3.i = yh[i__6].i 
			    + save1[i__3].i;
		    q__2.r = a[i__1].r * q__3.r - a[i__1].i * q__3.i, q__2.i =
			     a[i__1].r * q__3.i + a[i__1].i * q__3.r;
		    q__1.r = save2[i__2].r - q__2.r, q__1.i = save2[i__2].i - 
			    q__2.i;
		    save2[i__4].r = q__1.r, save2[i__4].i = q__1.i;
/* L260: */
		}
	    }
	} else if (*impl == 2) {
	    if (*evalfa) {
		(*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
		if (*n == 0) {
		    *jstate = 9;
		    return 0;
		}
	    } else {
		*evalfa = TRUE_;
	    }
	    i__5 = *n;
	    for (i__ = 1; i__ <= i__5; ++i__) {
/* L280: */
		i__7 = i__;
		i__4 = i__;
		q__2.r = *h__ * save2[i__4].r, q__2.i = *h__ * save2[i__4].i;
		i__2 = i__ + a_dim1;
		i__1 = i__ + (yh_dim1 << 1);
		i__6 = i__;
		q__4.r = yh[i__1].r + save1[i__6].r, q__4.i = yh[i__1].i + 
			save1[i__6].i;
		q__3.r = a[i__2].r * q__4.r - a[i__2].i * q__4.i, q__3.i = a[
			i__2].r * q__4.i + a[i__2].i * q__4.r;
		q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
		save2[i__7].r = q__1.r, save2[i__7].i = q__1.i;
	    }
	} else if (*impl == 3) {
	    if (*evalfa) {
		(*fa)(n, t, &y[1], &a[*ml + 1 + a_dim1], matdim, ml, mu, nde);
		if (*n == 0) {
		    *jstate = 9;
		    return 0;
		}
	    } else {
		*evalfa = TRUE_;
	    }
	    i__7 = *n;
	    for (i__ = 1; i__ <= i__7; ++i__) {
/* L270: */
		i__4 = i__;
		i__2 = i__;
		q__1.r = *h__ * save2[i__2].r, q__1.i = *h__ * save2[i__2].i;
		save2[i__4].r = q__1.r, save2[i__4].i = q__1.i;
	    }
	    mw = *ml + 1 + *mu;
	    i__4 = *nde;
	    for (j = 1; j <= i__4; ++j) {
/* Computing MAX */
		i__2 = *ml + 1, i__7 = mw + 1 - j;
/* Computing MIN */
		i__6 = mw + *nde - j, i__5 = mw + *ml;
		i__1 = min(i__6,i__5);
		for (i__ = max(i__2,i__7); i__ <= i__1; ++i__) {
		    i__2 = i__ + j - mw;
		    i__7 = i__ + j - mw;
		    i__6 = i__ + j * a_dim1;
		    i__5 = j + (yh_dim1 << 1);
		    i__3 = j;
		    q__3.r = yh[i__5].r + save1[i__3].r, q__3.i = yh[i__5].i 
			    + save1[i__3].i;
		    q__2.r = a[i__6].r * q__3.r - a[i__6].i * q__3.i, q__2.i =
			     a[i__6].r * q__3.i + a[i__6].i * q__3.r;
		    q__1.r = save2[i__7].r - q__2.r, q__1.i = save2[i__7].i - 
			    q__2.i;
		    save2[i__2].r = q__1.r, save2[i__2].i = q__1.i;
/* L290: */
		}
	    }
	}
	cgbsl_(&dfdy[dfdy_offset], matdim, n, ml, mu, &ipvt[1], &save2[1], &
		c__0);
	if (*ierror == 1 || *ierror == 5) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__4 = i__;
		i__2 = i__;
		i__7 = i__;
		q__1.r = save1[i__2].r + save2[i__7].r, q__1.i = save1[i__2]
			.i + save2[i__7].i;
		save1[i__4].r = q__1.r, save1[i__4].i = q__1.i;
/* L300: */
		i__4 = i__;
		c_div(&q__1, &save2[i__], &ywt[i__]);
		save2[i__4].r = q__1.r, save2[i__4].i = q__1.i;
	    }
	} else {
	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
		i__1 = i__;
		i__2 = i__;
		i__7 = i__;
		q__1.r = save1[i__2].r + save2[i__7].r, q__1.i = save1[i__2]
			.i + save2[i__7].i;
		save1[i__1].r = q__1.r, save1[i__1].i = q__1.i;
/* L305: */
		i__1 = i__;
		i__2 = i__;
/* Computing MAX */
		r__2 = c_abs(&y[i__]), r__3 = c_abs(&ywt[i__]);
		r__1 = dmax(r__2,r__3);
		q__1.r = save2[i__2].r / r__1, q__1.i = save2[i__2].i / r__1;
		save2[i__1].r = q__1.r, save2[i__1].i = q__1.i;
	    }
	}
	*d__ = scnrm2_(n, &save2[1], &c__1) / sqrt((real) (*n));
    } else if (*miter == 3) {
	iflag = 2;
	(*users)(&y[1], &yh[(yh_dim1 << 1) + 1], &ywt[1], &save1[1], &save2[1]
		, t, h__, &el[*nq * 13 + 1], impl, n, nde, &iflag);
	if (*n == 0) {
	    *jstate = 10;
	    return 0;
	}
	if (*ierror == 1 || *ierror == 5) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = i__;
		i__4 = i__;
		i__7 = i__;
		q__1.r = save1[i__4].r + save2[i__7].r, q__1.i = save1[i__4]
			.i + save2[i__7].i;
		save1[i__2].r = q__1.r, save1[i__2].i = q__1.i;
/* L320: */
		i__2 = i__;
		c_div(&q__1, &save2[i__], &ywt[i__]);
		save2[i__2].r = q__1.r, save2[i__2].i = q__1.i;
	    }
	} else {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		i__1 = i__;
		i__4 = i__;
		i__7 = i__;
		q__1.r = save1[i__4].r + save2[i__7].r, q__1.i = save1[i__4]
			.i + save2[i__7].i;
		save1[i__1].r = q__1.r, save1[i__1].i = q__1.i;
/* L325: */
		i__1 = i__;
		i__4 = i__;
/* Computing MAX */
		r__2 = c_abs(&y[i__]), r__3 = c_abs(&ywt[i__]);
		r__1 = dmax(r__2,r__3);
		q__1.r = save2[i__4].r / r__1, q__1.i = save2[i__4].i / r__1;
		save2[i__1].r = q__1.r, save2[i__1].i = q__1.i;
	    }
	}
	*d__ = scnrm2_(n, &save2[1], &c__1) / sqrt((real) (*n));
    }
    return 0;
} /* cdcor_ */

