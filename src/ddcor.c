/* ddcor.f -- translated by f2c (version 12.02.01).
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

/* DECK DDCOR */
/* Subroutine */ int ddcor_(doublereal *dfdy, doublereal *el, S_fp fa, 
	doublereal *h__, integer *ierror, integer *impl, integer *ipvt, 
	integer *matdim, integer *miter, integer *ml, integer *mu, integer *n,
	 integer *nde, integer *nq, doublereal *t, S_fp users, doublereal *y, 
	doublereal *yh, doublereal *ywt, logical *evalfa, doublereal *save1, 
	doublereal *save2, doublereal *a, doublereal *d__, integer *jstate)
{
    /* System generated locals */
    integer a_dim1, a_offset, dfdy_dim1, dfdy_offset, yh_dim1, yh_offset, 
	    i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, mw;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer iflag;
    extern /* Subroutine */ int dgbsl_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *), dgesl_(
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *);

/* ***BEGIN PROLOGUE  DDCOR */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine DDCOR computes corrections to the Y array. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      DOUBLE PRECISION (SDCOR-S, DDCOR-D, CDCOR-C) */
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

/* ***ROUTINES CALLED  DGBSL, DGESL, DNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  DDCOR */
/* ***FIRST EXECUTABLE STATEMENT  DDCOR */
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
		save1[i__] = (*h__ * save2[i__] - yh[i__ + (yh_dim1 << 1)] - 
			save1[i__]) / ywt[i__];
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
		d__2 = (d__1 = y[i__], abs(d__1)), d__3 = ywt[i__];
		save1[i__] = (*h__ * save2[i__] - yh[i__ + (yh_dim1 << 1)] - 
			save1[i__]) / max(d__2,d__3);
/* L102: */
	    }
	}
	*d__ = dnrm2_(n, &save1[1], &c__1) / sqrt((doublereal) (*n));
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L105: */
	    save1[i__] = *h__ * save2[i__] - yh[i__ + (yh_dim1 << 1)];
	}
    } else if (*miter == 1 || *miter == 2) {
	if (*impl == 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L130: */
		save2[i__] = *h__ * save2[i__] - yh[i__ + (yh_dim1 << 1)] - 
			save1[i__];
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
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L150: */
		save2[i__] = *h__ * save2[i__];
	    }
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
		i__2 = *n;
		for (i__ = 1; i__ <= i__2; ++i__) {
/* L160: */
		    save2[i__] -= a[i__ + j * a_dim1] * (yh[j + (yh_dim1 << 1)
			    ] + save1[j]);
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
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L180: */
		save2[i__] = *h__ * save2[i__] - a[i__ + a_dim1] * (yh[i__ + (
			yh_dim1 << 1)] + save1[i__]);
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
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L140: */
		save2[i__] = *h__ * save2[i__];
	    }
	    i__2 = *nde;
	    for (j = 1; j <= i__2; ++j) {
		i__1 = *nde;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* L170: */
		    save2[i__] -= a[i__ + j * a_dim1] * (yh[j + (yh_dim1 << 1)
			    ] + save1[j]);
		}
	    }
	}
	dgesl_(&dfdy[dfdy_offset], matdim, n, &ipvt[1], &save2[1], &c__0);
	if (*ierror == 1 || *ierror == 5) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		save1[i__] += save2[i__];
/* L200: */
		save2[i__] /= ywt[i__];
	    }
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		save1[i__] += save2[i__];
/* L205: */
/* Computing MAX */
		d__2 = (d__1 = y[i__], abs(d__1)), d__3 = ywt[i__];
		save2[i__] /= max(d__2,d__3);
	    }
	}
	*d__ = dnrm2_(n, &save2[1], &c__1) / sqrt((doublereal) (*n));
    } else if (*miter == 4 || *miter == 5) {
	if (*impl == 0) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L230: */
		save2[i__] = *h__ * save2[i__] - yh[i__ + (yh_dim1 << 1)] - 
			save1[i__];
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
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L250: */
		save2[i__] = *h__ * save2[i__];
	    }
	    mw = *ml + 1 + *mu;
	    i__1 = *n;
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
		i__2 = *ml + 1, i__3 = mw + 1 - j;
/* Computing MIN */
		i__5 = mw + *n - j, i__6 = mw + *ml;
		i__4 = min(i__5,i__6);
		for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
		    save2[i__ + j - mw] -= a[i__ + j * a_dim1] * (yh[j + (
			    yh_dim1 << 1)] + save1[j]);
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
	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/* L280: */
		save2[i__] = *h__ * save2[i__] - a[i__ + a_dim1] * (yh[i__ + (
			yh_dim1 << 1)] + save1[i__]);
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
	    i__4 = *n;
	    for (i__ = 1; i__ <= i__4; ++i__) {
/* L270: */
		save2[i__] = *h__ * save2[i__];
	    }
	    mw = *ml + 1 + *mu;
	    i__4 = *nde;
	    for (j = 1; j <= i__4; ++j) {
/* Computing MAX */
		i__1 = *ml + 1, i__2 = mw + 1 - j;
/* Computing MIN */
		i__5 = mw + *nde - j, i__6 = mw + *ml;
		i__3 = min(i__5,i__6);
		for (i__ = max(i__1,i__2); i__ <= i__3; ++i__) {
		    save2[i__ + j - mw] -= a[i__ + j * a_dim1] * (yh[j + (
			    yh_dim1 << 1)] + save1[j]);
/* L290: */
		}
	    }
	}
	dgbsl_(&dfdy[dfdy_offset], matdim, n, ml, mu, &ipvt[1], &save2[1], &
		c__0);
	if (*ierror == 1 || *ierror == 5) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		save1[i__] += save2[i__];
/* L300: */
		save2[i__] /= ywt[i__];
	    }
	} else {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		save1[i__] += save2[i__];
/* L305: */
/* Computing MAX */
		d__2 = (d__1 = y[i__], abs(d__1)), d__3 = ywt[i__];
		save2[i__] /= max(d__2,d__3);
	    }
	}
	*d__ = dnrm2_(n, &save2[1], &c__1) / sqrt((doublereal) (*n));
    } else if (*miter == 3) {
	iflag = 2;
	(*users)(&y[1], &yh[(yh_dim1 << 1) + 1], &ywt[1], &save1[1], &save2[1]
		, t, h__, &el[*nq * 13 + 1], impl, n, nde, &iflag);
	if (*n == 0) {
	    *jstate = 10;
	    return 0;
	}
	if (*ierror == 1 || *ierror == 5) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		save1[i__] += save2[i__];
/* L320: */
		save2[i__] /= ywt[i__];
	    }
	} else {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		save1[i__] += save2[i__];
/* L325: */
/* Computing MAX */
		d__2 = (d__1 = y[i__], abs(d__1)), d__3 = ywt[i__];
		save2[i__] /= max(d__2,d__3);
	    }
	}
	*d__ = dnrm2_(n, &save2[1], &c__1) / sqrt((doublereal) (*n));
    }
    return 0;
} /* ddcor_ */

