/* cdntl.f -- translated by f2c (version 12.02.01).
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

static integer c__0 = 0;
static integer c__1 = 1;

/* DECK CDNTL */
/* Subroutine */ int cdntl_(real *eps, S_fp f, S_fp fa, real *hmax, real *
	hold, integer *impl, integer *jtask, integer *matdim, integer *maxord,
	 integer *mint, integer *miter, integer *ml, integer *mu, integer *n, 
	integer *nde, complex *save1, real *t, real *uround, S_fp users, 
	complex *y, complex *ywt, real *h__, integer *mntold, integer *mtrold,
	 integer *nfe, real *rc, complex *yh, complex *a, logical *convrg, 
	real *el, complex *fac, logical *ier, integer *ipvt, integer *nq, 
	integer *nwait, real *rh, real *rmax, complex *save2, real *tq, real *
	trend, integer *iswflg, integer *jstate)
{
    /* System generated locals */
    integer a_dim1, a_offset, yh_dim1, yh_offset, i__1, i__2, i__3;
    real r__1, r__2, r__3;
    complex q__1;

    /* Local variables */
    static integer i__;
    static real sum;
    static integer info;
    static real oldl0;
    extern /* Subroutine */ int cgbfa_(complex *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), cgefa_(complex *, 
	    integer *, integer *, integer *, integer *);
    static integer iflag;
    extern /* Subroutine */ int cdscl_(real *, integer *, integer *, real *, 
	    real *, real *, real *, complex *), cgbsl_(complex *, integer *, 
	    integer *, integer *, integer *, integer *, complex *, integer *),
	     cgesl_(complex *, integer *, integer *, integer *, complex *, 
	    integer *), cdcst_(integer *, integer *, integer *, real *, real *
	    );
    extern doublereal scnrm2_(integer *, complex *, integer *);

/* ***BEGIN PROLOGUE  CDNTL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine CDNTL is called to set parameters on the first */
/*            call to CDSTP, on an internal restart, or when the user has */
/*            altered MINT, MITER, and/or H. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      COMPLEX (SDNTL-S, DDNTL-D, CDNTL-C) */
/* ***AUTHOR  Kahaner, D. K., (NIST) */
/*             National Institute of Standards and Technology */
/*             Gaithersburg, MD  20899 */
/*           Sutherland, C. D., (LANL) */
/*             Mail Stop D466 */
/*             Los Alamos National Laboratory */
/*             Los Alamos, NM  87545 */
/* ***DESCRIPTION */

/*  On the first call, the order is set to 1 and the initial derivatives */
/*  are calculated.  RMAX is the maximum ratio by which H can be */
/*  increased in one step.  It is initially RMINIT to compensate */
/*  for the small initial H, but then is normally equal to RMNORM. */
/*  If a failure occurs (in corrector convergence or error test), RMAX */
/*  is set at RMFAIL for the next increase. */
/*  If the caller has changed MINT, or if JTASK = 0, CDCST is called */
/*  to set the coefficients of the method.  If the caller has changed H, */
/*  YH must be rescaled.  If H or MINT has been changed, NWAIT is */
/*  reset to NQ + 2 to prevent further increases in H for that many */
/*  steps.  Also, RC is reset.  RC is the ratio of new to old values of */
/*  the coefficient L(0)*H.  If the caller has changed MITER, RC is */
/*  set to 0 to force the partials to be updated, if partials are used. */

/* ***ROUTINES CALLED  CDCST, CDSCL, CGBFA, CGBSL, CGEFA, CGESL, SCNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  CDNTL */
/* ***FIRST EXECUTABLE STATEMENT  CDNTL */
    /* Parameter adjustments */
    a_dim1 = *matdim;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    yh_dim1 = *n;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --save1;
    --y;
    --ywt;
    el -= 14;
    --fac;
    --ipvt;
    --save2;
    tq -= 4;

    /* Function Body */
    *ier = FALSE_;
    if (*jtask >= 0) {
	if (*jtask == 0) {
	    cdcst_(maxord, mint, iswflg, &el[14], &tq[4]);
	    *rmax = 1e4f;
	}
	*rc = 0.f;
	*convrg = FALSE_;
	*trend = 1.f;
	*nq = 1;
	*nwait = 3;
	(*f)(n, t, &y[1], &save2[1]);
	if (*n == 0) {
	    *jstate = 6;
	    return 0;
	}
	++(*nfe);
	if (*impl != 0) {
	    if (*miter == 3) {
		iflag = 0;
		(*users)(&y[1], &yh[yh_offset], &ywt[1], &save1[1], &save2[1],
			 t, h__, &el[14], impl, n, nde, &iflag);
		if (iflag == -1) {
		    *ier = TRUE_;
		    return 0;
		}
		if (*n == 0) {
		    *jstate = 10;
		    return 0;
		}
	    } else if (*impl == 1) {
		if (*miter == 1 || *miter == 2) {
		    (*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
		    if (*n == 0) {
			*jstate = 9;
			return 0;
		    }
		    cgefa_(&a[a_offset], matdim, n, &ipvt[1], &info);
		    if (info != 0) {
			*ier = TRUE_;
			return 0;
		    }
		    cgesl_(&a[a_offset], matdim, n, &ipvt[1], &save2[1], &
			    c__0);
		} else if (*miter == 4 || *miter == 5) {
		    (*fa)(n, t, &y[1], &a[*ml + 1 + a_dim1], matdim, ml, mu, 
			    nde);
		    if (*n == 0) {
			*jstate = 9;
			return 0;
		    }
		    cgbfa_(&a[a_offset], matdim, n, ml, mu, &ipvt[1], &info);
		    if (info != 0) {
			*ier = TRUE_;
			return 0;
		    }
		    cgbsl_(&a[a_offset], matdim, n, ml, mu, &ipvt[1], &save2[
			    1], &c__0);
		}
	    } else if (*impl == 2) {
		(*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
		if (*n == 0) {
		    *jstate = 9;
		    return 0;
		}
		i__1 = *nde;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    i__2 = i__ + a_dim1;
		    if (a[i__2].r == 0.f && a[i__2].i == 0.f) {
			*ier = TRUE_;
			return 0;
		    } else {
			i__2 = i__;
			c_div(&q__1, &save2[i__], &a[i__ + a_dim1]);
			save2[i__2].r = q__1.r, save2[i__2].i = q__1.i;
		    }
/* L150: */
		}
		i__1 = *n;
		for (i__ = *nde + 1; i__ <= i__1; ++i__) {
/* L155: */
		    i__2 = i__ + a_dim1;
		    a[i__2].r = 0.f, a[i__2].i = 0.f;
		}
	    } else if (*impl == 3) {
		if (*miter == 1 || *miter == 2) {
		    (*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
		    if (*n == 0) {
			*jstate = 9;
			return 0;
		    }
		    cgefa_(&a[a_offset], matdim, nde, &ipvt[1], &info);
		    if (info != 0) {
			*ier = TRUE_;
			return 0;
		    }
		    cgesl_(&a[a_offset], matdim, nde, &ipvt[1], &save2[1], &
			    c__0);
		} else if (*miter == 4 || *miter == 5) {
		    (*fa)(n, t, &y[1], &a[*ml + 1 + a_dim1], matdim, ml, mu, 
			    nde);
		    if (*n == 0) {
			*jstate = 9;
			return 0;
		    }
		    cgbfa_(&a[a_offset], matdim, nde, ml, mu, &ipvt[1], &info)
			    ;
		    if (info != 0) {
			*ier = TRUE_;
			return 0;
		    }
		    cgbsl_(&a[a_offset], matdim, nde, ml, mu, &ipvt[1], &
			    save2[1], &c__0);
		}
	    }
	}
	i__2 = *nde;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L170: */
	    i__1 = i__;
	    i__3 = i__;
/* Computing MAX */
	    r__2 = 1.f, r__3 = c_abs(&ywt[i__]);
	    r__1 = dmax(r__2,r__3);
	    q__1.r = save2[i__3].r / r__1, q__1.i = save2[i__3].i / r__1;
	    save1[i__1].r = q__1.r, save1[i__1].i = q__1.i;
	}
	sum = scnrm2_(nde, &save1[1], &c__1) / sqrt((real) (*nde));
	if (sum > *eps / dabs(*h__)) {
	    r__1 = *eps / sum;
	    *h__ = r_sign(&r__1, h__);
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L180: */
	    i__3 = i__ + (yh_dim1 << 1);
	    i__2 = i__;
	    q__1.r = *h__ * save2[i__2].r, q__1.i = *h__ * save2[i__2].i;
	    yh[i__3].r = q__1.r, yh[i__3].i = q__1.i;
	}
	if (*miter == 2 || *miter == 5 || *iswflg == 3) {
	    i__3 = *n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L20: */
		i__2 = i__;
		r__1 = sqrt(*uround);
		fac[i__2].r = r__1, fac[i__2].i = 0.f;
	    }
	}
    } else {
	if (*miter != *mtrold) {
	    *mtrold = *miter;
	    *rc = 0.f;
	    *convrg = FALSE_;
	}
	if (*mint != *mntold) {
	    *mntold = *mint;
	    oldl0 = el[*nq * 13 + 1];
	    cdcst_(maxord, mint, iswflg, &el[14], &tq[4]);
	    *rc = *rc * el[*nq * 13 + 1] / oldl0;
	    *nwait = *nq + 2;
	}
	if (*h__ != *hold) {
	    *nwait = *nq + 2;
	    *rh = *h__ / *hold;
	    cdscl_(hmax, n, nq, rmax, hold, rc, rh, &yh[yh_offset]);
	}
    }
    return 0;
} /* cdntl_ */

