/* sdntl.f -- translated by f2c (version 12.02.01).
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

/* DECK SDNTL */
/* Subroutine */ int sdntl_(real *eps, S_fp f, S_fp fa, real *hmax, real *
	hold, integer *impl, integer *jtask, integer *matdim, integer *maxord,
	 integer *mint, integer *miter, integer *ml, integer *mu, integer *n, 
	integer *nde, real *save1, real *t, real *uround, S_fp users, real *y,
	 real *ywt, real *h__, integer *mntold, integer *mtrold, integer *nfe,
	 real *rc, real *yh, real *a, logical *convrg, real *el, real *fac, 
	logical *ier, integer *ipvt, integer *nq, integer *nwait, real *rh, 
	real *rmax, real *save2, real *tq, real *trend, integer *iswflg, 
	integer *jstate)
{
    /* System generated locals */
    integer a_dim1, a_offset, yh_dim1, yh_offset, i__1;
    real r__1, r__2;

    /* Local variables */
    static integer i__;
    static real sum;
    static integer info;
    static real oldl0;
    extern doublereal snrm2_(integer *, real *, integer *);
    static integer iflag;
    extern /* Subroutine */ int sgbfa_(real *, integer *, integer *, integer *
	    , integer *, integer *, integer *), sgefa_(real *, integer *, 
	    integer *, integer *, integer *), sdscl_(real *, integer *, 
	    integer *, real *, real *, real *, real *, real *), sgbsl_(real *,
	     integer *, integer *, integer *, integer *, integer *, real *, 
	    integer *), sgesl_(real *, integer *, integer *, integer *, real *
	    , integer *), sdcst_(integer *, integer *, integer *, real *, 
	    real *);

/* ***BEGIN PROLOGUE  SDNTL */
/* ***SUBSIDIARY */
/* ***PURPOSE  Subroutine SDNTL is called to set parameters on the first */
/*            call to SDSTP, on an internal restart, or when the user has */
/*            altered MINT, MITER, and/or H. */
/* ***LIBRARY   SLATEC (SDRIVE) */
/* ***TYPE      SINGLE PRECISION (SDNTL-S, DDNTL-D, CDNTL-C) */
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
/*  If the caller has changed MINT, or if JTASK = 0, SDCST is called */
/*  to set the coefficients of the method.  If the caller has changed H, */
/*  YH must be rescaled.  If H or MINT has been changed, NWAIT is */
/*  reset to NQ + 2 to prevent further increases in H for that many */
/*  steps.  Also, RC is reset.  RC is the ratio of new to old values of */
/*  the coefficient L(0)*H.  If the caller has changed MITER, RC is */
/*  set to 0 to force the partials to be updated, if partials are used. */
/* ***ROUTINES CALLED  SDCST, SDSCL, SGBFA, SGBSL, SGEFA, SGESL, SNRM2 */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790601  DATE WRITTEN */
/*   900329  Initial submission to SLATEC. */
/* ***END PROLOGUE  SDNTL */
/* ***FIRST EXECUTABLE STATEMENT  SDNTL */
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
	    sdcst_(maxord, mint, iswflg, &el[14], &tq[4]);
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
		    sgefa_(&a[a_offset], matdim, n, &ipvt[1], &info);
		    if (info != 0) {
			*ier = TRUE_;
			return 0;
		    }
		    sgesl_(&a[a_offset], matdim, n, &ipvt[1], &save2[1], &
			    c__0);
		} else if (*miter == 4 || *miter == 5) {
		    (*fa)(n, t, &y[1], &a[*ml + 1 + a_dim1], matdim, ml, mu, 
			    nde);
		    if (*n == 0) {
			*jstate = 9;
			return 0;
		    }
		    sgbfa_(&a[a_offset], matdim, n, ml, mu, &ipvt[1], &info);
		    if (info != 0) {
			*ier = TRUE_;
			return 0;
		    }
		    sgbsl_(&a[a_offset], matdim, n, ml, mu, &ipvt[1], &save2[
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
		    if (a[i__ + a_dim1] == 0.f) {
			*ier = TRUE_;
			return 0;
		    } else {
			save2[i__] /= a[i__ + a_dim1];
		    }
/* L150: */
		}
		i__1 = *n;
		for (i__ = *nde + 1; i__ <= i__1; ++i__) {
/* L155: */
		    a[i__ + a_dim1] = 0.f;
		}
	    } else if (*impl == 3) {
		if (*miter == 1 || *miter == 2) {
		    (*fa)(n, t, &y[1], &a[a_offset], matdim, ml, mu, nde);
		    if (*n == 0) {
			*jstate = 9;
			return 0;
		    }
		    sgefa_(&a[a_offset], matdim, nde, &ipvt[1], &info);
		    if (info != 0) {
			*ier = TRUE_;
			return 0;
		    }
		    sgesl_(&a[a_offset], matdim, nde, &ipvt[1], &save2[1], &
			    c__0);
		} else if (*miter == 4 || *miter == 5) {
		    (*fa)(n, t, &y[1], &a[*ml + 1 + a_dim1], matdim, ml, mu, 
			    nde);
		    if (*n == 0) {
			*jstate = 9;
			return 0;
		    }
		    sgbfa_(&a[a_offset], matdim, nde, ml, mu, &ipvt[1], &info)
			    ;
		    if (info != 0) {
			*ier = TRUE_;
			return 0;
		    }
		    sgbsl_(&a[a_offset], matdim, nde, ml, mu, &ipvt[1], &
			    save2[1], &c__0);
		}
	    }
	}
	i__1 = *nde;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L170: */
/* Computing MAX */
	    r__1 = 1.f, r__2 = ywt[i__];
	    save1[i__] = save2[i__] / dmax(r__1,r__2);
	}
	sum = snrm2_(nde, &save1[1], &c__1) / sqrt((real) (*nde));
	if (sum > *eps / dabs(*h__)) {
	    r__1 = *eps / sum;
	    *h__ = r_sign(&r__1, h__);
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L180: */
	    yh[i__ + (yh_dim1 << 1)] = *h__ * save2[i__];
	}
	if (*miter == 2 || *miter == 5 || *iswflg == 3) {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
		fac[i__] = sqrt(*uround);
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
	    sdcst_(maxord, mint, iswflg, &el[14], &tq[4]);
	    *rc = *rc * el[*nq * 13 + 1] / oldl0;
	    *nwait = *nq + 2;
	}
	if (*h__ != *hold) {
	    *nwait = *nq + 2;
	    *rh = *h__ / *hold;
	    sdscl_(hmax, n, nq, rmax, hold, rc, rh, &yh[yh_offset]);
	}
    }
    return 0;
} /* sdntl_ */

